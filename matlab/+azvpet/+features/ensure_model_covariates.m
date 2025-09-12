function [T, meta] = ensure_model_covariates(T, opts)
  % Zajistí kovariáty pro model typu:
  %   DV_SUL_LOG ~ AgeR* + Sex + BMI + lTime + lDose + Sex:cAge + lDose:lTime + ...
  %                logVoxelVol + logAcqDur_s + HasTOF + HasPSF + Subsets + FilterFWHM_mm +
  %                MMI_to_MNIGMS + NCC_to_MNIGMS + GlobalRefPC1_z + (1|UNIS)
  %
  % Použití:
  %   [T, meta] = azvpet.features.ensure_model_covariates(T, struct('ageKnotsFile','./age_knots.mat', 'lbmVersion','James', 'refs', ["Median_Cerebellum_Gray_Matter_Right","Median_Cerebellum_Gray_Matter_Left","Median_Brainstem_Right","Median_Brainstem_Left"]));
  %
  % Závislosti:
  %   - azvpet.features.ensure_metadata_vars
  %   - azvpet.features.ensureSUL_LOG
  %   - azvpet.features.ensure_GlobalRefPC1_z (volitelné; jinak fallback uvnitř)
  
  if nargin<2, opts=struct(); end
  if ~isfield(opts,'lbmVersion'),    opts.lbmVersion = 'James'; end
  if ~isfield(opts,'ageKnotsFile'),  opts.ageKnotsFile = '';    end
  if ~isfield(opts,'refs'),          opts.refs = strings(0,1);  end   % regiony pro GlobalRefPC1_z
  if ~isfield(opts,'doseVar'),       opts.doseVar = '';         end
  if ~isfield(opts,'doseMultiplier'),opts.doseMultiplier = [];  end   % Bq->MBq = 1e-6 default
  
  vn0 = string(T.Properties.VariableNames);
  meta = struct();
  
  %% 0) Metadata (dose/height/LBM, přejmenování Sex, atd.)
  [T, infoMeta] = azvpet.features.ensure_metadata_vars(T, struct('lbmVersion',opts.lbmVersion));
  meta.infoMeta = infoMeta;
  
  % sjednoť Sex na categorical{'M','F'} pokud možno
  if ismember('Sex', T.Properties.VariableNames)
      if ~iscategorical(T.Sex), T.Sex = categorical(string(T.Sex)); end
      % omez na M/F, ostatní nech jako missing
      s = upper(string(T.Sex));
      s(~ismember(s,["M","F"])) = missing;
      T.Sex = categorical(s, ["M","F"]);
  elseif ismember('PatientSex', T.Properties.VariableNames)
      s = string(T.PatientSex); s = upper(s);
      s = replace(s,"MALE","M"); s = replace(s,"FEMALE","F");
      s(~ismember(s,["M","F"])) = missing;
      T.Sex = categorical(s,["M","F"]);
  else
      T.Sex = categorical(repmat(missing,height(T),1),["M","F"]);
  end
  
  %% 1) Věk + spliny i polynomy
  % Age v letech (datetime → double), fallback NaN
  if ismember('Study Date',vn0) && ismember('PatientBirthDate',vn0)
      T.Age = years(T.("Study Date") - T.PatientBirthDate);
  elseif ~ismember('Age',vn0)
      T.Age = nan(height(T),1);
  end
  T.Age = double(T.Age);
  muAge = mean(T.Age,'omitnan');
  T.cAge  = T.Age - muAge;
  T.cAge2 = T.cAge.^2;
  T.cAge3 = T.cAge.^3;
  
  % Age RCS (fixní uzly – když jsou k dispozici)
  if ~isempty(opts.ageKnotsFile) && isfile(opts.ageKnotsFile)
      S = load(opts.ageKnotsFile);   % očekává proměnnou 'knots'
      if isfield(S,'knots') && numel(S.knots)>=4
          [T, rcsi] = addAgeRCS_local(T,'Age','AgeR',numel(S.knots),S.knots);
          meta.AgeR_knots = rcsi.knots;
      else
          [T, rcsi] = addAgeRCS_local(T,'Age','AgeR',5,[]);
          meta.AgeR_knots = rcsi.knots;
      end
  else
      [T, rcsi] = addAgeRCS_local(T,'Age','AgeR',5,[]);
      meta.AgeR_knots = rcsi.knots;
  end
  
  %% 2) BMI (kg/m^2)
  if ismember('PatientWeight_kg', string(T.Properties.VariableNames)) && ...
     ismember('PatientSize_m',    string(T.Properties.VariableNames))
      T.BMI = double(T.PatientWeight_kg) ./ (double(T.PatientSize_m).^2);
  elseif ~ismember('BMI', string(T.Properties.VariableNames))
      T.BMI = nan(height(T),1);
  end
  
  %% 3) lTime (log uptake time v minutách)
  if ismember('UptakeTime_min', string(T.Properties.VariableNames))
      T.lTime = safelog_local(T.UptakeTime_min);
  elseif ~ismember('lTime', string(T.Properties.VariableNames))
      T.lTime = nan(height(T),1);
  end
  
  %% 4) lDose (log dávky – vezmi InjectedDose_MBq pokud je, jinak RadionuclideTotalDose)
  if ismember('InjectedDose_MBq', string(T.Properties.VariableNames))
      T.lDose = safelog_local(T.InjectedDose_MBq);
  elseif ismember('RadionuclideTotalDose', string(T.Properties.VariableNames))
      % RadionuclideTotalDose je v Bq → log(Bq); pro konzistenci můžeš raději převést na MBq
      T.lDose = safelog_local(T.RadionuclideTotalDose);
  else
      T.lDose = nan(height(T),1);
  end
  
  %% 5) logVoxelVol
  % Přednostně VoxelVolume_mm3; jinak ze spacingu (mm)
  candVol = "VoxelVolume_mm3";
  if ismember(candVol, string(T.Properties.VariableNames))
      T.logVoxelVol = safelog_local(double(T.(candVol)));
  else
      % pokus: SpacingX_mm, SpacingY_mm, SpacingZ_mm → mm^3
      have = ismember(["SpacingX_mm","SpacingY_mm","SpacingZ_mm"], string(T.Properties.VariableNames));
      if all(have)
          vv = double(T.SpacingX_mm).*double(T.SpacingY_mm).*double(T.SpacingZ_mm);
          T.logVoxelVol = safelog_local(vv);
      elseif ~ismember('logVoxelVol', string(T.Properties.VariableNames))
          T.logVoxelVol = nan(height(T),1);
      end
  end
  
  %% 6) logAcqDur_s
  if ismember('AcquisitionDuration_s', string(T.Properties.VariableNames))
      T.logAcqDur_s = safelog_local(double(T.AcquisitionDuration_s));
  elseif ismember('SeriesDuration_s', string(T.Properties.VariableNames))
      T.logAcqDur_s = safelog_local(double(T.SeriesDuration_s));
  elseif ~ismember('logAcqDur_s', string(T.Properties.VariableNames))
      T.logAcqDur_s = nan(height(T),1);
  end
  
  %% 7) HasTOF / HasPSF / Subsets / FilterFWHM_mm (best-effort aliasing)
  % binární flagem: převeď na {0,1} i ze stringu
  T.HasTOF = ensure_flag_local(T, ["HasTOF","TOF"]);
  T.HasPSF = ensure_flag_local(T, ["HasPSF","PSF","PSFEnabled"]);
  % Subsets
  if ~ismember('Subsets', string(T.Properties.VariableNames))
      T.Subsets = pick_first_numeric_local(T, ["Subsets","OSEM_Subsets","RecoSubsets"]);
  end
  % FilterFWHM_mm
  if ~ismember('FilterFWHM_mm', string(T.Properties.VariableNames))
      T.FilterFWHM_mm = pick_first_numeric_local(T, ["FilterFWHM_mm","GaussFWHM_mm","PostFilterFWHM_mm"]);
  end
  
  %% 8) Registrace/QA posuny (ponech jak jsou; když chybí → NaN)
  for nm = ["MMI_to_MNIGMS","NCC_to_MNIGMS"]
      if ~ismember(nm, string(T.Properties.VariableNames))
          T.(nm) = nan(height(T),1);
      end
  end
  
  %% 9) GlobalRefPC1_z
  % Preferuj existující implementaci:
  try
      % if ~ismember('GlobalRefPC1_z', string(T.Properties.VariableNames))
      %     [T, ~] = azvpet.features.ensure_GlobalRefPC1_z(T, './models/_globals'); % přizpůsob cestu dle projektu
      % end
      if ~ismember('GlobalRefPC1_z', string(T.Properties.VariableNames))
        % použij configem daný adresář (pokud existuje)
        pdir = './models/_globals';
        if isfield(opts,'paramsDir') && ~isempty(opts.paramsDir)
            pdir = opts.paramsDir;
        end
        [T, ~] = azvpet.features.ensure_GlobalRefPC1_z(T, pdir);
      end
  catch
      % Fallback: z daných referenčních regionů (SUL_LOG), pokud jsou; jinak NaN
      refs = string(opts.refs(:));
      refs = refs(refs~="");
      if ~isempty(refs)
          % zajisti SUL_LOG pro refy
          for r = refs.'
              [T, ~] = azvpet.features.ensureSUL_LOG(T, char(r), struct('doseVar',opts.doseVar,'multiplier',opts.doseMultiplier,'lbmVersion',opts.lbmVersion));
          end
          refLog = refs + "_SUL_LOG";
          refLog = refLog(ismember(refLog, string(T.Properties.VariableNames)));
          if numel(refLog) >= 2
              X = [];
              for c = refLog(:).'
                  x = double(T.(char(c))); x(~isfinite(x)) = NaN;
                  mu = mean(x,'omitnan'); sd = std(x,'omitnan');
                  if isfinite(sd) && sd>0, X = [X, (x-mu)./sd]; end %#ok<AGROW>
              end
              T.GlobalRefPC1_z = nan(height(T),1);
              if size(X,2) >= 2
                  ok = all(isfinite(X),2);
                  if nnz(ok) >= 20
                      [~,score] = pca(X(ok,:), 'Centered', false, 'NumComponents',1);
                      pc1 = nan(height(T),1); pc1(ok) = score(:,1);
                      muPC = mean(pc1(ok),'omitnan'); sdPC = std(pc1(ok),0,'omitnan');
                      T.GlobalRefPC1_z = (pc1 - muPC)./sdPC;
                  end
              end
          else
              if ~ismember('GlobalRefPC1_z', string(T.Properties.VariableNames))
                  T.GlobalRefPC1_z = nan(height(T),1);
              end
          end
      else
          if ~ismember('GlobalRefPC1_z', string(T.Properties.VariableNames))
              T.GlobalRefPC1_z = nan(height(T),1);
          end
      end
  end
  
  %% 10) UNIS (random effect ID) – aliasing a typ
  if ismember('UNIS', string(T.Properties.VariableNames))
      if ~iscategorical(T.UNIS), T.UNIS = categorical(string(T.UNIS)); end
  else
      id = pick_first_text_local(T, ["UNIS","SiteID","Prefix","PatientID"]);
      if isempty(id), id = repmat("UNK", height(T),1); end
      T.UNIS = categorical(id);
  end
  
  %% 11) Interakce – nic se nepřipravuje, fitlme je sestaví z termů
  % Jen jistota, že cAge existuje (už jsme vytvořili)
  % lDose:lTime vytvoří model, my jen držíme lDose,lTime
  
  meta.muAge = muAge;
  meta.AgeR_names = string(T.Properties.VariableNames(startsWith(string(T.Properties.VariableNames),'AgeR')));
  end % === main ===
  
  
  %% --- lokální pomůcky ---
  function y = safelog_local(x)
      x = double(x);
      y = log(max(eps, x));
  end
  
  function [T, info] = addAgeRCS_local(T, ageVar, prefix, nKnots, knots)
      if ~ismember(ageVar, string(T.Properties.VariableNames))
          error('addAgeRCS_local: missing %s', ageVar);
      end
      x = double(T.(ageVar));
      x(~isfinite(x)) = NaN;
      if nargin<5 || isempty(knots)
          p = [5 27.5 50 72.5 95];
          knots = prctile(x(isfinite(x)), p);
          knots = unique(knots(:)');
          if numel(knots) < max(4, nKnots)
              kmin = min(x); kmax = max(x);
              knots = linspace(kmin, kmax, nKnots);
          end
      else
          knots = sort(knots(:)');
      end
      K = numel(knots);
      if K < 4, error('addAgeRCS_local: potřebuji ≥4 uzly, mám %d.', K); end
  
      B = rcsDesign_local(x, knots);  % [N × (K-1)]
      for j = 1:size(B,2)
          T.(sprintf('%s%d', prefix, j)) = B(:,j);
      end
      info = struct('knots', knots, 'K', K);
  end
  
  function B = rcsDesign_local(x, knots)
      x = x(:);
      K = numel(knots); kK = knots(K); kKm = knots(K-1);
      d = @(u,c) ((u - c).*(u > c)).^3;
      Z = zeros(numel(x), K-2);
      for j = 1:(K-2)
          kj = knots(j);
          Z(:,j) = d(x,kj) - d(x,kK)*(kK-kj)/(kK-kKm) + d(x,kKm)*(kKm-kj)/(kK-kKm);
      end
      B = [x, Z]; % 1. sloupec = lineární věk
  end
  
  function v = pick_first_numeric_local(T, names)
      v = nan(height(T),1);
      for nm = string(names)
          if ismember(nm, string(T.Properties.VariableNames))
              x = double(T.(nm));
              if any(isfinite(x)), v = x; return; end
          end
      end
  end
  
  function v = pick_first_text_local(T, names)
      v = strings(height(T),1);
      for nm = string(names)
          if ismember(nm, string(T.Properties.VariableNames))
              x = string(T.(nm));
              if any(~ismissing(x)), v = x; return; end
          end
      end
  end
  
  function z = ensure_flag_local(T, cands)
      z = nan(height(T),1);
      for nm = string(cands)
          if ismember(nm, string(T.Properties.VariableNames))
              col = T.(nm);
              if isnumeric(col) || islogical(col)
                  x = double(col);
                  if any(~isnan(x)), z = double(x~=0); return; end
              else
                  s = lower(string(col));
                  tf = ismember(s, ["1","true","t","yes","y","on"]);
                  ff = ismember(s, ["0","false","f","no","n","off"]);
                  out = nan(height(T),1);
                  out(tf) = 1; out(ff) = 0;
                  if any(~isnan(out)), z = out; return; end
              end
          end
      end
      if all(isnan(z)), z = zeros(height(T),1); end  % default 0 (bez TOF/PSF)
  end
  