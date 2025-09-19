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
    if ~isfield(opts,'muAge'),        opts.muAge = []; end
    if ~isfield(opts,'AgeR_knots'),   opts.AgeR_knots = []; end  
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
  
  % %% 1) Věk + spliny i polynomy
  % % Age v letech (datetime → double), fallback NaN
  % if ismember('Study Date',vn0) && ismember('PatientBirthDate',vn0)
  %     T.Age = years(T.("Study Date") - T.PatientBirthDate);
  % elseif ~ismember('Age',vn0)
  %     T.Age = nan(height(T),1);
  % end
  % T.Age = double(T.Age);
  % muAge = mean(T.Age,'omitnan');
  % T.cAge  = T.Age - muAge;
  % T.cAge2 = T.cAge.^2;
  % T.cAge3 = T.cAge.^3;
  % 
  % % Age RCS (fixní uzly – když jsou k dispozici)
  % if ~isempty(opts.ageKnotsFile) && isfile(opts.ageKnotsFile)
  %     S = load(opts.ageKnotsFile);   % očekává proměnnou 'knots'
  %     if isfield(S,'knots') && numel(S.knots)>=4
  %         [T, rcsi] = addAgeRCS_local(T,'Age','AgeR',numel(S.knots),S.knots);
  %         meta.AgeR_knots = rcsi.knots;
  %     else
  %         [T, rcsi] = addAgeRCS_local(T,'Age','AgeR',5,[]);
  %         meta.AgeR_knots = rcsi.knots;
  %     end
  % else
  %     [T, rcsi] = addAgeRCS_local(T,'Age','AgeR',5,[]);
  %     meta.AgeR_knots = rcsi.knots;
  % end
%% 1) Věk + spliny i polynomy
  % Najdi sloupce s datem vyšetření a narození včetně aliasů
  vn = string(T.Properties.VariableNames);

  studyAliases = ["Study Date","StudyDate","Study_Date","StudyDateTime","SeriesDate","Series_Date"];
  birthAliases = ["PatientBirthDate","Patient_Birth_Date","PatientBirthDateTime","BirthDate","Birth_Date"];

  hasAge = ismember("Age", vn);
  studyCol = studyAliases(ismember(studyAliases, vn));
  birthCol = birthAliases(ismember(birthAliases, vn));

  % Pokud Age chybí NEBO je celý NaN, zkus ho spočítat
  needComputeAge = ~hasAge || all(~isfinite(double(T.("Age"))));

  if needComputeAge && ~isempty(studyCol) && ~isempty(birthCol)
      scol = studyCol(1); bcol = birthCol(1);
      % převod na datetime (zkus text/číslo)
      Sdt = toDatetimeFlexible(T.(scol));
      Bdt = toDatetimeFlexible(T.(bcol));
      AgeY = years(Sdt - Bdt);
      T.Age = double(AgeY);
  elseif ~hasAge
      % fallback: někdy už existuje „AgeYears" nebo podobné
      cand = ["AgeYears","Age_Years","AgeY"];
      hit = cand(ismember(cand, vn));
      if ~isempty(hit)
          T.Age = double(T.(hit(1)));
      else
          T.Age = nan(height(T),1);
      end
  else
      % Age existuje – jen se ujisti, že je double
      T.Age = double(T.Age);
  end

  % if isfield(opts,'muAge') && isfinite(opts.muAge)
  %   muAge = opts.muAge;          % použij tréninkový průměr
  %   else
  %       muAge = mean(T.Age,'omitnan'); % fallback, když se netrénuje
  %   end
  %   T.cAge  = T.Age - muAge;
  %   T.cAge2 = T.cAge.^2;
  %   T.cAge3 = T.cAge.^3;
  %   meta.muAge = muAge;
  %     % Age RCS (fixní uzly – když jsou k dispozici)
  %     if ~isempty(opts.ageKnotsFile) && isfile(opts.ageKnotsFile)
  %         S = load(opts.ageKnotsFile);   % očekává proměnnou 'knots'
  %         if isfield(S,'knots') && numel(S.knots)>=4
  %             [T, rcsi] = addAgeRCS_local(T,'Age','AgeR',numel(S.knots),S.knots);
  %             meta.AgeR_knots = rcsi.knots;
  %         else
  %             [T, rcsi] = addAgeRCS_local(T,'Age','AgeR',5,[]);
  %             meta.AgeR_knots = rcsi.knots;
  %         end
  %     else
  %         [T, rcsi] = addAgeRCS_local(T,'Age','AgeR',5,[]);
  %         meta.AgeR_knots = rcsi.knots;
  %     end
  % použij tréninkový průměr, pokud je k dispozici
    if ~isempty(opts.muAge) && isfinite(opts.muAge)
        muAge = opts.muAge;
    else
        muAge = mean(T.Age,'omitnan');   % fallback (při tréninku)
    end
    T.cAge  = T.Age - muAge;
    T.cAge2 = T.cAge.^2;
    T.cAge3 = T.cAge.^3;
    meta.muAge = muAge;
    
    % --- Age spliny (RCS) – preferuj předané uzly ---
    if ~isempty(opts.AgeR_knots)
        knots = sanitize_knots_local(opts.AgeR_knots);
        [T, rcsi] = addAgeRCS_local(T,'Age','AgeR',numel(knots),knots);
        meta.AgeR_knots = rcsi.knots;
    elseif ~isempty(opts.ageKnotsFile) && isfile(opts.ageKnotsFile)
        S = load(opts.ageKnotsFile);
        if isfield(S,'knots') && numel(S.knots)>=4
            knots = sanitize_knots_local(S.knots);
            [T, rcsi] = addAgeRCS_local(T,'Age','AgeR',numel(knots),knots);
            meta.AgeR_knots = rcsi.knots;
        else
            [T, rcsi] = addAgeRCS_local(T,'Age','AgeR',5,[]);
            meta.AgeR_knots = rcsi.knots;
        end
    else
        % jen při tréninku; pro single-row inference je to nevhodné
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
% try
%     % 9a) Nejdřív zajisti *_SUL_LOG pro referenční regiony
%     refs = string(opts.refs(:));
%     refs = refs(refs~="");
%     if ~isempty(refs)
%         doseOpts = struct('doseVar',opts.doseVar, ...
%                           'multiplier',opts.doseMultiplier, ...
%                           'lbmVersion',opts.lbmVersion);
%         for r = refs.'
%             % r je base name (např. "Mean_Cerebellum_Gray_Matter_Right")
%             % ensureSUL_LOG vytvoří r+"_SUL_LOG", pokud zdrojové sloupce existují
%             [T, ~] = azvpet.features.ensureSUL_LOG(T, char(r), doseOpts);
%         end
%     end
% 
%     % 9b) Spočítej GlobalRefPC1_z pomocí uloženého projektoru (bez refitování)
%     if ~ismember('GlobalRefPC1_z', string(T.Properties.VariableNames))
%         pdir = './models/_globals';
%         % sjednoť název volby: preferuj opts.params_dir, ale toleruj i opts.paramsDir
%         if isfield(opts,'params_dir') && ~isempty(opts.params_dir)
%             pdir = opts.params_dir;
%         elseif isfield(opts,'paramsDir') && ~isempty(opts.paramsDir)
%             pdir = opts.paramsDir;
%         end
%         [T, ~] = azvpet.features.ensure_GlobalRefPC1_z(T, pdir);
%     end
% catch ME
%     % fallback – nic nerefituj, raději nech NaN (a dej vědět)
%     warning('ensure_model_covariates: GlobalRefPC1_z fallback (%s)', ME.message);
%     if ~ismember('GlobalRefPC1_z', string(T.Properties.VariableNames))
%         T.GlobalRefPC1_z = nan(height(T),1);
%     end
% end
%% 9) GlobalRefPC1_z
refs = string(opts.refs(:));
refs = refs(refs~="");

try
    % --- vždy nejdřív zajisti *_SUL_LOG pro referenční regiony
    if ~isempty(refs)
        doseOpts = struct('doseVar',opts.doseVar, ...
                          'multiplier',opts.doseMultiplier, ...
                          'lbmVersion',opts.lbmVersion);
        for r = refs.'
            [T, ~] = azvpet.features.ensureSUL_LOG(T, char(r), doseOpts);
        end
    end

    % --- kam ukládat / odkud číst projektor
    pdir = './models/_globals';
    if isfield(opts,'params_dir') && ~isempty(opts.params_dir)
        pdir = opts.params_dir;
    elseif isfield(opts,'paramsDir') && ~isempty(opts.paramsDir)
        pdir = opts.paramsDir;
    end

    if ~ismember('GlobalRefPC1_z', string(T.Properties.VariableNames))
        if isfield(opts,'mode') && strcmpi(opts.mode,'train')
            % >>> první běh: spočítej PCA a ulož projektor
            [T, ~] = azvpet.features.ensure_GlobalRefPC1_z(T, pdir, 'fit', refs);
        else
            % >>> predikce: použij uložený projektor
            [T, ~] = azvpet.features.ensure_GlobalRefPC1_z(T, pdir, 'predict', refs);
        end
    end
catch ME
    warning('ensure_model_covariates: GlobalRefPC1_z fallback (%s)', ME.message);
    if ~ismember('GlobalRefPC1_z', string(T.Properties.VariableNames))
        T.GlobalRefPC1_z = nan(height(T),1);
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
  
 % --- local helper uvnitř souboru ---
  function dt = toDatetimeFlexible(x)
      % Přijme string/char/numeric/datetime a vrátí datetime (NaT, když to nejde)
      if isdatetime(x)
          dt = x;
          return
      end
      if iscellstr(x) || isstring(x)
          try
              dt = datetime(x, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss', 'TimeZone','local');
          catch
              try
                  dt = datetime(x, 'InputFormat','yyyy-MM-dd', 'TimeZone','local');
              catch
                  try
                      dt = datetime(x); % last resort – nech to odhadnout
                  catch
                      dt = NaT(size(x));
                  end
              end
          end
          return
      end
      if isnumeric(x)
          % často bývá yyyymmdd (např. 20240131)
          try
              s = string(x);
              dt = datetime(s, 'InputFormat','yyyyMMdd', 'TimeZone','local');
          catch
              dt = NaT(size(x));
          end
          return
      end
      % jinak
      try
          dt = datetime(x);
      catch
          dt = NaT(size(x));
      end
  end

  function knots = sanitize_knots_local(knots)
    knots = knots(:).';
    knots = knots(isfinite(knots));
    knots = sort(unique(knots,'stable'));
    if numel(knots)<4
        error('AgeRCS: potřebuji ≥4 unikátní uzly, mám %d.', numel(knots));
    end
    % přísně rostoucí (ochrana proti numerické shodě posledních dvou)
    for i=2:numel(knots)
        if knots(i) <= knots(i-1)
            knots(i) = nextafter(knots(i-1), inf);
        end
    end
end
