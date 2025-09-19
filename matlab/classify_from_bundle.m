function S = classify_from_bundle(new_csv, out_root)
% Klasifikace nových vyšetření bez retrainingu – čistě s uloženým bundlem.

azvpet.util.check_requirements();

% --- cesty & config
paths = jsondecode(fileread('./config/paths.local.json'));
cfg   = jsondecode(fileread('./config/model_config.json'));

Bfile = './models/_globals/trained_bundle.mat';
if ~isfile(Bfile), error('Nenalezen %s. Nejdřív spusť train_pipeline.', Bfile); end
% B = load(Bfile, 'M','info','cal','nameMap');  % nameMap: ORIG<->CLEAN z tréninku
% M = B.M; info = B.info; cal = B.cal; nameMap = B.nameMap;

% --- načti nová data, zachovej ORIG hlavičky
Tnew0 = readtable(new_csv, 'VariableNamingRule','preserve');
% %%
% % % --- stejné kovariáty jako u tréninku (AGE splines, cAge, GlobalRefPC1_z, …)
% % %     ensure_model_covariates očekává, že PCA parametry najde v ./models/_globals
% % covOpts = struct( ...
% %   'lbmVersion',   'James', ...
% %   'ageKnotsFile', cfg.age_splines.knots_file, ...
% %   'refs',         string(cfg.global_ref.refs), ...
% %   'doseVar',      'InjectedDose_MBq', ...
% %   'doseMultiplier', 1, ...
% %   'paramsDir',   './models/_globals' ... % <<< bylo params_dir, musí být paramsDir
% % );
% % [Tnew0, ~] = azvpet.features.ensure_model_covariates(Tnew0, covOpts);
% B = load('./models/_globals/trained_bundle.mat','M','info','cal','nameMap','muAge');
% covOpts = struct( ...
%   'lbmVersion','James', ...
%   'ageKnotsFile', cfg.age_splines.knots_file, ...
%   'refs', string(cfg.global_ref.refs), ...
%   'doseVar','InjectedDose_MBq', ...
%   'doseMultiplier',1, ...
%   'paramsDir','./models/_globals', ... % pro PCA
%   'muAge', B.muAge ...                  % <<< důležité
% );
% [Tnew0, ~] = azvpet.features.ensure_model_covariates(Tnew0, covOpts);
B = load('./models/_globals/trained_bundle.mat','M','info','cal','nameMap','metaCov');
M = B.M; info = B.info; cal = B.cal; nameMap = B.nameMap;
muAge_tr   = B.metaCov.muAge;
AgeR_knots = B.metaCov.AgeR_knots;

% ... načti Tnew0 ...

covOpts = struct( ...
  'lbmVersion',   'James', ...
  'ageKnotsFile', '', ...                 % NEspoléhat na soubor
  'AgeR_knots',   AgeR_knots, ...         % <<< předat přímo
  'muAge',        muAge_tr, ...           % <<< stejné centrování jako v tréninku
  'refs',         string(cfg.global_ref.refs), ...
  'doseVar',      'InjectedDose_MBq', ...
  'doseMultiplier', 1, ...
  'params_dir',   './models/_globals', ... % pro GlobalRefPC1_z
  'mode','predict'...
);
[Tnew0, ~] = azvpet.features.ensure_model_covariates(Tnew0, covOpts);


head(Tnew0(:, {'Age','cAge','AgeR1','AgeR2','AgeR3','AgeR4','GlobalRefPC1_z','Sex','BMI','lTime','lDose'}))
%%
% Tnew0 = removevars(Tnew0, "GlobalRefPC1_z"); % trochu hloupej trik

% --- zajisti *_SUL_LOG pro VŠECHNY modelované regiony (na ORIG názvech!)
resp_clean   = string(info.responses(:));           % CLEAN: ..._SUL_LOG
bases_clean  = regexprep(resp_clean, '_SUL_LOG$','');
[tf,loc]     = ismember(bases_clean, string(nameMap.clean));
bases_orig   = bases_clean;
bases_orig(tf) = string(nameMap.orig(loc(tf)));     % ORIG: např. I-V místo I_V

doseOpts = struct('doseVar','InjectedDose_MBq','multiplier',1,'lbmVersion','James');
[Tnew0, ~] = azvpet.features.ensureSUL_LOG(Tnew0, cellstr(unique(bases_orig)), doseOpts);

% pokud referenční regiony v cfg běží také na _SUL_LOG, dopočítej je též
if isfield(cfg,'global_ref') && isfield(cfg.global_ref,'metric_suffix') ...
   && strcmpi(cfg.global_ref.metric_suffix,'_SUL_LOG')
  for r = string(cfg.global_ref.refs).'
    [Tnew0, ~] = azvpet.features.ensureSUL_LOG(Tnew0, char(r), doseOpts);
  end
end
%%
% % --- stejné kovariáty jako u tréninku (AGE splines, cAge, GlobalRefPC1_z, …)
% %     ensure_model_covariates očekává, že PCA parametry najde v ./models/_globals
% covOpts = struct( ...
%   'lbmVersion',   'James', ...
%   'ageKnotsFile', cfg.age_splines.knots_file, ...
%   'refs',         string(cfg.global_ref.refs), ...
%   'doseVar',      'InjectedDose_MBq', ...
%   'doseMultiplier', 1, ...
%   'paramsDir',   './models/_globals' ... % <<< bylo params_dir, musí být paramsDir
% );
% [Tnew0, ~] = azvpet.features.ensure_model_covariates(Tnew0, covOpts);
%%
% --- sjednoť názvy na CLEAN stejně jako při tréninku (responses + grouping)
opts = struct('grouping','UNIS','response_list',{cellstr(resp_clean)});
[Tnew, ~, opts] = azvpet.util.ensure_valid_varnames(Tnew0, opts);
opts.name_map = nameMap;  %#ok<NASGU> % (kdyby ho později chtěl někdo použít dál)

% % --- UNIS & kategorie: drobné dorovnání na trénink
% if ~ismember('UNIS', Tnew.Properties.VariableNames)
%   Tnew.UNIS = categorical(repmat("<undefined>",height(Tnew),1));
% end
% if ~iscategorical(Tnew.UNIS), Tnew.UNIS = categorical(string(Tnew.UNIS)); end

% % zarovnej hladiny kategorií podle prvního modelu (sex, scanner, …)
% km = keys(M);
% firstKey = km{1};
% Lany  = M(firstKey);
% Train = Lany.Variables;
% vnames = intersect(Tnew.Properties.VariableNames, Train.Properties.VariableNames);
% for j = 1:numel(vnames)
%   v = vnames{j};
%   if iscategorical(Train.(v))
%     levTrain = categories(Train.(v));
%     if ~iscategorical(Tnew.(v)), Tnew.(v) = categorical(string(Tnew.(v))); end
%     % přidej případné nové hladiny, ale na konci je seřaď jako v tréninku
%     Tnew.(v) = setcats(Tnew.(v), union(levTrain, categories(Tnew.(v))));
%     Tnew.(v) = categorical(Tnew.(v), levTrain);
%   end
% end
% 
% Tnew.UNIS=Tnew0.UNIS;

% zarovnej hladiny kategorií podle prvního modelu (sex, scanner, …)
km = keys(M);
firstKey = km{1};
Lany  = M(firstKey);
Train = Lany.Variables;

vnames = intersect(Tnew.Properties.VariableNames, Train.Properties.VariableNames);
for j = 1:numel(vnames)
  v = vnames{j};
  if iscategorical(Train.(v))
    levTrain = categories(Train.(v));
    if ~iscategorical(Tnew.(v))
        Tnew.(v) = categorical(string(Tnew.(v)));
    end
    % 1) povol všechny hladiny: tréninkové i nové
    allLevels = union(levTrain, categories(Tnew.(v)), 'stable');
    Tnew.(v) = setcats(Tnew.(v), allLevels);
    % 2) seřaď tak, aby tréninkové byly první, nové až za nimi
    newOnly   = setdiff(categories(Tnew.(v)), levTrain, 'stable');
    Tnew.(v)  = reordercats(Tnew.(v), [levTrain; newOnly]);
  end
end


% --- výstupní složka
if ~exist(out_root,'dir'), mkdir(out_root); end

% head(Tnew0(:, {'Age','cAge','GlobalRefPC1_z','Sex','BMI','lTime','lDose'}))
head(Tnew(:, {'AgeR1','AgeR2','AgeR3','AgeR4','BMI','FilterFWHM_mm','GlobalRefPC1_z','HasPSF','HasTOF','MMI_to_MNIGMS','NCC_to_MNIGMS','Sex','Subsets','UNIS','cAge','lTime','lDose','logAcqDur_s','logVoxelVol'}))
% --- klasifikuj po řádcích (MODEL-boxplot varianta; Ttrain nepotřebujeme)
S = struct('cases',[]);
for i = 1:height(Tnew)
  row = Tnew(i,:);
  % robustní ID: když UNIS chybí -> row_####
  sid = "<undefined>";
  if ismember('UNIS', row.Properties.VariableNames)
    sid = string(row.UNIS); sid = sid(1);
  end
  if ismissing(sid) || strlength(sid)==0, sid = sprintf("row_%04d", i); end
  outdir = fullfile(out_root, safe_fname(sid));
  if ~exist(outdir,'dir'), mkdir(outdir); end

  S_case = report_new_patient_all_regions([], row, info, M, cal, outdir);
  S.cases = [S.cases; struct('id',char(sid), 'summary',S_case.summary_table, ...
              'csv',S_case.csv,'json',S_case.json,'zplot',S_case.z_plot)];
end

fprintf('[OK] Classified %d new rows -> %s\n', height(Tnew), out_root);
end

% --- lokální util – aby to nespadlo na zvláštních znacích v UNIS
function s = safe_fname(txt)
s = regexprep(char(string(txt)),'[^A-Za-z0-9\-]+','_');
end
