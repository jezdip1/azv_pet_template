% End-to-end: load → features → LMEM → LOO → calibration → report → artifacts
azvpet.util.check_requirements();
paths = jsondecode(fileread('./config/paths.local.json'));
featCfg = jsondecode(fileread('./config/features_config.json'));
cfg  = jsondecode(fileread('./config/model_config.json'));
T0 = azvpet.io.load_table(paths.raw_table);
% cfg = load_model_config('model_config.json');

% 0) metadata + kovariáty (předávej věci z configu)
covOpts = struct( ...
    'lbmVersion',     'James', ...
    'ageKnotsFile',   cfg.age_splines.knots_file, ...
    'refs',           string(cfg.global_ref.refs), ...
    'doseVar',        'InjectedDose_MBq', ...   % máme z ensure_metadata_vars
    'doseMultiplier', 1 ...
);
[T0, metaCov] = azvpet.features.ensure_model_covariates(T0, covOpts);

% 1) seznam regionů podle configu
baseNames = azvpet.features.list_region_bases(T0, 'Prefix', string(cfg.region_prefix));

% 2) zajisti *_SUL_LOG pro všechny regiony (NaN-tolerantně)
doseOpts = struct('doseVar','InjectedDose_MBq','multiplier',1,'lbmVersion','James');
[T0, ~] = azvpet.features.ensureSUL_LOG(T0, baseNames, doseOpts);

% 3) pokud global_ref.metric_suffix je '_SUL_LOG', zajisti totéž i pro refy:
if isfield(cfg,'global_ref') && isfield(cfg.global_ref,'metric_suffix') ...
   && strcmpi(cfg.global_ref.metric_suffix,'_SUL_LOG')
    for r = string(cfg.global_ref.refs).'
        [T0, ~] = azvpet.features.ensureSUL_LOG(T0, char(r), doseOpts);
    end
end

% 4) volitelně: předej params_dir implementaci, která čte uložené PCA parametry
% (pokud ji máš – já v ensure_model_covariates už volám ensure_GlobalRefPC1_z s './models/_globals')

[form, opts] = azvpet.model.define_formula(cfg, T0);

[Tc, nameMap, optsC] = azvpet.util.ensure_valid_varnames(T0, opts);

[M, info]    = azvpet.model.train_lmem(Tc, form, opts);


% cv = azvpet.model.loo_cv(Tc, form, opts); %% pozor delal jsme jeste
% zasahy do te paralelni verze, takze bude potreba overit ta neparalelni

cv = azvpet.model.loo_cv_par(Tc, form, optsC);

cal = azvpet.model.calibrate_from_cv(M, cv);
% cal = azvpet.model.calibrate_model(Tc, M, info, cv);
azvpet.model.report_model(Tc, M, info, cv, cal, paths);
azvpet.io.save_artifacts(M, info, cal, paths, featCfg, mdlCfg);
