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
%%
load('/home/jezdip1/ISARGUBU/media/shared_storage/motol/AZV_PET/tools/azv_pet_template/matlab/test_PET_08092025.mat')

[T0, metaCov] = azvpet.features.ensure_model_covariates(T0, covOpts);
[Tnew, metaCov_new] = azvpet.features.ensure_model_covariates(Tnew, covOpts);

% 1) seznam regionů podle configu
baseNames = azvpet.features.list_region_bases(T0, 'Prefix', string(cfg.region_prefix));
baseNames_new = azvpet.features.list_region_bases(Tnew, 'Prefix', string(cfg.region_prefix));

% 2) zajisti *_SUL_LOG pro všechny regiony (NaN-tolerantně)
doseOpts = struct('doseVar','InjectedDose_MBq','multiplier',1,'lbmVersion','James');
[T0, ~] = azvpet.features.ensureSUL_LOG(T0, baseNames, doseOpts);
[Tnew, ~] = azvpet.features.ensureSUL_LOG(Tnew, baseNames_new, doseOpts);

% 3) pokud global_ref.metric_suffix je '_SUL_LOG', zajisti totéž i pro refy:
if isfield(cfg,'global_ref') && isfield(cfg.global_ref,'metric_suffix') ...
   && strcmpi(cfg.global_ref.metric_suffix,'_SUL_LOG')
    for r = string(cfg.global_ref.refs).'
        [T0, ~] = azvpet.features.ensureSUL_LOG(T0, char(r), doseOpts);
        [Tnew, ~] = azvpet.features.ensureSUL_LOG(Tnew, char(r), doseOpts);
    end
end

% 4) volitelně: předej params_dir implementaci, která čte uložené PCA parametry
% (pokud ji máš – já v ensure_model_covariates už volám ensure_GlobalRefPC1_z s './models/_globals')

[form, opts] = azvpet.model.define_formula(cfg, T0);
[form_new, opts_new] = azvpet.model.define_formula(cfg, Tnew);





[Tc, nameMap, optsC] = azvpet.util.ensure_valid_varnames(T0, opts);
[Tnew, nameMap_new, optsC_new] = azvpet.util.ensure_valid_varnames(Tnew, opts_new);

%%

one = Tnew(1, :);

outdir = fullfile('reports','new_exam_SOME_ID');
S = report_new_patient_all_regions(Tc, one, info, M, cal, outdir);

disp(S.summary_table(1:10, {'Region','Observed','Pred','PI_lo','PI_hi','z','isOut95','isOut99'}))
