%% run_train_and_classify_new_exam.m
% End-to-end: TRAIN -> LOO -> CALIB -> REPORT -> NEW EXAM CLASSIFICATION

clc; clear;

% --- 0) sanity / cesty / config ---
azvpet.util.check_requirements();

paths  = jsondecode(fileread('./config/paths.local.json'));
cfg    = jsondecode(fileread('./config/model_config.json'));
featCfg= jsondecode(fileread('./config/features_config.json')); %#ok<NASGU>

reports_root = fullfile('reports');
if ~exist(reports_root,'dir'), mkdir(reports_root); end

% --- 1) načtení tréninkové tabulky (RAW) ---
T0 = azvpet.io.load_table(paths.raw_table);

% --- 2) kovariáty přesně dle configu (Age splines, lDose, lTime, PCA ref, ...) ---
covOpts = struct( ...
    'lbmVersion',     'James', ...
    'ageKnotsFile',   cfg.age_splines.knots_file, ...
    'refs',           string(cfg.global_ref.refs), ...
    'doseVar',        'InjectedDose_MBq', ...
    'doseMultiplier', 1 ...
);
[T0, ~] = azvpet.features.ensure_model_covariates(T0, covOpts);

% --- 3) seznam regionů a výpočet *_SUL_LOG (NaN-tolerantně) ---
baseNames = azvpet.features.list_region_bases(T0, 'Prefix', string(cfg.region_prefix));
doseOpts  = struct('doseVar','InjectedDose_MBq','multiplier',1,'lbmVersion','James');
[T0, ~]   = azvpet.features.ensureSUL_LOG(T0, baseNames, doseOpts);

% pokud je global ref ve stejné metrice, taky dopočti
if isfield(cfg,'global_ref') && isfield(cfg.global_ref,'metric_suffix') ...
   && strcmpi(cfg.global_ref.metric_suffix,'_SUL_LOG')
    for r = string(cfg.global_ref.refs).'
        [T0, ~] = azvpet.features.ensureSUL_LOG(T0, char(r), doseOpts);
    end
end

% --- 4) modelová formule + sanitace názvů tak, jak je bude chtít fitlme ---
[form, opts] = azvpet.model.define_formula(cfg, T0);
[Tc, nameMap, optsC] = azvpet.util.ensure_valid_varnames(T0, opts);
optsC.name_map = nameMap;   % důležité pro mapování názvů i mimo trénink

% --- 5) TRÉNINK LMEM ---
[M, info] = azvpet.model.train_lmem(Tc, form, optsC);

% --- 6) Leave-One-Group-Out CV (paralelně) ---
cv = azvpet.model.loo_cv_par(Tc, form, optsC);

% --- 7) KALIBRACE z LOO (alpha, beta, c, smear) ---
cal = azvpet.model.calibrate_from_cv_full(M, cv);

% --- 8) ULOŽENÍ balíčku pro inference ---
if ~exist('./models/_globals','dir'), mkdir('./models/_globals'); end
save('./models/_globals/trained_bundle.mat', 'M','info','cv','cal','nameMap','-v7.3');

% --- 9) REPORT z tréninku (SUL originální škála) ---
rep_opts = struct('y_scale',"SUL", 'sanitize',true, ...
                  'add_partial_resid',true, 'add_loo',true);
azvpet.model.report_model(Tc, M, info, cv, cal, paths, rep_opts);

fprintf('[OK] Training + CV + Calibration hotovo. Report uložen do %s\n', fullfile(paths.reports_dir,'latest'));

%% ========== NOVÝ PACIENT Z CSV ==========

% --- 10) Načti CSV s novým vyšetřením (musí mít stejné surové sloupce jako T0) ---
new_csv = '/mnt/data/test_pet_regions_with_mmi_JOINED_08092025.csv';  % <- uprav dle sebe
TnewRaw = readtable(new_csv);

% --- 11) Zarovnej nový řádek na trénink (doplní kovariáty, SUL_LOG, sjednotí kategorie, názvy) ---
if exist('+azvpet/+util/prepare_new_exam.m','file')
    Tnew = azvpet.util.prepare_new_exam(TnewRaw(1,:), Tc, cfg, optsC);
else
    % fallback: lokální helper definovaný níže v tomto skriptu
    Tnew = prepare_new_exam_local(TnewRaw(1,:), Tc, cfg, optsC);
end

% --- 12) Klasifikační report přes všechny regiony (obě varianty boxplotu) ---
outdir = fullfile(reports_root, 'new_exam_from_csv');
S = report_new_patient_all_regions(Tc, Tnew, info, M, cal, outdir);

% --- 13) rychlý náhled do výsledků ---
disp(S.summary_table(1:min(15,height(S.summary_table)), ...
    {'Region','Observed','Pred','PI_lo','PI_hi','z','p','p_outside','isOut95','isOut99'}));

fprintf('[OK] New-exam report hotový: %s\n', outdir);

%% --------------- lokální helper, pokud není v +azvpet/+util ----------------
function TnewC = prepare_new_exam_local(TnewRaw, Ttrain, cfg, optsC)
    Tnew = TnewRaw;

    % stejné kovariáty
    covOpts = struct( ...
        'lbmVersion',     'James', ...
        'ageKnotsFile',   cfg.age_splines.knots_file, ...
        'refs',           string(cfg.global_ref.refs), ...
        'doseVar',        'InjectedDose_MBq', ...
        'doseMultiplier', 1 ...
    );
    [Tnew, ~] = azvpet.features.ensure_model_covariates(Tnew, covOpts);

    % stejné regiony a metriky
    baseNames = azvpet.features.list_region_bases(Ttrain, 'Prefix', string(cfg.region_prefix));
    doseOpts  = struct('doseVar','InjectedDose_MBq','multiplier',1,'lbmVersion','James');
    [Tnew, ~] = azvpet.features.ensureSUL_LOG(Tnew, baseNames, doseOpts);

    if isfield(cfg,'global_ref') && isfield(cfg.global_ref,'metric_suffix') ...
       && strcmpi(cfg.global_ref.metric_suffix,'_SUL_LOG')
        for r = string(cfg.global_ref.refs).'
            [Tnew, ~] = azvpet.features.ensureSUL_LOG(Tnew, char(r), doseOpts);
        end
    end

    % sjednocení kategorií (levels) s tréninkem
    vnames = string(Ttrain.Properties.VariableNames);
    for j = 1:numel(vnames)
        vn = vnames(j);
        if ismember(vn, string(Tnew.Properties.VariableNames)) && iscategorical(Ttrain.(vn))
            cats = categories(Ttrain.(vn));
            if ~iscategorical(Tnew.(vn)), Tnew.(vn) = categorical(string(Tnew.(vn))); end
            rawVals = string(Tnew.(vn));
            Tnew.(vn) = setcats(categorical(Tnew.(vn)), cats);
            bad = ~ismember(rawVals, cats);
            if any(bad)
                warning('New exam: var %s má neznámé kategorie: %s', vn, strjoin(unique(rawVals(bad))',', '));
            end
        end
    end

    % sanitace názvů přes nameMap
    if isfield(optsC,'name_map') && ~isempty(optsC.name_map)
        keys = optsC.name_map.keys;
        for k = 1:numel(keys)
            orig = keys{k}; safe = optsC.name_map(orig);
            if ismember(orig, Tnew.Properties.VariableNames) && ~ismember(safe, Tnew.Properties.VariableNames)
                Tnew.(safe) = Tnew.(orig);
            end
        end
    end

    TnewC = Tnew;
end
