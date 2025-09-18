%% run_classify_new_exams.m
% INFERENCE ONLY: načtení natrénovaných modelů a klasifikace nových vyšetření

clc; clear;

% --- konfigurace & cesty ---
paths = jsondecode(fileread('./config/paths.local.json'));
cfg   = jsondecode(fileread('./config/model_config.json'));

bundle_file = './models/_globals/trained_bundle.mat';   % kde je uložen M, info, cal, nameMap
new_csv     = './data/to_test/test_pet_regions_with_mmi_JOINED_08092025.csv';  % << uprav dle sebe
reports_root= fullfile('reports','new_exams_from_csv');
if ~exist(reports_root,'dir'), mkdir(reports_root); end

% --- načti balíček ---
S = load(bundle_file, 'M','info','cal','nameMap');
M = S.M; info = S.info; cal = S.cal; nameMap = S.nameMap;

% --- (volitelně) načti trénink pro empirické boxploty ---
Ttrain = [];
try
    if isfield(paths,'raw_table') && exist(paths.raw_table,'file')
        T0 = azvpet.io.load_table(paths.raw_table);
        % jen kvůli empirickým boxplotům potřebujeme mít stejné featury jako při tréninku:
        covOpts = struct('lbmVersion','James', ...
                         'ageKnotsFile', cfg.age_splines.knots_file, ...
                         'refs', string(cfg.global_ref.refs), ...
                         'doseVar','InjectedDose_MBq','doseMultiplier',1);
        [T0, ~] = azvpet.features.ensure_model_covariates(T0, covOpts);
        baseNames = azvpet.features.list_region_bases(T0, 'Prefix', string(cfg.region_prefix));
        doseOpts  = struct('doseVar','InjectedDose_MBq','multiplier',1,'lbmVersion','James');
        [T0, ~]   = azvpet.features.ensureSUL_LOG(T0, baseNames, doseOpts);
        if isfield(cfg,'global_ref') && isfield(cfg.global_ref,'metric_suffix') ...
           && strcmpi(cfg.global_ref.metric_suffix,'_SUL_LOG')
            for r = string(cfg.global_ref.refs).'
                [T0, ~] = azvpet.features.ensureSUL_LOG(T0, char(r), doseOpts);
            end
        end
        % sjednoť validní názvy jako při tréninku
        [Ttrain, ~] = azvpet.util.ensure_valid_varnames(T0, struct('name_map',nameMap));
    end
catch ME
    warning('Nepodařilo se načíst/sjednotit trénink pro empirické boxploty: %s', ME.message);
    Ttrain = [];
end

% --- načti nové vyšetření/í z CSV ---
Tall = readtable(new_csv, 'VariableNamingRule','preserve');
nNew = height(Tall);
if nNew==0, error('V %s nejsou žádné řádky.', new_csv); end

% detekuj sloupec s ID (UNIS / PatientID)
idCol = intersect(string(Tall.Properties.VariableNames), ["UNIS","PatientID","ID"]);
if isempty(idCol), idCol = "row"; Tall.row = (1:nNew)'; end
idCol = char(idCol(1));

% helper pro přípravu nové tabulky do stejného formátu jako při tréninku
prep = @(Trow) prepare_new_exam_local(Trow, Ttrain, cfg, struct('name_map',nameMap));

% --- smyčka přes nová vyšetření ---
for i = 1:nNew
    Trow = Tall(i,:);
    % připrav
    try
        Tnew = prep(Trow);
    catch ME
        warning('Řádek %d přeskočen (příprava selhala): %s', i, ME.message);
        continue;
    end

    % slož jméno výstupního adresáře
    if ismember(idCol, string(Tall.Properties.VariableNames))
        pid = string(Tall.(idCol)(i));
        if ~isfinite(pid), pid = "unknown"; end
    else
        pid = "unknown";
    end
    outdir = fullfile(reports_root, sprintf('exam_%s_%d', char(pid), i));
    if ~exist(outdir,'dir'), mkdir(outdir); end

    % report (obě boxplot varianty se uloží uvnitř funkce)
    try
        if isempty(Ttrain)
            % nemáme empirické rozdělení – report to zvládne, uloží prediktivní box
            Srep = report_new_patient_all_regions(Tnew, Tnew, info, M, cal, outdir);
        else
            Srep = report_new_patient_all_regions(Ttrain, Tnew, info, M, cal, outdir);
        end

        % ukázka do konzole
        disp(Srep.summary_table(1:min(12,height(Srep.summary_table)), ...
            {'Region','Observed','Pred','PI_lo','PI_hi','z','p','p_outside','isOut95','isOut99'}));
        fprintf('[OK] Report pro %s -> %s\n', char(pid), outdir);
    catch ME
        warning('Report pro %s selhal: %s', char(pid), ME.message);
    end
end

%% ---------- lokální helper: příprava nového vyšetření ----------
function TnewC = prepare_new_exam_local(TnewRaw, Tref, cfg, optsC)
    % Tref: může být prázdné; používáme jen pro kategorie. Když není,
    % držíme se aspoň mapy názvů a dopočtu featur.
    Tnew = TnewRaw;

    % 1) dopočti stejné kovariáty
    covOpts = struct('lbmVersion','James', ...
                     'ageKnotsFile', cfg.age_splines.knots_file, ...
                     'refs', string(cfg.global_ref.refs), ...
                     'doseVar','InjectedDose_MBq','doseMultiplier',1);
    [Tnew, ~] = azvpet.features.ensure_model_covariates(Tnew, covOpts);

    % 2) stejné regiony + metrika
    if ~isempty(Tref)
        baseNames = azvpet.features.list_region_bases(Tref, 'Prefix', string(cfg.region_prefix));
    else
        % fallback: vezmi z názvů sloupců nové tabulky (pref. Median_* bez sufixu)
        baseNames = azvpet.features.list_region_bases(Tnew, 'Prefix', string(cfg.region_prefix));
    end
    doseOpts  = struct('doseVar','InjectedDose_MBq','multiplier',1,'lbmVersion','James');
    [Tnew, ~] = azvpet.features.ensureSUL_LOG(Tnew, baseNames, doseOpts);

    if isfield(cfg,'global_ref') && isfield(cfg.global_ref,'metric_suffix') ...
       && strcmpi(cfg.global_ref.metric_suffix,'_SUL_LOG')
        for r = string(cfg.global_ref.refs).'
            [Tnew, ~] = azvpet.features.ensureSUL_LOG(Tnew, char(r), doseOpts);
        end
    end

    % 3) sjednoť kategorie s tréninkem (pokud máme referenci)
    if ~isempty(Tref)
        vnames = string(Tref.Properties.VariableNames);
        for j = 1:numel(vnames)
            vn = vnames(j);
            if ismember(vn, string(Tnew.Properties.VariableNames)) && iscategorical(Tref.(vn))
                baseCats = categories(Tref.(vn));
                if ~iscategorical(Tnew.(vn)), Tnew.(vn) = categorical(string(Tnew.(vn))); end
                rawVals = string(Tnew.(vn));
    
                % rozšiř kategorie o nové hodnoty z Tnew (včetně nového UNIS)
                addCats = setdiff(unique(rawVals), baseCats);
                allCats = [baseCats; addCats];
    
                Tnew.(vn) = setcats(categorical(rawVals), allCats);
    
                if ~isempty(addCats)
                    warning('New exam: var %s – přidány nové kategorie: %s', vn, strjoin(cellstr(addCats'),', '));
                end
            end
        end
    end


    % 4) sanitace názvů přes mapu z tréninku
    if isfield(optsC,'name_map') && ~isempty(optsC.name_map)
        K = keys(optsC.name_map);   % <-- správně pro containers.Map
        for k = 1:numel(K)
            orig = K{k};
            safe = optsC.name_map(orig);
            if ismember(orig, Tnew.Properties.VariableNames) && ~ismember(safe, Tnew.Properties.VariableNames)
                Tnew.(safe) = Tnew.(orig);
            end
        end
    end

    TnewC = Tnew;
end
