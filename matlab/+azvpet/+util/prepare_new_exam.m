function [TnewC] = prepare_new_exam(TnewRaw, Ttrain, cfg, optsC)
% Zarovná nový řádek na trénink:
% - dopočítá stejné kovariáty (splines, glob.ref, lDose, lTime, …)
% - zajistí *_SUL_LOG pro stejné regiony
% - sjednotí kategorie s Ttrain
% - sjednotí názvy sloupců přes nameMap (sanitize)

    arguments
        TnewRaw table
        Ttrain  table
        cfg     struct
        optsC   struct
    end

    Tnew = TnewRaw;

    % --- 1) stejné kovariáty jako v tréninku ---
    covOpts = struct( ...
        'lbmVersion',     'James', ...
        'ageKnotsFile',   cfg.age_splines.knots_file, ...
        'refs',           string(cfg.global_ref.refs), ...
        'doseVar',        'InjectedDose_MBq', ...
        'doseMultiplier', 1 ...
    );
    [Tnew, ~] = azvpet.features.ensure_model_covariates(Tnew, covOpts);

    % --- 2) stejné regiony a metriky (SUL_LOG) ---
    baseNames = azvpet.features.list_region_bases(Ttrain, 'Prefix', string(cfg.region_prefix));
    doseOpts  = struct('doseVar','InjectedDose_MBq','multiplier',1,'lbmVersion','James');
    [Tnew, ~] = azvpet.features.ensureSUL_LOG(Tnew, baseNames, doseOpts);

    if isfield(cfg,'global_ref') && isfield(cfg.global_ref,'metric_suffix') ...
       && strcmpi(cfg.global_ref.metric_suffix,'_SUL_LOG')
        for r = string(cfg.global_ref.refs).'
            [Tnew, ~] = azvpet.features.ensureSUL_LOG(Tnew, char(r), doseOpts);
        end
    end

    % --- 3) sjednotit kategorie přesně na levels z tréninku ---
    vnames = string(Ttrain.Properties.VariableNames);
    for j = 1:numel(vnames)
        vn = vnames(j);
        if ismember(vn, string(Tnew.Properties.VariableNames))
            if iscategorical(Ttrain.(vn))
                cats = categories(Ttrain.(vn));         % levels z tréninku
                % pokud Tnew není categorical, udělej ho
                if ~iscategorical(Tnew.(vn))
                    Tnew.(vn) = categorical(string(Tnew.(vn)));
                end
                % přemapuj mimo-slovník na <missing>, a fixní pořadí
                Tnew.(vn) = setcats(categorical(Tnew.(vn)), cats);
                % vše mimo cats spadne na <missing> -> predikce pak selže,
                % takže radši warn:
                bad = ~ismember(string(TnewRaw.(vn)), cats);
                if any(bad)
                    warning('New exam: var %s has unseen categories: %s', vn, strjoin(unique(string(TnewRaw.(vn)(bad)))',', '));
                end
            end
        end
    end

    % --- 4) sanitace názvů stejně jako u tréninku (optsC.name_map) ---
    if isfield(optsC,'name_map') && ~isempty(optsC.name_map)
        % name_map je containers.Map(orig -> safe)
        keys = optsC.name_map.keys;
        for k = 1:numel(keys)
            orig = keys{k}; safe = optsC.name_map(orig);
            if ismember(orig, Tnew.Properties.VariableNames) && ~ismember(safe, Tnew.Properties.VariableNames)
                Tnew.(safe) = Tnew.(orig);
            end
        end
    end

    % --- 5) finále: ponech jen sloupce, které model umí (volitelně) ---
    % (většinou není nutné; hlavní je, aby *vše potřebné* existovalo a nebylo NaN.)

    TnewC = Tnew;
end
