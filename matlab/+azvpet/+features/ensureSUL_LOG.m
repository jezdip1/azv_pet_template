function [T, info] = ensureSUL_LOG(T, base, doseOpts)
%ENSURESUL_LOG Zajistí existenci sloupců *_SUL_LOG pro zadané regiony.
%   [T, info] = azvpet.features.ensureSUL_LOG(T, base, doseOpts)
%
%   base      : char | string | string array | cellstr (jeden nebo více názvů regionů)
%   doseOpts  : struct s volitelnými poli:
%               - doseVar: název sloupce s dávkou (např. 'InjectedDose_MBq' nebo 'Radiopharm_TotalDose_MBq')
%               - multiplier: škálování dávky (např. 1e-6, pokud je dávka v Bq)
%               - lbmVersion: 'James' (M/Ž) | 'Janmahasatian' | atd. (pokud není LBM sloupec hotový)
%
%   Funkce:
%     - pokud existuje <base>_SUL_LOG -> nedělá nic
%     - pokud existuje <base>_SUL     -> dopočítá log a uloží <base>_SUL_LOG
%     - jinak spočítá SUL z <base> při znalosti dávky a LBM:
%           SUL = conc_MBq_per_ml ./ (Dose_MBq ./ LBM_kg)
%       a uloží log(SUL) do <base>_SUL_LOG
%
%   Pozn.: Pro numerickou stabilitu používá eps při logu.
%
%   Vrací:
%     info: struct array se stavy pro každý region (‘created’ | ‘from_SUL’ | ‘exists’).

    if nargin < 3 || isempty(doseOpts), doseOpts = struct(); end

    % --- sjednocení vstupu 'base' na cellstr ---
    if ischar(base)
        bases = {base};
    elseif isstring(base)
        bases = cellstr(base);
    elseif iscellstr(base)
        bases = base;
    else
        error('ensureSUL_LOG:InvalidBaseType', ...
              'Parametr "base" musí být char/string/cellstr, dostal jsem %s.', class(base));
    end

    vn = T.Properties.VariableNames;
    info = struct('base', {}, 'status', {}, 'message', {});

    % --- připrav dávku a LBM (vektory přes řádky tabulky), pokud budou potřeba ---
    % Dávka (MBq)
    Dose_MBq = [];
    lbmKnown = ismember('LBM_kg', vn);
    if ~lbmKnown
        % Budeme případně dopočítávat LBM z hmotnosti/výšky, když není přímo LBM_kg
        needAnthro = any(ismember(bases, vn)); %#ok<NASGU> % jen placeholder, skutečný výpočet níže
    end

    % Příprava dávky lazily (jen pokud ji budeme potřebovat)
    function D = getDoseMBq()
        if ismember('InjectedDose_MBq', vn)
            D = T.InjectedDose_MBq; return;
        end
        if ~isempty(Dose_MBq), D = Dose_MBq; return; end
        % Preferuj explicitní volbu
        if isfield(doseOpts, 'doseVar') && ismember(doseOpts.doseVar, vn)
            D = T.(doseOpts.doseVar);
        else
            % Heuristika: zkus známé názvy
            cand = {'InjectedDose_MBq','Radiopharm_TotalDose_MBq','Dose_MBq','Total_Dose_MBq', ...
                    'InjectedDose_Bq','Radiopharm_TotalDose_Bq','Dose_Bq','Total_Dose_Bq'};
            hit = cand(ismember(cand, vn));
            if isempty(hit)
                error('ensureSUL_LOG:MissingDose', ...
                    'Nenašel jsem sloupec s dávkou (MBq/Bq). Zadej doseOpts.doseVar nebo přidej dávku do tabulky.');
            end
            D = T.(hit{1});
            % Pokud je v Bq, přeškáluj
            if endsWith(hit{1}, '_Bq')
                if isfield(doseOpts, 'multiplier') && ~isempty(doseOpts.multiplier)
                    D = D * doseOpts.multiplier;
                else
                    % default: Bq -> MBq
                    D = D * 1e-6;
                end
            end
        end
        Dose_MBq = D;
    end

    % LBM (kg)
    function L = getLBMkg()
        if ismember('LBM_kg', vn)
            L = T.('LBM_kg');
            return;
        end
        % Dopočet LBM z pohlaví, výšky, hmotnosti
        need = {'PatientWeight_kg','PatientHeight_cm','Sex'};
        have = ismember(need, vn);
        if ~all(have)
            missing = strjoin(need(~have), ', ');
            error('ensureSUL_LOG:MissingLBM', ...
                ['Chybí LBM_kg a zároveň chybí vstupy pro dopočet LBM (%s). ', ...
                 'Doplň LBM_kg nebo potřebná antropometrická data.'], missing);
        end
        wt = T.('PatientWeight_kg');      % kg
        ht = T.('PatientHeight_cm');      % cm
        sx = T.('Sex');                   % 'M'/'F' nebo logická/0-1
        if ~isstring(sx) && ~iscellstr(sx)
            % pokus o převod: 1 -> 'M', 0 -> 'F'
            sx = string(sx);
            sx = replace(sx, "1", "M");
            sx = replace(sx, "0", "F");
        else
            sx = string(sx);
        end

        % Volba vzorce
        version = 'James';
        if isfield(doseOpts, 'lbmVersion') && ~isempty(doseOpts.lbmVersion)
            version = string(doseOpts.lbmVersion);
        end

        % Výpočet
        L = nan(height(T),1);
        isM = (upper(sx)=="M");
        isF = (upper(sx)=="F");

        switch lower(version)
            case "janmahasatian"
                % Janmahasatian et al. Clin Pharmacokinet. 2005.
                % Muži:   LBM = 9270*wt / (6680 + 216*BMI)
                % Ženy:   LBM = 9270*wt / (8780 + 244*BMI)
                BMI = wt ./ ( (ht/100).^2 );
                L(isM) = 9270*wt(isM) ./ (6680 + 216*BMI(isM));
                L(isF) = 9270*wt(isF) ./ (8780 + 244*BMI(isF));
            otherwise % "James"
                % James (1976)
                % Muži: LBM = 1.10*wt - 128*(wt/ht)^2  [wt kg, ht cm]
                % Ženy: LBM = 1.07*wt - 148*(wt/ht)^2
                ratio = wt ./ ht; % kg/cm
                L(isM) = 1.10*wt(isM) - 128*(ratio(isM).^2);
                L(isF) = 1.07*wt(isF) - 148*(ratio(isF).^2);
        end
    end

    % Pomocné: bezpečný log
    function y = safelog(x)
        x = double(x);
        y = log(x + eps);
    end
    % Pomocné: převod koncentrace na MBq/ml s ohledem na Units po řádcích
    function concMBqml = toMBqml(conc)
        concMBqml = double(conc);
        if ismember('Units', vn)
            u = string(T.Units);
            scale = ones(height(T),1);
            scale(strcmpi(u,'BQML'))  = 1e-6;   % Bq/ml -> MBq/ml
            scale(strcmpi(u,'MBQML')) = 1;      % už je MBq/ml
            % ostatní jednotky (neznámé): ponecháme scale=1, ať to nespadne
            concMBqml = concMBqml .* scale;
        end
    end
 % --- hlavní smyčka přes regiony ---
    for k = 1:numel(bases)
        b = strtrim(bases{k});
        dstLog = [b '_SUL_LOG'];
        row.base = b; row.status = ''; row.message = '';

        if ismember(dstLog, vn)
            row.status = 'exists'; info(end+1) = row; %#ok<AGROW>
            continue;
        end

        srcSUL = [b '_SUL'];
        if ismember(srcSUL, vn)
            x = T.(srcSUL);
            y = nan(size(x));
            valid = isfinite(x) & x > 0;
            y(valid) = log(double(x(valid)));
            T.(dstLog) = y;
            vn = T.Properties.VariableNames;
            row.status = 'from_SUL';
            row.message = sprintf('valid=%d, NaN/else=%d', nnz(valid), numel(valid)-nnz(valid));
            info(end+1) = row; %#ok<AGROW>
            continue;
        end

        % potřebujeme zdrojový sloupec b
        if ~ismember(b, vn)
            row.status = 'error';
            row.message = sprintf('Chybí zdrojový sloupec "%s".', b);
            info(end+1) = row; %#ok<AGROW>
            error('ensureSUL_LOG:MissingSource', row.message);
        end

        % koncentrace -> MBq/ml (NaN-tolerant)
        conc = toMBqml(T.(b));

        % dávka a LBM: proměnné musí existovat, ale mohou obsahovat NaN
        D = []; LBM = [];
        % 1) dávka v MBq
        if ismember('InjectedDose_MBq', vn)
            D = T.InjectedDose_MBq;
        else
            D = getDoseMBq(); % tvoje původní lazy-funkce; jen neházej error na NaN
        end
        % 2) LBM v kg
        if ismember('LBM_kg', vn)
            LBM = T.LBM_kg;
        else
            LBM = getLBMkg(); % může vyhodit error jen pokud zcela chybí vstupní proměnné
        end

        % Předpočty + valid maska (NaN-tolerant)
        n = height(T);
        y = nan(n,1);

        % povolíme NaN v libovolné komponentě; vybereme jen řádky, kde jde výpočet provést
        valid = isfinite(conc) & isfinite(D) & isfinite(LBM) & (D > 0) & (LBM > 0);

        if any(valid)
            SUL = conc(valid) ./ (D(valid) ./ LBM(valid));    % bezrozměrné
            pos = isfinite(SUL) & (SUL > 0);
            tmp = nan(sum(valid),1);
            tmp(pos) = log(double(SUL(pos)));                 % log jen z kladných, jinak NaN
            y(valid) = tmp;
        end

        % Zápis výsledku (celý vektor včetně NaN)
        T.(dstLog) = y;
        vn = T.Properties.VariableNames;

        row.status  = 'created';
        row.message = sprintf('valid=%d, NaN/else=%d', nnz(valid), n - nnz(valid));
        info(end+1) = row; %#ok<AGROW>
    end
end
