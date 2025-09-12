function T = ensureSUL_LOG(T, base, opts)
    % ENSURE_SUL_LOG  Vytvoří <base>_SUL_LOG deterministicky.
    %
    %  - Pokud už <base>_SUL_LOG existuje, nedělá nic.
    %  - Jinak se pokusí spočítat SUL:
    %       1) Nejprve zkusí vytvořit SUVbw (pokud není přímo k dispozici).
    %          SUVbw = (C_tissue [kBq/mL]) * BW[kg] / Dose[kBq]
    %       2) SUL = SUVbw * (BW / LBM)  (tedy SUL = (C * BW^2) / (Dose * LBM))
    %       3) SUL_LOG = log(SUL)
    %
    %  - Jednotky & autodetekce:
    %       * C_tissue = T.(base) — očekává se kBq/mL; pokud je v Bq/mL, použij opts.activityMultiplier = 1e-3
    %       * Dose — pokusí se najít sloupec podle opts.doseVar nebo běžných názvů (viz detectDose)
    %         a zkonvertuje na kBq (MBq * 1000 nebo kBq přímo).
    %       * BW (kg) — pokusí se najít ve sloupcích: "Weight", "Weight_kg", "PatientWeight", "weight_kg"
    %       * Height (cm) — "Height", "Height_cm", "PatientSize_cm"; (pokud chybí, lze dopočítat z BMI)
    %       * BMI — "BMI", "bmi"
    %       * Sex — "sex", "Sex" (očekává 'M'/'F' nebo 1/0)
    %
    %  - LBM verze: opts.lbmVersion ∈ {'janmahasatian','james'}, default 'janmahasatian'
    %       Janmahasatian (2005):
    %         LBMm = 9270*W / (6680 + 216*BMI)
    %         LBMf = 9270*W / (8780 + 244*BMI)
    %       James:
    %         LBMm = 1.1*W - 128*(W/H)^2
    %         LBMf = 1.07*W - 148*(W/H)^2
    %
    %  - Fallback:
    %       * Pokud chybí klíčové veličiny a nelze spočítat SUL, zkusí:
    %           - má-li <base>_SUL → SUL_LOG = log(<base>_SUL)
    %           - jinak SUL_LOG = log(<base>)  (poslední fallback, aby pipeline běžela)
    %
    % Příklad:
    %   T = ensureSUL_LOG(T, "Mean_Caudate_Left", struct('lbmVersion','janmahasatian','activityMultiplier',1e-3));
    %
    % Pozn.: Volitelný vstup opts může být [], všechno se autodetekuje z tabulky.
    
    if nargin < 3, opts = struct(); end
    if isstring(base), base = char(base); end
    
    dst = [base '_SUL_LOG'];
    vn  = string(T.Properties.VariableNames);
    if ismember(dst, vn)
        return; % už existuje
    end
    assert(ismember(base, vn), 'ensureSUL_LOG: base column not found: %s', base);
    
    % ---------- 0) načtení a normalizace voleb ----------
    lbmVersion = getOpt(opts, 'lbmVersion', 'janmahasatian');       % 'janmahasatian' | 'james'
    actMult    = getOpt(opts, 'activityMultiplier', 1);             % 1 = už v kBq/mL; 1e-3 = Bq/mL -> kBq/mL
    doseVar    = getOpt(opts, 'doseVar', '');                       % název sloupce s dávkou (MBq/kBq), prázdné = autodetekce
    
    % ---------- 1) DATA: activity (C_tissue), dose (kBq), BW, H (cm), BMI, Sex ----------
    C = double(T.(base));             % očekáváno v kBq/mL; případně upravíme actMult
    C = C .* actMult;
    
    [dose_kBq, doseName] = detectDose(T, doseVar);   %#ok<ASGLU>
    [BW, Hcm, BMI, SEX]  = detectAnthro(T);
    
    % ---------- 2) LBM ----------
    LBM = computeLBM(BW, Hcm, BMI, SEX, lbmVersion);
    
    % ---------- 3) Pokus o výpočet SUL ----------
    canSUV = all(isfinite(C)) & all(isfinite(BW)) & all(isfinite(dose_kBq));
    canSUL = canSUV & all(isfinite(LBM));
    
    SUL = nan(height(T),1);
    if canSUL
        % SUVbw = C*kBq/mL * BW / dose_kBq  → SUL = SUVbw * BW/LBM
        SUVbw = (C .* BW) ./ dose_kBq;
        SUL   = SUVbw .* (BW ./ LBM);
        T.(dst) = log(SUL);
        return
    end
    
    % ---------- 4) fallbacky ----------
    % 4a) přímo SUL → SUL_LOG
    sulCol = base + "_SUL";
    if ismember(sulCol, vn)
        T.(dst) = log(double(T.(sulCol)));
        return
    end
    
    % 4b) přímo SUVbw → SUL z BW/LBM (pokud LBM jde spočítat)
    suvCol = base + "_SUVbw";
    if ismember(suvCol, vn) && all(isfinite(BW)) && all(isfinite(LBM))
        SUVbw = double(T.(suvCol));
        SUL   = SUVbw .* (BW ./ LBM);
        T.(dst) = log(SUL);
        return
    end
    
    % 4c) poslední možnost – log z raw intensity (explicitně, aby pipeline běžela)
    warning('ensureSUL_LOG: using fallback log(%s) due to missing dose/LBM inputs.', base);
    T.(dst) = log(double(T.(base)));
    end
    
    % ===== helpers ==============================================================
    function v = getOpt(s, name, default)
    if isfield(s, name), v = s.(name); else, v = default; end
    end
    
    function [dose_kBq, usedName] = detectDose(T, prefer)
    % Vrátí dávku v kBq; hledá typické názvy.
    names = string(T.Properties.VariableNames);
    cand = [string(prefer), ...
            "InjectedDose_MBq","Radiopharmaceutical_TotalDose_MBq","TotalDose_MBq", ...
            "InjectedDose_kBq","Radiopharmaceutical_TotalDose_kBq","TotalDose_kBq"];
    cand = cand(cand~="");
    usedName = "";
    dose_kBq = nan(height(T),1);
    
    for c = cand
        if ismember(c, names)
            x = double(T.(c));
            if endsWith(c,"_MBq"), dose_kBq = x * 1000; usedName = c; return; end
            if endsWith(c,"_kBq"), dose_kBq = x;         usedName = c; return; end
            % neznámé – tip: pokud čísla vypadají jako desítky/stovky → MBq
            if median(x,'omitnan') < 1000
                dose_kBq = x * 1000; usedName = c; return;
            else
                dose_kBq = x;         usedName = c; return;
            end
        end
    end
    
    % Pokud nic nenalezeno, zůstane NaN (zabere fallback)
    end
    
    function [BW, Hcm, BMI, SEX] = detectAnthro(T)
    names = string(T.Properties.VariableNames);
    
    % BW (kg)
    BW = pickFirst(T, names, ["Weight","Weight_kg","PatientWeight","weight_kg","BW_kg","BW"]);
    
    % Height (cm) – povoľ i metry s převodem
    Hcm = pickFirst(T, names, ["Height","Height_cm","PatientSize_cm","height_cm","Height_mm","PatientSize"]);
    if ~all(isfinite(Hcm))
        % Někdy bývá výška v metrech; zkus přepočet při podezřele malých číslech
        Hcm_guess = Hcm;
        mask = isfinite(Hcm_guess) & (Hcm_guess < 3.5); % zjevně metry
        Hcm_guess(mask) = Hcm_guess(mask) * 100;
        Hcm = Hcm_guess;
    end
    
    % BMI
    BMI = pickFirst(T, names, ["BMI","bmi","BodyMassIndex"]);
    
    % Sex (M/F or 1/0)
    SEXraw = pickFirst(T, names, ["sex","Sex","Gender","gender"]);
    SEX = nan(height(T),1);
    if iscellstr(SEXraw) || isstring(SEXraw) || iscategorical(SEXraw)
        s = string(SEXraw);
        SEX = double(upper(s)=="M");         % M→1, F→0
    else
        % číselné: předpoklad 1=male, 0=female (pokud jsou 1/2, mapni 2→0)
        SEX = double(SEXraw);
        if any(SEX==2), SEX(SEX==2) = 0; end
    end
    
    % Pokud BMI chybí, ale je BW & výška → BMI = W / (H[m])^2
    if ~all(isfinite(BMI)) && all(isfinite(BW)) && all(isfinite(Hcm))
        Hm = Hcm/100;
        BMI = BW ./ (Hm.^2);
    end
    end
    
    function v = pickFirst(T, names, candidates)
    v = nan(height(T),1);
    for c = string(candidates)
        if ismember(c, names)
            v = double(T.(c));
            return
        end
    end
    end
    
    function LBM = computeLBM(W, Hcm, BMI, SEX, version)
    % SEX: 1=male, 0=female
    W = double(W); Hcm = double(Hcm); BMI = double(BMI); SEX = double(SEX);
    LBM = nan(size(W));
    
    switch lower(string(version))
        case "james"
            % James (H v cm)
            male   = (SEX==1);
            female = (SEX==0);
            % ochrany proti dělení nulou
            Hcm(Hcm<=0 | ~isfinite(Hcm)) = NaN;
            LBM(male)   = 1.1.*W(male)   - 128.*(W(male)./Hcm(male)).^2;
            LBM(female) = 1.07.*W(female)- 148.*(W(female)./Hcm(female)).^2;
    
        otherwise % 'janmahasatian' (default)
            % Janmahasatian přes BMI
            male   = (SEX==1);
            female = (SEX==0);
            LBM(male)   = 9270.*W(male)   ./ (6680 + 216.*BMI(male));
            LBM(female) = 9270.*W(female) ./ (8780 + 244.*BMI(female));
    end
    
    % rozumné minimum (abychom nešli do záporných/0 hodnot)
    LBM(LBM<=0 | ~isfinite(LBM)) = NaN;
    end
    