function report_model(T, M, info, cv, cal, paths, opts)
% REPORT pro všechny regiony v M (klíče mapy).
% Ukládá PNG + JSON se souhrnem metrik.
%
% Požaduje:
%  - T ........ tabulka (s validními názvy proměnných)
%  - M ........ containers.Map (resp -> LinearMixedModel)
%  - info ..... (nepovinné, jen kvůli responses; lze ignorovat)
%  - cv ....... struct z loo_cv_par (pole resp se sloupci: Age, y_true, y_pred, CI_lo_link_mean, CI_hi_link_mean, sd_pred_link, resid)
%  - cal ...... struct z calibrate_from_cv_full (alpha,beta,c,smear,localSD_fun,...)
%  - paths .... struct s poli: reports_dir (adresář pro výstupy)
%
% Výstupy:
%  - ./reports/latest/<resp>_fixed_forest.png
%  - ./reports/latest/<resp>_age_effect.png
%  - ./reports/latest/<resp>_calibration.png
%  - ./reports/latest/<resp>_zscore_vs_age.png
%  - ./reports/latest/summary.json
    if nargin < 7 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'y_scale'),   opts.y_scale = "SUL"; end   % "SUL" | "RAW" | "LINK"
    if ~isfield(opts,'sanitize'),  opts.sanitize = true;  end
    out = fullfile(paths.reports_dir,'latest');
    if ~exist(out,'dir'), mkdir(out); end

    % pokus o načtení uzlů pro AgeR spliny (není povinné)
    ageKnots = [];
    try
        if exist('./models/_globals/age_knots.mat','file'), L = load('./models/_globals/age_knots.mat'); end
        if ~exist('L','var') && exist('./age_knots.mat','file'), L = load('./age_knots.mat'); end
        if exist('L','var') && isfield(L,'knots'), ageKnots = L.knots(:)'; end
    catch
        ageKnots = [];
    end

    % klíče = regiony, které jsme opravdu fitovali
    ks = keys(M);

    % souhrn (pro JSON)
    SUM = struct();

    for i = 1:numel(ks)
        resp = ks{i};
        mdl  = M(resp);

        % ---------- 0) metriky ----------
        try
            [R2m, R2c] = azvpet.util.stats_utils.r2_lmm(mdl);
        catch
            R2m = NaN; R2c = NaN;
        end
        AIC = mdl.ModelCriterion.AIC;
        BIC = mdl.ModelCriterion.BIC;

        SUM.(resp) = struct('AIC',AIC,'BIC',BIC,'R2m',R2m,'R2c',R2c);

        % ---------- 1) Fixed-effects forest ----------
        try
            C = coeftable_compat(mdl.Coefficients);  % robustní převod na table
            if ~isstring(C.Name), C.Name = string(C.Name); end
            ok = isfinite(C.Estimate) & isfinite(C.SE) & isfinite(C.DF);
            C = C(ok,:);
            C = C(~strcmpi(C.Name,'(Intercept)'),:);
            if ~isempty(C)
                tcrit = tinv(0.975, C.DF);
                lo = C.Estimate - tcrit.*C.SE;
                hi = C.Estimate + tcrit.*C.SE;

                f = figure('Visible','off','Color','w'); hold on
                y = 1:height(C);
                for j=1:numel(y)
                    plot([lo(j) hi(j)], [y(j) y(j)], '-'); hold on
                end
                plot(C.Estimate, y, 'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','k');
                xline(0,'k:');
                set(gca,'YTick',y,'YTickLabel',cellstr(strrep(C.Name,'_','-')),'YDir','reverse');
                set_xlabel(gca, 'Estimate (link)', opts.sanitize); set_title(gca, [nice_label(resp) ' — Fixed effects (95% CI)'], opts.sanitize);
                close(f);
            end
        catch ME
            warning('Forest plot (%s): %s', resp, ME.message);
        end

        % ---------- 2) Age effect (Sex strata) s CI+PI na orig škále ----------
        try
            if isfield(cal, resp)
                Cc = cal.(resp);
            else
                Cc = struct('alpha',0,'beta',1,'c',1,'smear',NaN,'localSD_fun',[]);
            end

            hasAge = ismember('Age', T.Properties.VariableNames);
            hasSex = ismember('Sex', T.Properties.VariableNames);
            if hasAge
                % trén. věky kvůli rozsahu
                ageTrain = double(T.Age); ageTrain = ageTrain(isfinite(ageTrain));
                aLo = prctile(ageTrain,5); aHi = prctile(ageTrain,95);
                age = linspace(aLo,aHi,100)';

                % referenční řádek: mediány/moda
                ref = refRowLike(mdl, T, resp);

                % sex křivky
                nCurves = 1; sexCats = string({});
                if hasSex
                    if ~iscategorical(T.Sex), T.Sex = categorical(string(T.Sex)); end
                    sexCats = string(categories(removecats(T.Sex)));
                    nCurves = numel(sexCats);
                end

                % připrav predikce
                mu = nan(numel(age), nCurves);
                loM= mu; hiM = mu;    % CI mean (link -> orig)
                loP= mu; hiP = mu;    % PI (link -> orig)

                % jak sestavit AgeR* / polynomy:
                feNames = string(mdl.Variables.Properties.VariableNames);
                useRCS  = any(startsWith(feNames,"AgeR"));
                usePoly = any(ismember(["cAge","cAge2","cAge3"], feNames));

                % back-transform
                isLogDV = endsWith(resp,'_LOG') || endsWith(resp,'_SUL_LOG');
                smear = Cc.smear; if ~isfinite(smear) || smear<=0, smear=1; end
                if isLogDV, tr = @(x) exp(x).*smear; else, tr = @(x) x; end

                for s = 1:nCurves
                    newT = repmat(ref, numel(age), 1);
                    if hasSex
                        newT.Sex = categorical(repmat(sexCats(s), numel(age), 1), categories(T.Sex));
                    end
                    newT.Age = age;

                    % --- doplň AgeR* nebo polynomy ---
                    if useRCS && ~isempty(ageKnots)
                        B = rcsDesign(age, ageKnots);  % 1. sloupec = Age
                        for j = 1:size(B,2)
                            nm = sprintf('AgeR%d', j);
                            if ismember(nm, newT.Properties.VariableNames)
                                newT.(nm) = B(:,j);
                            end
                        end
                    elseif usePoly
                        muAge = mean(ageTrain,'omitnan');
                        if ismember('cAge', newT.Properties.VariableNames),  newT.cAge  = age - muAge; end
                        if ismember('cAge2',newT.Properties.VariableNames), newT.cAge2 = (age - muAge).^2; end
                        if ismember('cAge3',newT.Properties.VariableNames), newT.cAge3 = (age - muAge).^3; end
                    else
                        % nic – nebude age efekt
                        continue;
                    end

                    % population-level predikce (link)
                    [yhat_l, yCI_l] = predict(mdl, newT, 'Conditional', false, 'Alpha',0.05);
                    z975 = norminv(0.975);
                    SEm  = (yCI_l(:,2)-yCI_l(:,1))/(2*z975);
                    s2res= mdl.MSE;
                    addV = addedREvariance(mdl, newT);
                    sd_l = sqrt(max(0, SEm.^2 + s2res + addV));

                    % kalibrace link (alpha + beta*y)
                    yhat_cal_l = Cc.alpha + Cc.beta.*yhat_l;

                    % CI mean (nekalibrujeme – je to inovace kolem mean)
                    mu(:,s)  = tr(yhat_cal_l);
                    loM(:,s) = tr(yCI_l(:,1));
                    hiM(:,s) = tr(yCI_l(:,2));

                    % PI (kalibrované c-faktorem)
                    loP(:,s) = tr(yhat_cal_l - z975*(Cc.c*sd_l));
                    hiP(:,s) = tr(yhat_cal_l + z975*(Cc.c*sd_l));
                end

                % kreslení
                f = figure('Visible','off','Color','w'); hold on
                cols = lines(nCurves);
                for s=1:nCurves
                    % PI pásmo
                    if all(isfinite(loP(:,s))) && all(isfinite(hiP(:,s)))
                        fill([age; flipud(age)], [loP(:,s); flipud(hiP(:,s))], cols(s,:), ...
                             'FaceAlpha',0.10,'EdgeColor','none'); hold on
                    end
                    % CI mean pásmo
                    if all(isfinite(loM(:,s))) && all(isfinite(hiM(:,s)))
                        fill([age; flipud(age)], [loM(:,s); flipud(hiM(:,s))], cols(s,:), ...
                             'FaceAlpha',0.20,'EdgeColor','none'); hold on
                    end
                    % střed
                    plot(age, mu(:,s), 'LineWidth',2, 'Color', cols(s,:));
                end
                % výběr škály a převodu
                scale = upper(string(opts.y_scale));   % "SUL" | "RAW" | "LINK"
                ylab = y_label_for(resp, scale);
                

                if scale == "RAW"
                    % odhadni mapování RAW ≈ a + b * SUL(orig)
                    map = raw_from_sul_mapping(T, resp);
                    if map.ok
                        mu  = map.a + map.b * mu;
                        loM = map.a + map.b * loM;   hiM = map.a + map.b * hiM;
                        loP = map.a + map.b * loP;   hiP = map.a + map.b * hiP;
                    else
                        warning('RAW scale requested, but RAW mapping not available -> using SUL.');
                        ylab = y_label_for(resp, "SUL");
                    end
                end

                
                % kreslení s čistými popisky
                f = figure('Visible','off','Color','w'); hold on
                cols = lines(nCurves);
                for s=1:nCurves
                    if all(isfinite(loP(:,s))) && all(isfinite(hiP(:,s)))
                        fill([age; flipud(age)], [loP(:,s); flipud(hiP(:,s))], cols(s,:), ...
                             'FaceAlpha',0.10,'EdgeColor','none');
                    end
                    if all(isfinite(loM(:,s))) && all(isfinite(hiM(:,s)))
                        fill([age; flipud(age)], [loM(:,s); flipud(hiM(:,s))], cols(s,:), ...
                             'FaceAlpha',0.20,'EdgeColor','none');
                    end
                    plot(age, mu(:,s), 'LineWidth',2, 'Color', cols(s,:));
                end
                set_ylabel(gca, ylab, opts.sanitize);
                set_xlabel(gca, 'Age [y]', opts.sanitize); grid on
                if hasSex
                    leg = {};
                    for s=1:nCurves
                        leg = [leg, {sprintf('PI | Sex=%s', string(sexCats(s)))}, ...
                                    {sprintf('CI | Sex=%s', string(sexCats(s)))}, ...
                                    {sprintf('Mean | Sex=%s', string(sexCats(s)))}];
                    end
                    set_legend(gca, {leg}, opts.sanitize);
                end
                set_title(gca, sprintf('%s — Age effect (%s, CI & PI)', nice_label(regexprep(resp,'_LOG$','')), upper(string(ylab(end-3:end-1)))), opts.sanitize);
                % saveas(f, fullfile(out,[resp '_age_effect.png'])); close(f);

                saveas(f, fullfile(out,[resp '_age_effect.png'])); close(f);
            end
        catch ME
            warning('Age effect (%s): %s', resp, ME.message);
        end

        % ---------- 3) Kalibrace (Obs vs Pred, orig scale) ----------
        try
            if isfield(cv, resp)
                TT = cv.(resp);
                % back-transform
                isLogDV = endsWith(resp,'_LOG') || endsWith(resp,'_SUL_LOG');
                smear = cal.(resp).smear; if ~isfinite(smear) || smear<=0, smear=1; end
                if isLogDV
                    Obs  = exp(double(TT.y_true))*smear;
                    Pred = exp(double(cal.(resp).alpha + cal.(resp).beta.*TT.y_pred))*smear;
                else
                    Obs  = double(TT.y_true);
                    Pred = double(cal.(resp).alpha + cal.(resp).beta.*TT.y_pred);
                end

                v = isfinite(Obs) & isfinite(Pred);
                if nnz(v) >= 5
                    % decile binning
                    Pc = Pred(v); Oc = Obs(v);
                    edges = prctile(Pc, 0:10:100);
                    edges(1)=min(Pc); edges(end)=max(Pc);
                    [~,~,bin] = histcounts(Pc, edges);
                    K = max(bin);
                    xC = nan(K,1); yC = nan(K,1);
                    for kbin=1:K
                        xC(kbin) = mean(Pc(bin==kbin),'omitnan');
                        yC(kbin) = mean(Oc(bin==kbin),'omitnan');
                    end

                    f = figure('Visible','off','Color','w'); hold on
                    plot(xC, yC, 'o-','LineWidth',1.5);
                    mm = [min(Pc) max(Pc)];
                    plot(mm, mm, 'k--'); % ideál
                    % kalibrovaná přímka na orig škále: aproximace přes link je nelineární u LOG,
                    % ale pro vizuál stačí spojit body (x, pred->orig) vs (x, obs->orig) z decilů.
                    title([resp ' — Calibration (10 bins)']); xlabel('Pred (orig)'); ylabel('Obs (orig)'); grid on
                    saveas(f, fullfile(out,[resp '_calibration.png'])); close(f);
                end
            end
        catch ME
            warning('Calibration plot (%s): %s', resp, ME.message);
        end

        % ---------- 4) Z-score (standardized residuals) vs Age ----------
        try
            if isfield(cv, resp)
                TT = cv.(resp);
                Cc = cal.(resp);
                y  = double(TT.y_true);
                p  = double(TT.y_pred);
                sd = double(TT.sd_pred_link);
                z  = (y - p) ./ (Cc.c * sd);
                A  = double(TT.Age);

                v = isfinite(z) & isfinite(A);
                f = figure('Visible','off','Color','w'); hold on
                plot(A(v), z(v), '.', 'MarkerSize',10);
                yline(0,'k:'); yline(2,'r--'); yline(3,'r-'); yline(-2,'r--'); yline(-3,'r-');
                xlabel('Age [y]'); set_ylabel(gca, 'z = (Obs - Pred) / (c * sd_pred)', opts.sanitize);
                set_title(gca, [nice_label(resp) ' — Standardized residuals'], opts.sanitize); grid on
                saveas(f, fullfile(out,[resp '_zscore_vs_age.png'])); close(f);
            end
        catch ME
            warning('Z-score plot (%s): %s', resp, ME.message);
        end
            % === Partial residual plot (Age) ===
        try
            y_true = double(T.(resp));
            y_fit  = fitted(mdl);  % celé predikce
        
            % najdi sloupce designMatrix odpovídající Age a Sex
            coefNames = mdl.CoefficientNames;
            keepIdx = contains(coefNames,'AgeR') | strcmp(coefNames,'Sex') | contains(coefNames,'Sex:Age');
            X = mdl.designMatrix;  % n x p
        
            beta = mdl.Coefficients.Estimate;
            y_age = X(:,keepIdx) * beta(keepIdx);
        
            % partial residuals
            y_part = y_true - (y_fit - y_age);
        
            % převed na požadovanou škálu
            if scale=="SUL"
                y_part = exp(y_part) * cal.(resp).smear;
            elseif scale=="RAW"
                map = raw_from_sul_mapping(T, resp);
                if map.ok
                    y_part = map.a + map.b * (exp(y_part)*cal.(resp).smear);
                end
            end
        
            f = figure('Visible','off','Color','w');
            hold on
            scatter(T.Age, y_part, 20, 'o', 'MarkerEdgeAlpha',0.3);
            plot(age, mu(:,1), 'r-', 'LineWidth',2); % model curve (už máš spočítanou výše)
            set_xlabel(gca, 'Age [y]', opts.sanitize);
            set_ylabel(gca, ylab, opts.sanitize);
            set_title(gca, sprintf('%s — Partial residuals (Age)', nice_label(resp)), opts.sanitize);
            saveas(f, fullfile(out,[resp '_partial_resid.png']));
            close(f);
        catch ME
            warning('Partial residual plot failed for %s: %s', resp, ME.message);
        end
        % === Predicted vs Observed (LOO) ===
        try
            cvt = cv.(resp);
            y_true = double(cvt.y_true);
            y_pred = double(cvt.y_pred);
        
            % převeď na požadovanou škálu
            if scale=="SUL"
                y_true = exp(y_true) * cal.(resp).smear;
                y_pred = exp(y_pred) * cal.(resp).smear;
            elseif scale=="RAW"
                map = raw_from_sul_mapping(T, resp);
                if map.ok
                    y_true = map.a + map.b * (exp(y_true)*cal.(resp).smear);
                    y_pred = map.a + map.b * (exp(y_pred)*cal.(resp).smear);
                end
            end
        
            % metriky
            R = corr(y_true, y_pred, 'rows','complete');
            R2 = R^2;
            MAE = mean(abs(y_true - y_pred),'omitnan');
            RMSE = sqrt(mean((y_true - y_pred).^2,'omitnan'));
        
            f = figure('Visible','off','Color','w'); hold on
            scatter(y_true, y_pred, 20, 'filled', 'MarkerFaceAlpha',0.3);
            lims = [min([y_true;y_pred]), max([y_true;y_pred])];
            plot(lims, lims, 'k--','LineWidth',1.5); % diagonála
        
            set_xlabel(gca, 'Observed', opts.sanitize);
            set_ylabel(gca, 'Predicted (LOO)', opts.sanitize);
            set_title(gca, sprintf('%s — Pred vs Obs (LOO)\nR²=%.3f, MAE=%.3f, RMSE=%.3f',...
                nice_label(resp), R2, MAE, RMSE), opts.sanitize);
            axis equal; xlim(lims); ylim(lims);
            saveas(f, fullfile(out,[resp '_pred_vs_obs.png']));
            close(f);
        catch ME
            warning('Pred vs Obs plot failed for %s: %s', resp, ME.message);
        end
    end
    % ---------- 5) summary.json ----------
    try
        fid = fopen(fullfile(out,'summary.json'),'w');
        fwrite(fid, jsonencode(SUM,'PrettyPrint',true));
        fclose(fid);
    catch
        warning('Nemohu zapsat summary.json');
    end
end

% ==== Helpers (lokální kopie, ať je funkce self-contained) ====
function TT = coeftable_compat(C)
    if istable(C), TT=C; return; end
    try
        TT = dataset2table(C,'ReadRowNames',false);
    catch
        try
            S = dataset2struct(C); TT = struct2table(S);
        catch
            TT = struct2table(orderfields(C));
        end
    end
    if ismember('Term', TT.Properties.VariableNames) && ~ismember('Name', TT.Properties.VariableNames)
        TT.Properties.VariableNames{strcmpi(TT.Properties.VariableNames,'Term')} = 'Name';
    end
    if ismember('Name', TT.Properties.VariableNames) && ~isstring(TT.Name)
        TT.Name = string(TT.Name);
    end
end

function ref = refRowLike(L, Tref, DVname)
    varNames = L.Variables.Properties.VariableNames;
    ref = Tref(1,:);
    for i=1:numel(varNames)
        nm = string(varNames{i});
        if strcmp(nm, DVname), continue; end
        if ~ismember(nm, string(Tref.Properties.VariableNames)), continue; end
        col = Tref.(nm);
        if isnumeric(col) || islogical(col)
            ref.(nm)(1) = median(double(col),'omitnan');
        elseif isduration(col)
            sec = seconds(col); ref.(nm)(1) = seconds(median(sec,'omitnan'));
        elseif isdatetime(col)
            ii = find(~ismissing(col),1,'first'); if ~isempty(ii), ref.(nm)(1)=col(ii); end
        elseif iscategorical(col)
            cc = removecats(col);
            if isempty(categories(cc))
                ref.(nm)(1) = categorical(missing, categories(col));
            else
                ref.(nm)(1) = mode(cc);
            end
        else
            ii = find(~ismissing(col),1,'first'); if ~isempty(ii), ref.(nm)(1)=col(ii); end
        end
    end
end

function B = rcsDesign(x, knots)
    x = x(:);
    K = numel(knots);
    kK  = knots(K); kKm = knots(K-1);
    d = @(u,c) ((u - c).*(u > c)).^3;
    Z = zeros(numel(x), K-2);
    for j=1:(K-2)
        kj = knots(j);
        Z(:,j) = d(x,kj) - d(x,kK)*(kK-kj)/(kK-kKm) + d(x,kKm)*(kKm-kj)/(kK-kKm);
    end
    B = [x, Z]; % AgeR1..AgeR(K-1)
end

function addVar = addedREvariance(L, TT)
    addVar = zeros(height(TT),1);
    try
        [~, Info] = covarianceParameters(L);
        if ~istable(Info), Info = dataset2table(Info,'ReadRowNames',false); end
        G  = string(Info.('Group'));
        Ty = string(Info.('Type'));
        Nm = string(Info.('Name1'));
        Est= double(Info.('Estimate'));
        mask = Ty=="std" & Nm=="(Intercept)";
        groups = unique(G(mask));
        Train = L.Variables;
        for g = groups'
            gn = char(g);
            if ~ismember(gn, TT.Properties.VariableNames) || ~ismember(gn, Train.Properties.VariableNames)
                continue;
            end
            s2g = Est(find(G==g & mask,1)).^2;
            isNew = ~ismember(TT.(gn), removecats(categorical(Train.(gn))));
            addVar = addVar + s2g .* double(isNew);
        end
    catch
    end
end

function s = nice_label(varname)
    % uživatelský titulek (bez _ a ^ jako indexů)
    s = strrep(varname, '_', '-');
    s = strrep(s, '^', ' ');
end

function set_title(h, txt, sanitize)
    if sanitize
        title(h, txt, 'Interpreter','none');
    else
        title(h, txt);
    end
end
function set_xlabel(h, txt, sanitize)
    if sanitize
        xlabel(h, txt, 'Interpreter','none');
    else
        xlabel(h, txt);
    end
end
function set_ylabel(h, txt, sanitize)
    if sanitize
        ylabel(h, txt, 'Interpreter','none');
    else
        ylabel(h, txt);
    end
end
function set_legend(h, args, sanitize)
    if sanitize
        legend(h, args{:}, 'Interpreter','none', 'Location','best');
    else
        legend(h, args{:}, 'Location','best');
    end
end

function ylab = y_label_for(resp, scale)
    base = regexprep(resp,'_SUL_LOG$','');      % bez suffixu
    switch upper(string(scale))
        case "SUL"
            ylab = [nice_label(base) ' (SUL)'];
        case "RAW"
            ylab = [nice_label(base) ' (RAW)'];
        case "LINK"
            ylab = [nice_label(resp) ' (link/log)'];
        otherwise
            ylab = nice_label(base);
    end
end

function f = raw_multiplier_from_T(Tin, resp)
% Vrátí mediánový násobek RAW/SUL pro daný resp (…_SUL_LOG).
% Když chybí _SUL, pokusí se ho spočítat přes ensureSUL_LOG.
    T = Tin;   % lokální kopie (nechceme měnit volající tabulku)
    base  = regexprep(resp,'_SUL_LOG$','');
    sulVar = [base '_SUL'];
    rawVar = base;

    % 1) Musíme mít raw intenzitu (BASE); bez ní RAW prezentace nejde
    if ~ismember(rawVar, T.Properties.VariableNames)
        f = NaN; return;
    end

    % 2) Když chybí _SUL, zkus ho vytvořit
    if ~ismember(sulVar, T.Properties.VariableNames)
        try
            % ensureSUL_LOG přidá BASE_SUL a BASE_SUL_LOG, je NaN-tolerantní
            [T, ~] = azvpet.features.ensureSUL_LOG(T, base);
        catch
            % necháme f = NaN níže
        end
    end

    % 3) Pokud už _SUL existuje, spočítej robustní poměr
    if ismember(sulVar, T.Properties.VariableNames)
        num = double(T.(rawVar));
        den = double(T.(sulVar));
        msk = isfinite(num) & isfinite(den) & den > 0;
        if nnz(msk) >= 10
            r = num(msk) ./ den(msk);
            f = median(r, 'omitnan');
            if ~isfinite(f) || f<=0, f = NaN; end
        else
            f = NaN;
        end
    else
        f = NaN;
    end
end
function map = raw_from_sul_mapping(Tin, resp)
% Vrátí mapování RAW ≈ a + b * SUL(orig) pro daný resp (..._SUL_LOG).
% Pokud chybí *_SUL, pokusí se ho vytvořit přes ensureSUL_LOG.
% Výstup: struct('a',a,'b',b,'ok',logical)

    T = Tin;  map = struct('a',NaN,'b',NaN,'ok',false);

    base   = regexprep(resp,'_SUL_LOG$','');
    sulVar = [base '_SUL'];
    rawVar = base;

    if ~ismember(rawVar, T.Properties.VariableNames)
        return; % RAW není k dispozici -> nelze mapovat
    end
    if ~ismember(sulVar, T.Properties.VariableNames)
        try
            [T, ~] = azvpet.features.ensureSUL_LOG(T, base); % dopočti SUL
        catch
            return;
        end
    end

    % SUL(orig) = exp(SUL_LOG) * smear  — smear si vezmeme z 'cal' při volání
    % Tady stačí „pro data“ použít log→exp bez smear (faktor se chytne do a,b)
    if ismember([sulVar '_LOG'], T.Properties.VariableNames)
        sul_orig = exp(double(T.([sulVar '_LOG']))); % smear se absorbne do a,b
    else
        sul_orig = double(T.(sulVar));
    end

    raw = double(T.(rawVar));
    msk = isfinite(raw) & isfinite(sul_orig);
    if nnz(msk) < 20, return; end

    try
        lm = fitlm(sul_orig(msk), raw(msk));    % raw ~ a + b*sul
        a = lm.Coefficients.Estimate(1);
        b = lm.Coefficients.Estimate(2);
        if isfinite(a) && isfinite(b)
            map.a = a; map.b = b; map.ok = true;
        end
    catch
        % fallback: robustfit (Statistics TB)
        try
            b_ab = robustfit(sul_orig(msk), raw(msk));
            a = b_ab(1); b = b_ab(2);
            if isfinite(a) && isfinite(b)
                map.a = a; map.b = b; map.ok = true;
            end
        catch
        end
    end
end
