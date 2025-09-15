function report_model(T, M, info, cv, cal, paths)
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
                xlabel('Estimate (link)'); title([resp ' — Fixed effects (95% CI)']); grid on
                saveas(f, fullfile(out,[resp '_fixed_forest.png'])); close(f);
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
                ylabel(strrep(regexprep(resp,'_LOG$',''), '_','-'));
                xlabel('Age [y]'); grid on
                if hasSex, legend(["PI","CI","Mean"] + " | Sex=" + sexCats, 'Location','best'); end
                title([resp ' — Age effect (orig scale, CI & PI)']);
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
                xlabel('Age [y]'); ylabel('z = (Obs - Pred) / (c * sd_{pred})');
                title([resp ' — Standardized residuals']); grid on
                saveas(f, fullfile(out,[resp '_zscore_vs_age.png'])); close(f);
            end
        catch ME
            warning('Z-score plot (%s): %s', resp, ME.message);
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
