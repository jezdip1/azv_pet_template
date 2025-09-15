function cal = calibrate_from_cv_full(M, cv)
% Pro každý response v 'cv' spočti:
%   - link-kalibraci: alpha, beta, R2, CI
%   - PI kalibrační faktor c (z robustních standardizovaných residuí)
%   - smearing (z in-sample residuí modelu, fallback z LOO residuí)
%   - lokální SD(Age) na link-škále (funkce přes binning+interp)
%
% Výstup: cal.(resp).{alpha,beta,R2,alphaCI,betaCI,c,smear,localSD_fun,localSD_support, n_loo_pairs, n_train}

    f = fieldnames(cv);
    f = f(~strcmp(f,'grouping'));
    cal = struct();

    for i = 1:numel(f)
        resp = f{i};
        TT = cv.(resp);

        y  = double(TT.y_true);
        p  = double(TT.y_pred);
        sd = double(TT.sd_pred_link);

        v = isfinite(y) & isfinite(p);
        n_loo = nnz(v);

        % --- alpha, beta, R2 (link) ---
        alpha=NaN; beta=NaN; R2=NaN; alphaCI=[NaN NaN]; betaCI=[NaN NaN];
        if n_loo>=3
            lm = fitlm(p(v), y(v));
            ci = coefCI(lm, 0.05);
            alpha   = lm.Coefficients.Estimate(1);
            beta    = lm.Coefficients.Estimate(2);
            alphaCI = ci(1,:);
            betaCI  = ci(2,:);
            R2      = lm.Rsquared.Ordinary;
        end

        % --- PI kalibrační faktor c (link) ---
        okSD = v & isfinite(sd) & sd>0;
        if nnz(okSD)>=10
            rstd = (y(okSD) - p(okSD)) ./ sd(okSD);
            c = max(quantile(abs(rstd),0.975)/1.96, 1);
        else
            c = 1;
        end

        % --- smearing (LOG) ---
        smear = NaN; n_train = NaN;
        if isKey(M, resp)
            L = M(resp);
            e = residuals(L,'Conditional',false); e=double(e(:));
            if endsWith(resp,'_LOG') || endsWith(resp,'_SUL_LOG')
                smear = mean(exp(e),'omitnan');
            end
            n_train = height(L.Variables);
        else
            r = double(TT.resid); r=r(isfinite(r));
            if endsWith(resp,'_LOG') || endsWith(resp,'_SUL_LOG')
                smear = mean(exp(r),'omitnan');
            end
        end

        % --- lokální SD(Age) na link-škále z LOO residuí ---
        Age = double(TT.Age);
        okA = isfinite(Age) & isfinite(y) & isfinite(p);
        resid_link = y(okA) - p(okA);
        age_ok     = Age(okA);
        localSD_fun = @(a) localSD_eval(a, age_ok, resid_link);   % viz helper níž
        localSD_support = [min(age_ok), max(age_ok)];

        cal.(resp) = struct( ...
            'alpha',alpha,'alphaCI',alphaCI, ...
            'beta',beta,'betaCI',betaCI,'R2_link',R2, ...
            'c',c,'smear',smear, ...
            'localSD_fun',localSD_fun, ...
            'localSD_support',localSD_support, ...
            'n_loo_pairs',n_loo,'n_train',n_train ...
        );
    end
end

function sdFun = localSD_eval(aQuery, ageVec, residLink)
    % stejný algoritmus jako dřív: binning podle kvantilů + pchip
    msk = isfinite(ageVec) & isfinite(residLink);
    x = ageVec(msk);
    y2= (residLink(msk)).^2;

    if numel(x) < 20 || numel(unique(x)) < 5
        sdFun = sqrt(max(1e-9, mean(y2,'omitnan'))) * ones(size(aQuery));
        return;
    end

    numBins = max(5, min(15, floor(numel(x)/30)));
    edges = quantile(x, linspace(0,1,numBins+1));
    edges = unique(edges);
    if numel(edges)<4
        sdFun = sqrt(max(1e-9, mean(y2,'omitnan'))) * ones(size(aQuery));
        return;
    end
    binIdx = discretize(x, edges);
    nb = max(binIdx);
    ybin = accumarray(binIdx(~isnan(binIdx)), y2(~isnan(binIdx)), [nb 1], @mean, NaN);
    ctrs = (edges(1:end-1)+edges(2:end))/2;

    v = isfinite(ybin) & isfinite(ctrs(:));
    if nnz(v)<3
        sdFun = sqrt(max(1e-9, mean(y2,'omitnan'))) * ones(size(aQuery));
        return;
    end
    ybin_s = movmean(ybin(v),3,'Endpoints','shrink');
    s2q = interp1(ctrs(v), ybin_s, aQuery, 'pchip', 'extrap');
    sdFun = sqrt(max(1e-9, s2q));
end
