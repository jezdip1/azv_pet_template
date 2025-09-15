function cal = calibrate_from_cv(M, cv)
% Kalibrace pro všechny responses, které jsou opravdu v cv.
% - M: containers.Map (respName -> LinearMixedModel), může chybět pro některé klíče
% - cv: struct s poli 'ResponseName' (tabulky s y_true, y_pred, resid) a polem 'grouping'
%
% Výstup:
%   cal.(resp) = struct: smear, sigma_std, sigma_rob, alpha, alphaCI, beta, betaCI, R2_link,
%                        n_loo_pairs, n_train (je-li model)

    % 0) vylistuj responses podle 'cv' (ignoruj pomocná pole)
    f = fieldnames(cv);
    f = f(~strcmp(f,'grouping'));  % pouze DV pole

    cal = struct();

    for i = 1:numel(f)
        resp = f{i};
        TT = cv.(resp);

        % -- valid pairs (link-škála) --
        y_true = double(TT.y_true);
        y_pred = double(TT.y_pred);
        v = isfinite(y_true) & isfinite(y_pred);
        n_loo = nnz(v);

        alpha = NaN; beta = NaN; R2 = NaN; alphaCI = [NaN NaN]; betaCI = [NaN NaN];

        if n_loo >= 3
            lm = fitlm(y_pred(v), y_true(v));
            ci = coefCI(lm, 0.05);
            alpha   = lm.Coefficients.Estimate(1);
            beta    = lm.Coefficients.Estimate(2);
            alphaCI = ci(1,:);
            betaCI  = ci(2,:);
            R2      = lm.Rsquared.Ordinary;
        end

        % -- rezidua (link) pro sigma a případně smearing --
        r = double(TT.resid);
        r = r(isfinite(r));
        sigma_std = std(r, 'omitnan');
        sigma_rob = 1.4826 * mad(r, 1);

        % -- smearing: preferuj in-sample residuals z mdl, jinak LOO residuals --
        smear = NaN;
        n_train = NaN;
        if isKey(M, resp)
            mdl = M(resp);
            e = residuals(mdl, 'Conditional', false);
            e = double(e(:));
            if endsWith(resp,'_LOG') || endsWith(resp,'_SUL_LOG')
                smear = mean(exp(e), 'omitnan');
            end
            n_train = height(mdl.Variables);
        else
            % fallback: LOO rezidua (link) → smearing je definovatelný stejně
            if endsWith(resp,'_LOG') || endsWith(resp,'_SUL_LOG')
                smear = mean(exp(r), 'omitnan');
            end
        end

        cal.(resp) = struct( ...
            'smear',      smear, ...
            'sigma_std',  sigma_std, ...
            'sigma_rob',  sigma_rob, ...
            'alpha',      alpha, ...
            'alphaCI',    alphaCI, ...
            'beta',       beta, ...
            'betaCI',     betaCI, ...
            'R2_link',    R2, ...
            'n_train',    n_train, ...
            'n_loo_pairs',n_loo ...
        );
    end
end
