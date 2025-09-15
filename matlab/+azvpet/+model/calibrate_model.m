function cal = calibrate_model(T, M, cv, nameMap)
% Kalibrace pro všechny modely v M.
% - T: tabulka s validními názvy (po ensure_valid_varnames)
% - M: containers.Map (key = název DV, value = LinearMixedModel)
% - cv: struct s poli pojmenovanými stejně jako klíče M (výstup z loo_cv/loo_cv_par)
% - nameMap (volitelně): table s 'orig','clean' z ensure_valid_varnames
%
% Výstup:
%   cal.(respClean) = struct s poli: smear, sigma_std, sigma_rob,
%                     alpha, alphaCI, beta, betaCI, R2_link, n_train, n_loo_pairs
%   Pokud je nameMap, přidá i cal.(respClean).orig s původním jménem.

    % 1) klíče z M = pravda o tom, co se opravdu fitnulo
    k = keys(M);
    cal = struct();

    for i = 1:numel(k)
        resp = k{i};                  % čisté jméno (validní identifikátor)
        mdl  = M(resp);

        % ==== in-sample residua na link-škále (population) ====
        e = residuals(mdl, 'Conditional', false);
        e = double(e(:));
        s_rob = 1.4826 * mad(e, 1);
        s_std = std(e, 'omitnan');

        % Smearing (jen pro *_LOG)
        isLogDV = endsWith(resp, '_SUL_LOG') || endsWith(resp, '_LOG');
        if isLogDV
            smear = mean(exp(e), 'omitnan');
        else
            smear = NaN;
        end

        % ==== LOO link-kalibrace: y_true ~ alpha + beta*y_pred ====
        alpha = NaN; beta = NaN; R2 = NaN;
        alphaCI = [NaN NaN]; betaCI = [NaN NaN]; n_loo = 0;

        if isstruct(cv) && isfield(cv, resp)
            TT = cv.(resp);
        else
            % případně zkus sanitizaci (kolize by neměly být, ale pro jistotu)
            fbackup = matlab.lang.makeValidName(resp,'ReplacementStyle','underscore');
            if isstruct(cv) && isfield(cv, fbackup)
                TT = cv.(fbackup);
            else
                TT = [];
            end
        end

        if ~isempty(TT)
            y_true = double(TT.y_true);
            y_pred = double(TT.y_pred);
            v = isfinite(y_true) & isfinite(y_pred);
            n_loo = nnz(v);
            if n_loo >= 3
                lm = fitlm(y_pred(v), y_true(v));
                ci = coefCI(lm, 0.05);
                alpha   = lm.Coefficients.Estimate(1);
                beta    = lm.Coefficients.Estimate(2);
                alphaCI = ci(1,:);
                betaCI  = ci(2,:);
                R2      = lm.Rsquared.Ordinary;
            end
        end

        S = struct( ...
            'smear',      smear, ...
            'sigma_std',  s_std, ...
            'sigma_rob',  s_rob, ...
            'alpha',      alpha, ...
            'alphaCI',    alphaCI, ...
            'beta',       beta, ...
            'betaCI',     betaCI, ...
            'R2_link',    R2, ...
            'n_train',    height(mdl.Variables), ...
            'n_loo_pairs',n_loo ...
        );

        % volitelně ulož i původní jméno (hezké pro reporty)
        if nargin >= 4 && ~isempty(nameMap) && istable(nameMap) ...
           && all(ismember({'orig','clean'}, nameMap.Properties.VariableNames))
            ix = find(string(nameMap.clean) == string(resp), 1);
            if ~isempty(ix)
                S.orig = char(nameMap.orig(ix));
            end
        end

        cal.(resp) = S;
    end
end
