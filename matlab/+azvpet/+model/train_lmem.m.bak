function [M, info] = train_lmem(T, formula, opts)
% TRAIN_LMEM  Fitne LME pro vícero response proměnných.
%   opts.response_list : char | string | string array | cellstr

    % --- normalizace seznamu responses na cellstr (sloupcový) ---
    rs = opts.response_list;
    if ischar(rs)
        responses = {rs};
    elseif isstring(rs)
        responses = cellstr(rs(:));
    elseif iscellstr(rs)
        responses = rs(:);
    else
        error('train_lmem:BadResponses', ...
              'opts.response_list musí být char/string/cellstr, dostal jsem %s.', class(rs));
    end
    n = numel(responses);

    % --- init výstupů ---
    M = containers.Map('KeyType','char','ValueType','any');
    info = struct();
    info.responses = string(responses(:));
    info.aic = nan(n,1); info.bic = nan(n,1);
    info.R2m = nan(n,1); info.R2c = nan(n,1);
    info.errors = strings(n,1);

    rhs = char(formula);  % fitlme má rád char; rhs už obsahuje i "(1|UNIS)"

    % --- hlavní smyčka přes response proměnné ---
    fprintf('\ntrain LMEM:    ');
    for i = 1:n
        fprintf('\b\b\b\b%3.0f%%', i/n*100);
        resp = responses{i};                 % char název DV
        frm  = sprintf('%s ~ %s', resp, rhs);

        try
            mdl = fitlme(T, frm, 'FitMethod','REML');
            M(resp) = mdl;

            info.aic(i) = mdl.ModelCriterion.AIC;
            info.bic(i) = mdl.ModelCriterion.BIC;

            [R2m, R2c] = r2_lmm(mdl);
            info.R2m(i) = R2m; info.R2c(i) = R2c;

        catch ME
            % neházej celou smyčku; ulož chybu a pokračuj
            info.errors(i) = string(ME.message);
            % volitelně: M(resp) = [];  % mapu lze nechat bez klíče při chybě
        end
    end
end
function [R2m, R2c, parts] = r2_lmm(mdl)
% Nakagawa & Schielzeth R2 pro LMM s random intercepty.
% Vrací i parts.s2_f, parts.s2_u, parts.s2_e pro debug.

    % 1) Řádky použité modelem
    Tuse = mdl.Variables;

    % 2) Variance fixed části (ŷ na link-škále, Conditional=false)
    yhat_fix = predict(mdl, Tuse, 'Conditional', false);
    s2_f = var(double(yhat_fix), 'omitnan');

    % 3) Random variance (součet var komponent pro intercepty)
    s2_u = 0;
    try
        C = covarianceParameters(mdl);   % R2020b+: table
        if istable(C)
            mask = strcmpi(C.Type,'std') & strcmpi(C.Name1,'(Intercept)');
            stds = double(C.Estimate(mask));
            s2_u = nansum(stds.^2);
        else
            % fallback na starší verze
            s2_u = double(mdl.Psi{1}(1));
        end
    catch
        % pokud nic nenajdeme, necháme 0 (žádný RE)
        s2_u = 0;
    end

    % 4) Reziduální variance
    s2_e = double(mdl.MSE);

    % 5) R2
    denom = s2_f + s2_u + s2_e;
    if ~isfinite(denom) || denom<=0
        R2m = NaN; R2c = NaN;
    else
        R2m = s2_f / denom;
        R2c = (s2_f + s2_u) / denom;
    end

    % 6) pro ladění:
    parts = struct('s2_f',s2_f,'s2_u',s2_u,'s2_e',s2_e,'den',denom);
end

