function out = predict_new_exam(Tnew, resp, M, cal)
% Jedno vyšetření (řádek) -> predikce a klasifikace.
% Vstup:
%  - Tnew : table se všemi potřebnými kovariátami (valid names!)
%  - resp : název DV (např. 'Median_..._SUL_LOG')
%  - M    : containers.Map s vytrénovanými modely
%  - cal  : struct z calibrate_from_cv_full (pro stejné resp)
%
% Výstup: struct s poli:
%  Pred_link, CI_link (lo,hi), sd_pred_link, PI_link (lo,hi),
%  Pred_orig, CI_orig (lo,hi), PI_orig (lo,hi),
%  z, p_two_sided, is_outlier_95, is_outlier_99

    % if ~isKey(M, resp), error('predict_new_exam: model "%s" not found.', resp); end
    safe = azvpet.util.safe_resp_name(resp);
    if isKey(M, resp)
        L = M(resp);
    elseif isKey(M, safe)
        L = M(safe);
    else
        error('predict_new_exam: model "%s" not found.', resp);
    end
    % L = M(resp);
    C = cal.(resp);

    % population-level (bez RE) – link-škála
    [yhat_link, yCI_link] = predict(L, Tnew, 'Conditional', false, 'Alpha', 0.05);
    z975 = norminv(0.975);
    SEm  = (yCI_link(:,2)-yCI_link(:,1)) / (2*z975);
    s2res = L.MSE;
    addVar = azvpet.util.addedREvariance(L, Tnew);     % stejné jako výše
    sd_pred_link = sqrt(max(0, SEm.^2 + s2res + addVar));

    % kalibrace link-škály (alpha + beta*y)
    yhat_cal_link = C.alpha + C.beta.*yhat_link;

    % PI (link) s kalibračním faktorem c
    PI_lo_link = yhat_cal_link - z975 * (C.c * sd_pred_link);
    PI_hi_link = yhat_cal_link + z975 * (C.c * sd_pred_link);

    % zpět na originální škálu
    isLog = endsWith(resp,'_LOG') || endsWith(resp,'_SUL_LOG');
    if isLog
        smear = C.smear; if ~isfinite(smear) || smear<=0, smear=1; end
        tr  = @(x) exp(x).*smear;
    else
        tr  = @(x) x;
    end
    Pred_orig = tr(yhat_cal_link);
    CI_orig   = tr(yCI_link);         % CI mean (nekalibrujeme)
    PI_orig   = [tr(PI_lo_link), tr(PI_hi_link)];

    % z-score (link) – pokud máme pozorovanou hodnotu (např. pro audit)
    z = NaN; p = NaN; is95=false; is99=false;
    if ismember(resp, Tnew.Properties.VariableNames)
        % yobs_link = double(Tnew.(resp));
        yobs_link = NaN;
        if ismember(resp, Tnew.Properties.VariableNames)
            yobs_link = double(Tnew.(resp));
        elseif ismember(safe, Tnew.Properties.VariableNames)
            yobs_link = double(Tnew.(safe));
        end
        if isfinite(yobs_link) && isfinite(sd_pred_link) && sd_pred_link>0
            z = (yobs_link - yhat_cal_link) ./ (C.c * sd_pred_link);
            p = 2*(1 - normcdf(abs(z)));
            is95 = abs(z) > 1.96;
            is99 = abs(z) > 2.576;
        end
    end

    out = struct( ...
        'Pred_link', yhat_cal_link, ...
        'CI_link', yCI_link, ...
        'sd_pred_link', sd_pred_link, ...
        'PI_link', [PI_lo_link, PI_hi_link], ...
        'Pred_orig', Pred_orig, ...
        'CI_orig', CI_orig, ...
        'PI_orig', PI_orig, ...
        'z', z, 'p_two_sided', p, ...
        'is_outlier_95', is95, 'is_outlier_99', is99 ...
    );
end
