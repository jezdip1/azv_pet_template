function cv = loo_cv(T, formulaRHS, opts)
    % --- responses -> cellstr (sloupcově) ---
    rs = opts.response_list;
    if ischar(rs), responses = {rs};
    elseif isstring(rs), responses = cellstr(rs(:));
    elseif iscellstr(rs), responses = rs(:);
    else, error('loo_cv: response_list must be char/string/cellstr');
    end
    nResp = numel(responses);

    % --- grouping ---
    if ~isfield(opts,'grouping'), error('loo_cv: missing opts.grouping'); end
    grp = string(opts.grouping);
    if ~ismember(grp, string(T.Properties.VariableNames))
        error('loo_cv: grouping var %s not in table', grp);
    end
    gid  = string(T.(grp));
    pids = unique(gid,'stable');

    % --- připrav mapu jmen: orig -> valid field ---
    origNames = string(responses(:));
    cand = arrayfun(@(s) matlab.lang.makeValidName(char(s), 'ReplacementStyle','underscore'), origNames, 'uni', 0);
    % ošetři kolize (kdyby dvě různá orig dala stejný field)
    validFields = matlab.lang.makeUniqueStrings(string(cand), {}, namelengthmax);
    % init structu a mapy
    cv = struct();
    varNames = [grp, "y_true","y_pred","resid","MAE","RMSE"];
    varTypes = [{'string'},{'double'},{'double'},{'double'},{'double'},{'double'}];
    for i = 1:nResp
        f = char(validFields(i));
        cv.(f) = table('Size',[0 numel(varNames)], ...
                       'VariableTypes', varTypes, ...
                       'VariableNames', cellstr(varNames));
    end
    % ulož mapu do structu (pro snadný přístup)
    cv.name_map = table(origNames, validFields, 'VariableNames', {'orig','field'});
    cv.grouping = string(grp);

    % --- LOO ---
    rhs = char(formulaRHS);
    fprintf('\nLOO:    ');
    for gi = 1:numel(pids)
        fprintf('\b\b\b\b%3.0f%%', gi/numel(pids)*100);
        mask_test = (gid == pids(gi));
        Ttr = T(~mask_test,:); 
        Tte = T(mask_test,:);

        for i = 1:nResp
            resp = char(origNames(i));
            f    = char(validFields(i));
            if ~ismember(resp, T.Properties.VariableNames), continue; end

            frm = sprintf('%s ~ %s', resp, rhs);
            try
                mdl    = fitlme(Ttr, frm, 'FitMethod','REML');
                y_pred = predict(mdl, Tte, 'Conditional', false);  % population-level
                y_true = double(Tte.(resp));
                resid  = y_true - y_pred;
                MAE  = mean(abs(resid), 'omitnan');
                RMSE = sqrt(mean(resid.^2, 'omitnan'));

                row = table(repmat(pids(gi), numel(y_true),1), y_true, y_pred, resid, ...
                            repmat(MAE,  numel(y_true),1), repmat(RMSE, numel(y_true),1), ...
                            'VariableNames', cellstr(varNames));
                cv.(f) = [cv.(f); row]; %#ok<AGROW>

            catch ME
                % fail-safe: vyplň NaN, ať nepřijdeš o fold
                y_true = double(Tte.(resp));
                nanrow = table(repmat(pids(gi), numel(y_true),1), y_true, ...
                               nan(size(y_true)), nan(size(y_true)), ...
                               nan(size(y_true)),  nan(size(y_true)), ...
                               'VariableNames', cellstr(varNames));
                cv.(f) = [cv.(f); nanrow]; %#ok<AGROW>
                warning('LOO %s: %s', resp, ME.message);
            end
        end
    end
end
