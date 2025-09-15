function cv = loo_cv_par(T, formulaRHS, opts)
    % --- responses (cellstr) ---
    rs = opts.response_list;
    if ischar(rs), responses = {rs};
    elseif isstring(rs), responses = cellstr(rs(:));
    elseif iscellstr(rs), responses = rs(:);
    else, error('loo_cv_par: response_list must be char/string/cellstr'); end
    nResp = numel(responses);

    % --- grouping ---
    if ~isfield(opts,'grouping'), error('loo_cv_par: missing opts.grouping'); end
    grp = string(opts.grouping);
    if ~ismember(grp, string(T.Properties.VariableNames))
        error('loo_cv_par: grouping var %s not in table', grp);
    end
    gid  = string(T.(grp));
    pids = unique(gid,'stable');
    nG   = numel(pids);

    rhs = char(formulaRHS);

    % --- šablona tabulky a kontejnery ---
    varNames = [grp, "y_true","y_pred","resid","MAE","RMSE"];
    varTypes = {'string','double','double','double','double','double'};
    emptyTab = table('Size',[0 numel(varNames)], ...
                     'VariableTypes',varTypes, 'VariableNames',cellstr(varNames));
    parts = cell(nResp, nG);

    % === klíčová oprava: použij Constant pro T i gid ===
    Tc   = parallel.pool.Constant(T);
    gidC = parallel.pool.Constant(gid);

    % --- PARFOR přes skupiny ---
    parfor gi = 1:nG
        Tloc   = Tc.Value;             % vždy tabulka
        gidloc = gidC.Value;           % vždy string vektor

        mask_test = (gidloc == pids(gi));
        Ttr = Tloc(~mask_test,:);
        Tte = Tloc(mask_test,:);

        localParts = cell(1, nResp);

        for ri = 1:nResp
            resp = responses{ri};
            tab  = emptyTab;

            if ismember(resp, Tloc.Properties.VariableNames)
                frm = sprintf('%s ~ %s', resp, rhs);
                try
                    mdl    = fitlme(Ttr, frm, 'FitMethod','REML');
                    y_pred = predict(mdl, Tte, 'Conditional', false);
                    y_true = double(Tte.(resp));
                    resid  = y_true - y_pred;
                    MAE  = mean(abs(resid), 'omitnan');
                    RMSE = sqrt(mean(resid.^2, 'omitnan'));

                    tab  = table(repmat(pids(gi), numel(y_true),1), y_true, y_pred, resid, ...
                                 repmat(MAE,  numel(y_true),1), repmat(RMSE, numel(y_true),1), ...
                                 'VariableNames', cellstr(varNames));
                catch
                    y_true = double(Tte.(resp));
                    tab  = table(repmat(pids(gi), numel(y_true),1), y_true, ...
                                 nan(size(y_true)), nan(size(y_true)), ...
                                 nan(size(y_true)), nan(size(y_true)), ...
                                 'VariableNames', cellstr(varNames));
                end
            end
            localParts{ri} = tab;
        end

        % zapsat do sdíleného „parts“
        for ri = 1:nResp
            parts{ri, gi} = localParts{ri};
        end
    end

    % --- složení výstupu ---
    cv = struct();
    for ri = 1:nResp
        resp = responses{ri};
        cv.(resp) = vertcat(parts{ri,:});
    end
    cv.grouping = char(grp);
end
