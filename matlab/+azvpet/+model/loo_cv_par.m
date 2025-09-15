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

       % === před parfor: předpočítej seznam kategoriálních proměnných a jejich levels ===
    vnames = string(T.Properties.VariableNames);
    isCat  = false(size(vnames));
    cats   = cell(size(vnames));
    for j = 1:numel(vnames)
        isCat(j) = iscategorical(T.(vnames(j)));
        if isCat(j)
            cats{j} = categories(T.(vnames(j)));
        end
    end
    isCatC = parallel.pool.Constant(isCat);
    catsC  = parallel.pool.Constant(cats);
    vnamesC= parallel.pool.Constant(vnames);
    
    % --- PARFOR přes skupiny ---
    parfor gi = 1:nG
        Tloc   = Tc.Value;
        gidloc = gidC.Value;
        isCatL = isCatC.Value;
        catsL  = catsC.Value;
        vnamesL= vnamesC.Value;
    
        mask_test = (gidloc == pids(gi));
        Ttr = Tloc(~mask_test,:);
        Tte = Tloc(mask_test,:);
    
        % 2a) Zarovnej kategorie pro všechny kategorické sloupce
        for jj = 1:numel(vnamesL)
            if isCatL(jj)
                nm = vnamesL(jj);
                Ttr.(nm) = categorical(Ttr.(nm), catsL{jj});
                Tte.(nm) = categorical(Tte.(nm), catsL{jj});
            end
        end
    
        localParts = cell(1, nResp);
    
        for ri = 1:nResp
            resp = responses{ri};
            tab  = emptyTab;
    
            if ismember(resp, Tloc.Properties.VariableNames)
                frm = sprintf('%s ~ %s', resp, rhs);
                try
                    mdl = fitlme(Ttr, frm, 'FitMethod','REML');
    
                    % 2b) Seznam prediktorů použité fixní částí
                    preds = string(mdl.Formula.PredictorNames);
    
                    % 2c) maska platných řádků v Tte pro všechny prediktory
                    if isempty(preds)
                        maskValid = true(height(Tte),1);
                    else
                        maskValid = true(height(Tte),1);
                        for pp = 1:numel(preds)
                            nm = preds(pp);
                            if ~ismember(nm, Tte.Properties.VariableNames), continue; end
                            col = Tte.(nm);
                            if isnumeric(col)
                                maskValid = maskValid & isfinite(col);
                            elseif islogical(col)
                                maskValid = maskValid & ~isnan(double(col));
                            elseif iscategorical(col)
                                maskValid = maskValid & ~ismissing(col);
                            elseif isduration(col)
                                maskValid = maskValid & isfinite(seconds(col));
                            elseif isdatetime(col)
                                maskValid = maskValid & ~ismissing(col);
                            else
                                % jiné typy: pokus o ~ismissing
                                try
                                    maskValid = maskValid & ~ismissing(col);
                                catch
                                    % když nevíme, radši vyřadit
                                    maskValid = maskValid & false;
                                end
                            end
                        end
                    end
    
                    if any(maskValid)
                        y_pred = predict(mdl, Tte(maskValid,:), 'Conditional', false);
                        y_true = double(Tte.(resp)(maskValid));
                        resid  = y_true - y_pred;
                        MAE  = mean(abs(resid), 'omitnan');
                        RMSE = sqrt(mean(resid.^2, 'omitnan'));
    
                        tab  = table(repmat(pids(gi), numel(y_true),1), y_true, y_pred, resid, ...
                                     repmat(MAE,  numel(y_true),1), repmat(RMSE, numel(y_true),1), ...
                                     'VariableNames', cellstr(varNames));
                    else
                        % žádný validní řádek v tomhle foldu → prázdná tabulka (ne NaNy)
                        tab = emptyTab;
                    end
    
                catch
                    % fold selhal → prázdná tabulka (ne NaNy)
                    tab = emptyTab;
                end
            end
    
            localParts{ri} = tab;
        end

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
