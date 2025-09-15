function cv = loo_cv_par(T, formulaRHS, opts)
% Parallel LOO po skupině (opts.grouping) pro víc responses.
% Vrací struct s poli = response; tabulka s: Group, Age, y_true, y_pred,
% CI_lo_link_mean, CI_hi_link_mean, sd_pred_link, resid.

    % --- responses ---
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

    % --- šablona výstupní tabulky ---
    varNames = [grp,"Age","y_true","y_pred","CI_lo_link_mean","CI_hi_link_mean","sd_pred_link","resid"];
    varTypes = {'string','double','double','double','double','double','double','double'};
    emptyTab = table('Size',[0 numel(varNames)], 'VariableTypes',varTypes, 'VariableNames',cellstr(varNames));

    parts = cell(nResp, nG);

    % --- Constant broadcasty ---
    Tc   = parallel.pool.Constant(T);
    gidC = parallel.pool.Constant(gid);

    % předpočítej info o kategoriích (levels sjednotíme v tr/test)
    vnames = string(T.Properties.VariableNames);
    isCat  = false(size(vnames));
    cats   = cell(size(vnames));
    for j = 1:numel(vnames)
        isCat(j) = iscategorical(T.(vnames(j)));
        if isCat(j), cats{j} = categories(T.(vnames(j))); end
    end
    isCatC = parallel.pool.Constant(isCat);
    catsC  = parallel.pool.Constant(cats);
    vnamesC= parallel.pool.Constant(vnames);

    parfor gi = 1:nG
        Tloc   = Tc.Value;
        gidloc = gidC.Value;
        isCatL = isCatC.Value;
        catsL  = catsC.Value;
        vnamesL= vnamesC.Value;

        mask_test = (gidloc == pids(gi));
        Ttr = Tloc(~mask_test,:);
        Tte = Tloc(mask_test,:);

        % sjednotit kategorie ve všech cat sloupcích
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
                    L = fitlme(Ttr, frm, 'FitMethod','REML');

                    % CI_mean (link): použij predict s CI
                    [y_pred, yCI] = predict(L, Tte, 'Conditional', false, 'Alpha', 0.05);
                    z975 = norminv(0.975);
                    SEm  = (yCI(:,2) - yCI(:,1)) / (2*z975);

                    % predikční SD (link) pro PI: SEm^2 + MSE + added RE variance (nové skupiny)
                    s2res = L.MSE;
                    addVar = addedREvariance(L, Tte);   % helper níž
                    sd_pred_link = sqrt(max(0, SEm.^2 + s2res + addVar));

                    y_true = double(Tte.(resp));
                    resid  = y_true - y_pred;

                    Age = nan(height(Tte),1);
                    if ismember('Age', Tte.Properties.VariableNames)
                        Age = double(Tte.Age);
                    end

                    tab = table(repmat(pids(gi), height(Tte),1), Age, y_true, y_pred, ...
                                yCI(:,1), yCI(:,2), sd_pred_link, resid, ...
                                'VariableNames', cellstr(varNames));
                catch
                    % nech prázdný fold (rychleji se skládá než NaNy)
                    tab = emptyTab;
                end
            end

            localParts{ri} = tab;
        end

        for ri = 1:nResp
            parts{ri, gi} = localParts{ri};
        end
    end

    % slož výstup
    cv = struct();
    for ri = 1:nResp
        resp = responses{ri};
        TT = vertcat(parts{ri,:});
        % odfiltruj řádky bez páru (když byl fold prázdný)
        v = isfinite(TT.y_true) & isfinite(TT.y_pred);
        cv.(resp) = TT(v,:);
    end
    cv.grouping = char(grp);
end

function addVar = addedREvariance(L, TT)
% Přidá var komponentu pro nové úrovně náhodných efektů (intercepty), pokud jsou v TT nové.
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
        % nic – když covarianceParameters selže, necháme 0
    end
end
