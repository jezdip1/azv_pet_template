function [R2m, R2c] = r2_lmm(mdl)
% Nakagawa & Schielzeth R^2 for LinearMixedModel (link scale).
% - R2m: var(FE) / (var(FE)+var(RE)+var(E))
% - R2c: (var(FE)+var(RE)) / (var(FE)+var(RE)+var(E))

    % Fixed-effects fitted values (bez náhodných efektů)
    try
        muF = predict(mdl, mdl.Variables, 'Conditional', false);
    catch
        % starší verze: fallback – aproximace
        muF = fitted(mdl) - random(mdl); %#ok<RAND>
    end
    varF = var(muF, 1, 'omitnan');

    % Reziduální variance
    varE = mdl.MSE;

    % Variance náhodných efektů (součet přes všechny RE komponenty)
    varR = 0;
    try
        [Psi, Info] = covarianceParameters(mdl);  % metoda, ne pole
        if istable(Info)
            % preferuj řádky typu 'std' a umocni na druhou
            if any(strcmpi(Info.Properties.VariableNames,'Type')) && any(strcmpi(Info.Properties.VariableNames,'Estimate'))
                isStd = strcmpi(string(Info.Type), 'std');
                est   = double(Info.Estimate(isStd));
                if ~isempty(est)
                    varR = nansum(est.^2);
                end
            end
        else
            % dataset/struct fallback
            try
                T = dataset2table(Info);
            catch
                T = struct2table(Info);
            end
            if ismember('Type', T.Properties.VariableNames) && ismember('Estimate', T.Properties.VariableNames)
                isStd = strcmpi(string(T.Type),'std');
                est   = double(T.Estimate(isStd));
                varR  = nansum(est.^2);
            end
        end

        % Když Info neobsahuje std, dopočti z Psi (součet diagonál)
        if ~(isfinite(varR) && varR>0)
            varR = 0;
            if iscell(Psi)
                for i = 1:numel(Psi)
                    S = Psi{i};
                    if isnumeric(S) && ~isempty(S)
                        varR = varR + nansum(diag(S));
                    end
                end
            elseif isnumeric(Psi) && ~isempty(Psi)
                varR = varR + nansum(diag(Psi));
            end
        end
    catch
        % nic – nech varR=0 (Conditional=false už RE nuloval)
        varR = 0;
    end

    denom = varF + varR + varE;
    if ~isfinite(denom) || denom<=0
        R2m = NaN; R2c = NaN;
    else
        R2m = varF / denom;
        R2c = (varF + varR) / denom;
    end
end
