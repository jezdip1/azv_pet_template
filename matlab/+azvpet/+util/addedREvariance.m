function addVar = addedREvariance(L, TT)
% Sečte var-komponenty pro RE skupiny, které jsou v TT „nové“ vůči tréninku v L.
% Bezpečný fallback: když cokoliv chybí -> vrátí nuly.

    n = height(TT);
    addVar = zeros(n,1);

    try
        % R2020b+: vrátí tabulku s odhady var komponent
        [~, Info] = covarianceParameters(L);
        if ~istable(Info)
            try
                Info = dataset2table(Info,'ReadRowNames',false);
            catch
                addVar = zeros(n,1); return;
            end
        end

        % Hledáme interceptové RE (std) pro grouping
        G  = string(Info.('Group'));
        Ty = string(Info.('Type'));
        Nm = string(Info.('Name1'));
        Est= double(Info.('Estimate'));

        mask = Ty=="std" & Nm=="(Intercept)";
        groups = unique(G(mask));

        Train = L.Variables;                     % tréninková data, která model viděl

        for g = groups.'
            gn = char(g);
            if ~ismember(gn, TT.Properties.VariableNames) || ~ismember(gn, Train.Properties.VariableNames)
                continue;
            end
            % odhad směrodatné odchylky -> variance
            s2g = Est(find(G==g & mask, 1)).^2;

            % nové úrovně groupingu (v TT, které nebyly při tréninku)
            catsTrain = removecats(categorical(Train.(gn)));
            isNew = ~ismember(categorical(TT.(gn)), categories(catsTrain));
            addVar = addVar + s2g .* double(isNew);
        end
    catch
        % Fallback: neznáme komponenty -> vrať nuly (Conditional=false už randomy nulují)
        addVar = zeros(n,1);
    end
end
