% function [T, info] = ensure_GlobalRefPC1_z(T, params_dir)
% % Načte projektor z params_dir/globalref_pca.mat a spočte GlobalRefPC1_z.
% % Pokud něco chybí, vrátí NaN.
% 
%     if nargin<2 || isempty(params_dir), params_dir = './models/_globals'; end
%     f = fullfile(params_dir, 'globalref_pca.mat');
%     info = struct('used',false,'file',f,'n_ok',0);
%     if ~isfile(f)
%         warning('ensure_GlobalRefPC1_z: nenalezen %s -> GlobalRefPC1_z = NaN', f);
%         if ~ismember('GlobalRefPC1_z', string(T.Properties.VariableNames))
%             T.GlobalRefPC1_z = nan(height(T),1);
%         end
%         return;
%     end
% 
%     S = load(f);  % obsahuje: refCols, mu, sd, coef1, pc1_mu, pc1_sd, result_var
%     refCols = string(S.refCols);
%     miss = refCols(~ismember(refCols, string(T.Properties.VariableNames)));
%     if ~isempty(miss)
%         warning('ensure_GlobalRefPC1_z: chybí sloupce: %s', strjoin(cellstr(miss),', '));
%         % vyrob prázdný sloupec, ať predikce nespadne – ale budou NaN
%         if ~ismember('GlobalRefPC1_z', string(T.Properties.VariableNames))
%             T.GlobalRefPC1_z = nan(height(T),1);
%         end
%         return;
%     end
% 
%     X = T{:, refCols};
%     ok = all(isfinite(X),2);
%     Z = nan(height(T),1);
%     if any(ok)
%         Xz = (X(ok,:) - S.mu) ./ S.sd;
%         pc1 = Xz * S.coef1;
%         Z(ok) = (pc1 - S.pc1_mu) ./ S.pc1_sd;
%         info.n_ok = nnz(ok);
%     end
% 
%     varname = 'GlobalRefPC1_z';   % sjednoťme název (bez podtržítka)
%     T.(varname) = Z;
%     info.used = true;
% end
function [T, info] = ensure_GlobalRefPC1_z(T, params_dir, mode, refBases)
% ensure_GlobalRefPC1_z
% - FIT  : ze zadaných referenčních *_SUL_LOG sloupců spočte PCA(PC1),
%          uloží projektor do params_dir/globalref_pca.mat a vrátí GlobalRefPC1_z.
% - PRED : načte uložený projektor a dopočítá GlobalRefPC1_z.
%
% Volání:
%   [T,info] = ensure_GlobalRefPC1_z(T, './models/_globals','fit',  string(cfg.global_ref.refs));
%   [T,info] = ensure_GlobalRefPC1_z(T, './models/_globals','predict');
%
% Pozn:
% - Výstupní sloupec se jmenuje vždy 'GlobalRefPC1_z' (sjednoceno).
% - Pokud něco chybí (columns/param soubor), funkce nepadá → naplní NaN a varuje.

    if nargin < 2 || isempty(params_dir), params_dir = './models/_globals'; end
    if nargin < 3 || isempty(mode),       mode = 'predict';                end
    if nargin < 4,                         refBases = strings(0,1);         end

    info = struct('mode',mode,'file','', 'refCols',string.empty(1,0), ...
                  'n_ok',0,'saved',false,'used',false);

    if ~isfolder(params_dir), mkdir(params_dir); end
    f = fullfile(params_dir, 'globalref_pca.mat');
    info.file = f;

    % sjednoť název výstupní proměnné
    outVar = 'GlobalRefPC1_z';
    if ~ismember(outVar, string(T.Properties.VariableNames))
        T.(outVar) = nan(height(T),1); % předvyplň
    end

    switch lower(mode)
        case 'fit'
            % ---- 1) vyber referenční *_SUL_LOG sloupce ----
            vn = string(T.Properties.VariableNames);

            % když nejsou refBases, pokusíme se rozumně odhadnout
            if isempty(refBases)
                % preferuj Mean_*_SUL_LOG; fallback i na Median_*_SUL_LOG
                colsMean   = vn(~cellfun(@isempty, regexp(vn,'^Mean_.*_SUL_LOG$','once')));
                colsMedian = vn(~cellfun(@isempty, regexp(vn,'^Median_.*_SUL_LOG$','once')));
                refCols = colsMean;
                if numel(refCols) < 2, refCols = colsMedian; end
            else
                refCols = refBases(:) + "_SUL_LOG";
                refCols = refCols(ismember(refCols, vn));
            end

            if numel(refCols) < 2
                warning('ensure_GlobalRefPC1_z[FIT]: nemám ≥2 referenční *_SUL_LOG sloupce → ponechávám NaN.');
                return;
            end
            info.refCols = refCols;

            % ---- 2) příprava matice a standardizace ----
            X = T{:, refCols};
            ok = all(isfinite(X),2);
            if nnz(ok) < 10
                warning('ensure_GlobalRefPC1_z[FIT]: málo platných řádků pro PCA (n=%d).', nnz(ok));
                return;
            end

            mu = mean(X(ok,:), 1, 'omitnan');
            sd = std (X(ok,:), 0, 'omitnan');
            % ochrana proti nulové varianci
            sd(~isfinite(sd) | sd==0) = 1;

            Zok = (X(ok,:) - mu) ./ sd;

            % ---- 3) PCA na standardizovaných datech (už necentruj) ----
            [coeff, score, ~, ~, ~] = pca(Zok, 'Centered', false);

            pc1      = score(:,1);
            coef1    = coeff(:,1);
            pc1_mu   = mean(pc1, 'omitnan');
            pc1_sd   = std (pc1, 0, 'omitnan'); if ~isfinite(pc1_sd) || pc1_sd==0, pc1_sd = 1; end

            % ---- 4) ulož parametry projektoru ----
            result_var = outVar;
            save(f, 'refCols','mu','sd','coef1','pc1_mu','pc1_sd','result_var');
            info.saved = true;

            % ---- 5) vypočti GlobalRefPC1_z i pro tréninková data ----
            Zall = nan(height(T),1);
            if any(ok)
                Zstd = (X(ok,:) - mu) ./ sd;
                pc1_all = Zstd * coef1;
                Zall(ok) = (pc1_all - pc1_mu) ./ pc1_sd;
                info.n_ok = nnz(ok);
            end
            T.(outVar) = Zall;
            info.used = true;

        case 'predict'
            % ---- predikční větev: načti projektor a aplikuj ----
            if ~isfile(f)
                warning('ensure_GlobalRefPC1_z[PRED]: nenalezen %s → %s = NaN.', f, outVar);
                return;
            end
            S = load(f);  % musí obsahovat: refCols, mu, sd, coef1, pc1_mu, pc1_sd, result_var
            req = {'refCols','mu','sd','coef1','pc1_mu','pc1_sd'};
            if any(~isfield(S, req))
                warning('ensure_GlobalRefPC1_z[PRED]: soubor %s neobsahuje potřebná pole → %s = NaN.', f, outVar);
                return;
            end

            refCols = string(S.refCols);
            miss = refCols(~ismember(refCols, string(T.Properties.VariableNames)));
            if ~isempty(miss)
                warning('ensure_GlobalRefPC1_z[PRED]: chybí sloupce: %s → %s = NaN.', strjoin(cellstr(miss),', '), outVar);
                return;
            end

            X  = T{:, refCols};
            ok = all(isfinite(X),2);
            Z  = nan(height(T),1);
            if any(ok)
                Xz = (X(ok,:) - S.mu) ./ S.sd;
                pc1 = Xz * S.coef1;
                Z(ok) = (pc1 - S.pc1_mu) ./ S.pc1_sd;
                info.n_ok = nnz(ok);
            end

            T.(outVar) = Z;
            info.refCols = refCols;
            info.used    = true;

        otherwise
            error('ensure_GlobalRefPC1_z: neznámý mode "%s" (použij "fit" nebo "predict").', mode);
    end
end
