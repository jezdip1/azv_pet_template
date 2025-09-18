function [T, info] = ensure_GlobalRefPC1_z(T, params_dir)
% Načte projektor z params_dir/globalref_pca.mat a spočte GlobalRefPC1_z.
% Pokud něco chybí, vrátí NaN.

    if nargin<2 || isempty(params_dir), params_dir = './models/_globals'; end
    f = fullfile(params_dir, 'globalref_pca.mat');
    info = struct('used',false,'file',f,'n_ok',0);
    if ~isfile(f)
        warning('ensure_GlobalRefPC1_z: nenalezen %s -> GlobalRefPC1_z = NaN', f);
        if ~ismember('GlobalRefPC1_z', string(T.Properties.VariableNames))
            T.GlobalRefPC1_z = nan(height(T),1);
        end
        return;
    end

    S = load(f);  % obsahuje: refCols, mu, sd, coef1, pc1_mu, pc1_sd, result_var
    refCols = string(S.refCols);
    miss = refCols(~ismember(refCols, string(T.Properties.VariableNames)));
    if ~isempty(miss)
        warning('ensure_GlobalRefPC1_z: chybí sloupce: %s', strjoin(cellstr(miss),', '));
        % vyrob prázdný sloupec, ať predikce nespadne – ale budou NaN
        if ~ismember('GlobalRefPC1_z', string(T.Properties.VariableNames))
            T.GlobalRefPC1_z = nan(height(T),1);
        end
        return;
    end

    X = T{:, refCols};
    ok = all(isfinite(X),2);
    Z = nan(height(T),1);
    if any(ok)
        Xz = (X(ok,:) - S.mu) ./ S.sd;
        pc1 = Xz * S.coef1;
        Z(ok) = (pc1 - S.pc1_mu) ./ S.pc1_sd;
        info.n_ok = nnz(ok);
    end

    varname = 'GlobalRefPC1_z';   % sjednoťme název (bez podtržítka)
    T.(varname) = Z;
    info.used = true;
end
