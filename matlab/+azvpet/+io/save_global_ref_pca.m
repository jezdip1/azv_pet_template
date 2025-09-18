function P = save_global_ref_pca(T, refs, outdir)
% Uloží projektor pro GlobalRefPC1_z na základě daných ref sloupců *_SUL_LOG
% P.result_var = 'GlobalRefPC1_z'; P.file uložen na disk

    if ~exist(outdir,'dir'), mkdir(outdir); end
    result_var = 'GlobalRefPC1_z';

    refCols = string(refs(:)) + "_SUL_LOG";
    refCols = refCols(ismember(refCols, string(T.Properties.VariableNames)));
    assert(numel(refCols) >= 2, 'Potřebuji aspoň 2 referenční sloupce.');

    % připrav X (řádky s kompletními hodnotami)
    X = T{:, refCols};
    ok = all(isfinite(X),2);
    Xok = X(ok,:);

    % standardizace po sloupcích (na tréninku!)
    mu = mean(Xok,1,'omitnan');
    sd = std(Xok,0,1,'omitnan'); sd(sd<=0 | ~isfinite(sd)) = 1;
    Xz = (Xok - mu) ./ sd;

    % PCA jen na PC1
    [coef, score] = pca(Xz, 'Centered', false, 'NumComponents', 1);
    pc1 = score(:,1);
    pc1_mu = mean(pc1,'omitnan');
    pc1_sd = std(pc1,0,'omitnan'); if ~isfinite(pc1_sd) || pc1_sd<=0, pc1_sd = 1; end

    P = struct();
    P.refCols      = cellstr(refCols);
    P.mu           = mu;
    P.sd           = sd;
    P.coef1        = coef(:,1);     % sloupcový vektor velikosti numel(refCols)
    P.pc1_mu       = pc1_mu;
    P.pc1_sd       = pc1_sd;
    P.centered     = false;         % už jsme standardizovali sami
    P.result_var   = result_var;

    save(fullfile(outdir,'globalref_pca.mat'), '-struct','P');

    % Volitelně vrať i sloupec do T (pro kontrolu)
    Z = nan(height(T),1);
    Xall = T{:, refCols};
    ok2 = all(isfinite(Xall),2);
    Xz2 = (Xall(ok2,:) - mu) ./ sd;
    pc1_all = Xz2 * P.coef1;
    Z(ok2) = (pc1_all - pc1_mu) ./ pc1_sd;
    T.(result_var) = Z;
end
