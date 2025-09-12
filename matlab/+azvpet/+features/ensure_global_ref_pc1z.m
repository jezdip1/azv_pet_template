function T = ensure_global_ref_pc1z(T)
% Compute PC1 across all Mean_*_SUL_LOG columns; store z-score as GlobalRef_PC1_z
vn = string(T.Properties.VariableNames);
cols = vn(~cellfun(@isempty, regexp(vn,'^Mean_.*_SUL_LOG$','once')));
if isempty(cols); return; end
X = T(:,cols);
% numeric matrix, omit rows with any NaN for PCA fit (simple baseline)
M = X{:,:};
rm = any(~isfinite(M),2);
[coeff, score, ~, ~, ~, mu] = pca(M(~rm,:), 'Centered',true);
pc1 = nan(height(T),1);
pc1(~rm) = score(:,1);
% z-score of PC1
z = (pc1 - mean(pc1(~rm),'omitnan')) ./ std(pc1(~rm),'omitnan');
T.GlobalRef_PC1_z = z;
end

