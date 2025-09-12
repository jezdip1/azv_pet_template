function T = ensureSUL_LOG(T, base)
% Deterministic baseline: if column <base> exists and <base>_SUL_LOG missing,
% create it as log(<base>). No side effects otherwise.
% Example: base="Mean_Accumbens_Area_Left" -> creates "Mean_Accumbens_Area_Left_SUL_LOG".
if isstring(base), base = char(base); end
src = base;
dst = [base '_SUL_LOG'];
vn = string(T.Properties.VariableNames);
if ismember(src, vn) && ~ismember(dst, vn)
    x = T.(src);
    T.(dst) = log(double(x));
end
end

