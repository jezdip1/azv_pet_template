function T = ensureSUL_LOG(T, base)
% Deterministic baseline: if column <base> exists and <base>_SUL_LOG missing,
% create it as log(<base>). No side effects otherwise.
% Example: base="Mean_Accumbens_Area_Left" -> creates "Mean_Accumbens_Area_Left_SUL_LOG".

% Allow a string vector or char matrix of names
base = string(base);
vn = string(T.Properties.VariableNames);

for src = base(:)'
    dst = src + "_SUL_LOG";
    if ismember(src, vn) && ~ismember(dst, vn)
        x = T.(src);
        T.(dst) = log(double(x));
        vn(end+1) = dst; % track newly created variable
    end
end
end

