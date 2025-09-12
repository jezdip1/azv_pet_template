function T = ensure_asymmetry(T)
% Create asymmetry indices for Left/Right pairs: (L-R)/(L+R)
vn = string(T.Properties.VariableNames);
lefts = vn(endsWith(vn,"_Left"));
for i = 1:numel(lefts)
    l = lefts(i);
    r = replace(l,"_Left","_Right");
    if ismember(r, vn)
        num = double(T.(l)) - double(T.(r));
        den = double(T.(l)) + double(T.(r));
        ai = num ./ den;
        name = "Asym_" + extractBefore(l,"_Left");
        T.(name) = ai;
    end
end
end

