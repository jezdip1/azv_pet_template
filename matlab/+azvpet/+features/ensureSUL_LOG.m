function T = ensureSUL_LOG(T, base)
  if isstring(base), base = char(base); end
  src = base;
  dst = [base '_SUL_LOG'];
  if ismember(src, T.Properties.VariableNames) && ~ismember(dst, T.Properties.VariableNames)
    x = double(T.(src));
    T.(dst) = log(x);
  end
end


