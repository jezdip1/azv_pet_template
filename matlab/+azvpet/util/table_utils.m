function T = drop_missing(T, required)
  mask = true(height(T),1);
  for v = string(required)
    mask = mask & ~ismissing(T.(v));
  end
  T = T(mask,:);
end
