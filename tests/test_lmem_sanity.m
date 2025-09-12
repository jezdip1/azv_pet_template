function tests = test_lmem_sanity; tests = functiontests(localfunctions); end
function test_fit(~)
  rng default
  N=50; pid = repelem(string(1:10)',5);
  age = randi([20,80],N,1); sex = categorical(randi([0,1],N,1));
  y = 0.02*age + (double(sex)-0.5)*0.1 + randn(N,1)*0.2 + randn(10,1);
  T = table(pid, age, sex, y, 'VariableNames',{'PatientID','age','sex','Mean_Test_SUL_LOG'});
  cfg = struct('by_region',true,'region_pattern','^Mean_.*_SUL_LOG$', 'fixed_effects',{{'age','sex'}}, ...
               'random_effects',struct('grouping','PatientID','terms',{{'1'}}),'categoricals',{{'sex'}});
  [form, opts] = azvpet.model.define_formula(cfg, T);
  [M, info] = azvpet.model.train_lmem(T, form, opts);
  assert(isKey(M,'Mean_Test_SUL_LOG'));
  assert(~isempty(info.R2m));
end
