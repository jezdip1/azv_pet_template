function tests = test_feature_engineering; tests = functiontests(localfunctions); end
function test_minimal(~)
  T = table(string({'p1';'p2'}), [1;2], 'VariableNames',{'PatientID','age'});
  cfg = struct('ensure', struct('SUL_LOG',false,'asymmetry',false,'global_ref_pc1z',false,'metadata_clean',true));
  T2 = azvpet.features.ensure_features(T, cfg);
  assert(height(T2)==2);
end
