function tests = test_feature_engineering; tests = functiontests(localfunctions); end
function test_minimal(~)
  T = table(string({'p1';'p2'}), [1;2], [10;20], [12;18], ...
      'VariableNames',{'PatientID','age','Mean_X_Left','Mean_X_Right'});
  cfg = struct('ensure', struct('SUL_LOG',true,'asymmetry',true,'global_ref_pc1z',true,'metadata_clean',true));
  % also include base column to trigger SUL_LOG
  T.Mean_X = [11; 19];
  T2 = azvpet.features.ensure_features(T, cfg);
  assert(all(ismember({'Asym_Mean_X','GlobalRef_PC1_z','Mean_X_SUL_LOG'}, string(T2.Properties.VariableNames))));
end

