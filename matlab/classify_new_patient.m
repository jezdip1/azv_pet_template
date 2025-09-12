% Vstup: CSV s jediným řádkem (nový pacient) nebo tabulka v workspace
paths = jsondecode(fileread('config/paths.local.json'));
featCfg = jsondecode(fileread('config/features_config.json'));
A = load(fullfile(paths.artifacts_dir,'lmem_model.mat')); M=A.M; info=A.info;
cal = jsondecode(fileread(fullfile(paths.artifacts_dir,'calibration.json')));
% Načti nového pacienta
P = azvpet.io.load_table('data/raw/new_patient.csv');
P = azvpet.features.ensure_features(P, featCfg);
out = table; responses = info.responses;
for r = responses'
  resp = char(r); mdl = M(resp);
  y_true = P.(resp);
  y_pred = predict(mdl, P);
  resid = y_true - y_pred;
  s = cal.(resp).sigma_rob; z = resid ./ s; p = 2*normcdf(-abs(z));
  out = [out; table(string(resp), y_true, y_pred, resid, z, p, 'VariableNames',{'Region','True','Pred','Resid','Z','P'})]; %#ok<AGROW>
end
[~, ~, ~, p_fdr] = fdr_bh(out.P);
out.P_fdr = p_fdr;
writetable(out, fullfile(paths.reports_dir,'latest','new_patient_results.csv'));
