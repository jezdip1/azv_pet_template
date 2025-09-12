% End-to-end: load → features → LMEM → LOO → calibration → report → artifacts
azvpet.util.check_requirements();
paths = jsondecode(fileread('config/paths.local.json'));
featCfg = jsondecode(fileread('config/features_config.json'));
mdlCfg  = jsondecode(fileread('config/model_config.json'));
T0 = azvpet.io.load_table(paths.raw_table);
T1 = azvpet.features.ensure_features(T0, featCfg);
[form, opts] = azvpet.model.define_formula(mdlCfg, T1);
[M, info]    = azvpet.model.train_lmem(T1, form, opts);
cv = azvpet.model.loo_cv(T1, form, opts);
cal = azvpet.model.calibrate_model(T1, M, info, cv);
azvpet.model.report_model(T1, M, info, cv, cal, paths);
azvpet.io.save_artifacts(M, info, cal, paths, featCfg, mdlCfg);
