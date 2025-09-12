function save_artifacts(M, info, cal, paths, featCfg, mdlCfg)
  if ~exist(paths.artifacts_dir,'dir'), mkdir(paths.artifacts_dir); end
  save(fullfile(paths.artifacts_dir,'lmem_model.mat'), 'M','info','-v7.3');
  save(fullfile(paths.artifacts_dir,'feature_map.mat'), 'featCfg','-v7.3');
  fid = fopen(fullfile(paths.artifacts_dir,'calibration.json'),'w'); fwrite(fid, jsonencode(cal,'PrettyPrint',true)); fclose(fid);
  fid = fopen(fullfile(paths.artifacts_dir,'model_config.json'),'w'); fwrite(fid, jsonencode(mdlCfg,'PrettyPrint',true)); fclose(fid);
end
