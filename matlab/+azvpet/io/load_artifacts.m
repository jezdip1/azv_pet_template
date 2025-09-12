function A = load_artifacts(paths)
  A.M = []; A.info = []; A.cal = [];
  if exist(fullfile(paths.artifacts_dir,'lmem_model.mat'),'file')
    tmp = load(fullfile(paths.artifacts_dir,'lmem_model.mat'));
    A.M = tmp.M; A.info = tmp.info;
  end
  if exist(fullfile(paths.artifacts_dir,'calibration.json'),'file')
    A.cal = jsondecode(fileread(fullfile(paths.artifacts_dir,'calibration.json')));
  end
end
