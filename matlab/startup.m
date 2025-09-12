function startup
  root = fileparts(fileparts(mfilename('fullpath')));
  addpath(fullfile(root,'matlab'));
  addpath(genpath(fullfile(root,'matlab','+azvpet')));
  fprintf('[startup] AZV-PET paths added.\n');
end
