function startup
  % Spouštěj z kořene repa: run('matlab/startup.m')
  thisDir = fileparts(mfilename('fullpath'));   % .../AZV-PET/matlab
  addpath(thisDir);
  rehash path; rehash toolboxcache;
  fprintf('[startup] AZV-PET: added %s\n', thisDir);
end
