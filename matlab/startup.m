function startup
  % --> Tohle je přesná cesta na složku 'matlab' (kde leží tento soubor)
  thisDir = fileparts(mfilename('fullpath'));   % .../AZV-PET/matlab

  % Přidej JEN tuto složku (parent balíčku), ne package samotný
  addpath(thisDir);

  % Ověření, že existuje package a klíčový soubor
  pkgDir = fullfile(thisDir, '+azvpet');
  assert(exist(pkgDir, 'dir')==7, 'Nenalezena složka balíčku: %s', pkgDir);

  mdlFile = fullfile(pkgDir, 'model', 'train_lmem.m');
  assert(exist(mdlFile, 'file')==2, 'Chybí soubor: %s', mdlFile);

  rehash path; rehash toolboxcache;
  fprintf('[startup] AZV-PET: added %s\n', thisDir);
end
