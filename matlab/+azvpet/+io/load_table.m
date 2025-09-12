function T = load_table(path)
  [~,~,ext] = fileparts(path);
  switch lower(ext)
    case '.csv'
      T = readtable(path, 'TextType','string');
    case {'.xlsx','.xls'}
      T = readtable(path, 'TextType','string');
    otherwise
      error('Unsupported table format: %s', ext);
  end
  if any(strcmpi('StudyDate', T.Properties.VariableNames))
    T.StudyDate = datetime(T.StudyDate, 'InputFormat','yyyyMMdd');
  end
end
