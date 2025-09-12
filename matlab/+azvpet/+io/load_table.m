function T = load_table(path)
  [~,~,ext] = fileparts(path);
  switch lower(ext)
    case '.csv'
      T = readtable(path, 'TextType','string', 'VariableNamingRule','preserve');
    case {'.xlsx','.xls'}
      T = readtable(path, 'TextType','string', 'VariableNamingRule','preserve');
    otherwise
      error('Unsupported table format: %s', ext);
  end
  if any(strcmpi('StudyDate', string(T.Properties.VariableNames)))
    try
      T.StudyDate = datetime(T.StudyDate, 'InputFormat','yyyyMMdd');
    catch
      % nech jak je, zpracujeme později
    end
  end
end
