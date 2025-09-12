function C = load_model_config(path)
  txt = fileread(path);
  C = jsondecode(txt);
end