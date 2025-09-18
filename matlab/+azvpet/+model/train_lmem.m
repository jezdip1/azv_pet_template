function [M, info] = train_lmem(T, formula, opts)
  responses = string(opts.response_list(:));
  M = containers.Map('KeyType','char','ValueType','any');
  info = struct('responses', responses, 'aic', [], 'bic', [], 'R2m', [], 'R2c', []);

  vnames = string(T.Properties.VariableNames);      % << jednotný typ
  for i = 1:numel(responses)
   resp = string(responses(i));                 % původní (hezky čitelný) název
    col  = azvpet.util.resolve_resp_column(resp, T, opts);   % skutečný sloupec v Tc
    
    vnames = string(T.Properties.VariableNames);
    if ~ismember(col, vnames)
        warning('train_lmem: skip %s -> response column missing in T.', char(resp));
        continue
    end
    
    % Fit vždy na 'col' (to je název v Tc); ukládej pod původním klíčem 'resp'
    frm = sprintf('%s ~ %s', char(col), char(formula));
    try
        mdl = fitlme(T, frm, 'FitMethod','REML');
    catch ME
        warning('train_lmem: skip %s -> %s', char(resp), ME.message);
        continue
    end
    M(char(resp)) = mdl;
    
    [R2m, R2c] = azvpet.util.stats_utils.r2_lmm(mdl);
    info.aic(end+1,1) = mdl.ModelCriterion.AIC;
    info.bic(end+1,1) = mdl.ModelCriterion.BIC;
    info.R2m(end+1,1) = R2m;
    info.R2c(end+1,1) = R2c;

  end
end
