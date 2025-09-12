function [M, info] = train_lmem(T, formula, opts)
  responses = opts.response_list;
  M = containers.Map('KeyType','char','ValueType','any');
  info = struct('responses', responses, 'aic', [], 'bic', [], 'R2m', [], 'R2c', []);
  for r = responses'
    resp = char(r);
    frm = sprintf('%s ~ %s', resp, formula);
    mdl = fitlme(T, frm, 'FitMethod','REML');
    M(resp) = mdl;
    aic = mdl.ModelCriterion.AIC; bic = mdl.ModelCriterion.BIC;
    [R2m, R2c] = azvpet.util.stats_utils.r2_lmm(mdl);
    info.aic(end+1,1) = aic; info.bic(end+1,1) = bic; info.R2m(end+1,1)=R2m; info.R2c(end+1,1)=R2c;
  end
end
