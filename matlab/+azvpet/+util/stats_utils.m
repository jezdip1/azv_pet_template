function [R2m, R2c] = r2_lmm(mdl)
  y = mdl.Variables{:, mdl.ResponseName};
  mu = fitted(mdl);
  varF = var(mu, 1, 'omitnan');
  varR = mdl.covarianceParameters.Estimated;
  if isnumeric(varR)
    varR = sum(varR(:),'omitnan');
  else
    try
      varR = sum(cellfun(@(x) sum(x(:),'omitnan'), varR));
    catch
      varR = NaN;
    end
  end
  varE = mdl.MSE;
  R2m = varF / (varF + varR + varE);
  R2c = (varF + varR) / (varF + varR + varE);
end
