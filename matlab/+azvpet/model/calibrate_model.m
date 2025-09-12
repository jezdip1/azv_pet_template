function cal = calibrate_model(T, M, info, cv) %#ok<INUSD>
  responses = info.responses;
  cal = struct();
  for r = responses'
    resp = char(r);
    mdl = M(resp);
    y = T.(resp);
    yhat = fitted(mdl);
    resid = y - yhat;
    s_rob = robustScale(resid);
    s_std = std(resid, 'omitnan');
    cal.(resp) = struct('sigma_std', s_std, 'sigma_rob', s_rob, 'n', numel(resid));
  end
end
function s = robustScale(resid)
  s = 1.4826*mad(resid,1);
end
