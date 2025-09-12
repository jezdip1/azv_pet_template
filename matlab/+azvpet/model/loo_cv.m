function cv = loo_cv(T, formula, opts)
  grp = string(opts.grouping);
  pid = string(T.(grp));
  pids = unique(pid);
  responses = opts.response_list;
  cv = struct();
  for r = responses'
    key = char(r);
    cv.(key) = table('Size',[0 6], 'VariableTypes',{'string','double','double','double','double','double'}, ...
                     'VariableNames',{'PatientID','y_true','y_pred','resid','MAE','RMSE'});
  end
  for i = 1:numel(pids)
    mask_test = pid==pids(i);
    Ttr = T(~mask_test,:); Tte = T(mask_test,:);
    for r = responses'
      resp = char(r);
      frm = sprintf('%s ~ %s', resp, formula);
      mdl = fitlme(Ttr, frm, 'FitMethod','REML');
      y_true = Tte.(resp);
      y_pred = predict(mdl, Tte);
      resid  = y_true - y_pred;
      MAE = mean(abs(resid),'omitnan'); RMSE = sqrt(mean(resid.^2,'omitnan'));
      cv.(resp) = [cv.(resp); table(pids(i), y_true, y_pred, resid, MAE, RMSE, 'VariableNames',cv.(resp).Properties.VariableNames)]; %#ok<AGROW>
    end
  end
end
