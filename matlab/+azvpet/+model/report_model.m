function report_model(T, M, info, cv, paths) %#ok<INUSD>
  out = fullfile(paths.reports_dir,'latest'); if ~exist(out,'dir'), mkdir(out); end
  responses = info.responses;
  S = struct();
  for r = responses'
    resp = char(r);
    mdl = M(resp);
    S.(resp) = struct();
    S.(resp).AIC = mdl.ModelCriterion.AIC; S.(resp).BIC = mdl.ModelCriterion.BIC;
    [R2m, R2c] = azvpet.util.stats_utils.r2_lmm(mdl); S.(resp).R2m=R2m; S.(resp).R2c=R2c;
    f = figure('Visible','off'); resid = mdl.residuals('Raw'); histogram(resid); title([resp ' residuals']); saveas(f, fullfile(out,[resp '_resid_hist.png'])); close(f);
    f = figure('Visible','off'); qqplot(resid); title([resp ' QQ-plot']); saveas(f, fullfile(out,[resp '_qq.png'])); close(f);
    cvt = cv.(resp);
    f = figure('Visible','off'); scatter(cvt.y_true, cvt.y_pred, '.'); lsline; xlabel('True'); ylabel('Pred'); title([resp ' LOO']); saveas(f, fullfile(out,[resp '_loo.png'])); close(f);
    S.(resp).LOO = struct('MAE', mean(cvt.MAE), 'RMSE', mean(cvt.RMSE));
  end
  fid = fopen(fullfile(out,'summary.json'),'w'); fwrite(fid, jsonencode(S,'PrettyPrint',true)); fclose(fid);
end
