function [M, info] = train_lmem(T, formula, opts)
  responses = string(opts.response_list(:));
  M = containers.Map('KeyType','char','ValueType','any');
  info = struct('responses', responses, 'aic', [], 'bic', [], 'R2m', [], 'R2c', []);
  for i = 1:numel(responses)
    resp = responses(i);
    safe = azvpet.util.safe_resp_name(resp);

    % pokud safe != resp, vytvoř v T duplikát pro fit
    if ~strcmp(safe, resp) && ~ismember(safe, T.Properties.VariableNames) ...
                            &&  ismember(resp, T.Properties.VariableNames)
        T.(safe) = T.(resp);
    end

    % fit na SAFE jméno v LHS:
    frm = sprintf('%s ~ %s', safe, formula);
    try
        mdl = fitlme(T, frm, 'FitMethod','REML');
    catch ME
        warning('train_lmem: skip %s -> %s', char(resp), ME.message);
        continue
    end

    % ulož pod PŮVODNÍM jménem jako klíčem:
    M(char(resp)) = mdl;

    aic = mdl.ModelCriterion.AIC; bic = mdl.ModelCriterion.BIC;
    [R2m, R2c] = azvpet.util.stats_utils.r2_lmm(mdl);
    % [R2m, R2c] = azvpet.util.r2_lmm(mdl);
    info.aic(end+1,1) = aic; info.bic(end+1,1) = bic;
    info.R2m(end+1,1)=R2m;   info.R2c(end+1,1)=R2c;
  end
end
