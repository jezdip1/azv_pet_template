function [formula, opts] = define_formula(cfg, T)
  if isfield(cfg, 'by_region') && cfg.by_region
    rgx = cfg.region_pattern;
    all = string(T.Properties.VariableNames);
    resp = all(~cellfun(@isempty, regexp(all, rgx)));
    opts.response_list = resp;
  else
    opts.response_list = string(cfg.response);
  end
  if isfield(cfg,'categoricals')
    for c = string(cfg.categoricals)
      if ismember(c, T.Properties.VariableNames)
        T.(c) = categorical(T.(c));
      end
    end
  end
  fx = strjoin(string(cfg.fixed_effects), ' + ');
  if isfield(cfg,'covariate_interactions') && ~isempty(cfg.covariate_interactions)
    fx = fx + " + " + strjoin(string(cfg.covariate_interactions), ' + ');
  end
  grp = cfg.random_effects.grouping; terms = string(cfg.random_effects.terms);
  re = join("(" + strjoin(terms,' + ') + "|" + grp + ")");
  formula = fx + " + " + re;
  opts.grouping = grp;
end
