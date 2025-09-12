function T = ensure_features(T, cfg)
  namesT = string(T.Properties.VariableNames);

  if cfg.ensure.metadata_clean
    T = azvpet.features.ensure_metadata(T);
    namesT = string(T.Properties.VariableNames);
  end

  % --- SUL_LOG z bazálních 'Mean_*' (pokud chybí) ---
  if isfield(cfg,'ensure') && isfield(cfg.ensure,'SUL_LOG') && cfg.ensure.SUL_LOG
    baseCols = namesT(startsWith(namesT,"Mean_") & ~endsWith(namesT,"_SUL_LOG"));
    for base = baseCols'
      dst = base + "_SUL_LOG";
      if ~ismember(dst, namesT)
        T = azvpet.features.ensureSUL_LOG(T, base);
      end
    end
    namesT = string(T.Properties.VariableNames);
  end

  if isfield(cfg.ensure,'asymmetry') && cfg.ensure.asymmetry
    T = azvpet.features.ensure_asymmetry(T);
  end

  if isfield(cfg.ensure,'global_ref_pc1z') && cfg.ensure.global_ref_pc1z
    T = azvpet.features.ensure_global_ref_pc1z(T);
  end
end
