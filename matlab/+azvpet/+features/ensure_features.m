function T = ensure_features(T, cfg)
  namesT = string(T.Properties.VariableNames);
  if cfg.ensure.metadata_clean
    T = azvpet.features.ensure_metadata(T);
  end
  if cfg.ensure.SUL_LOG
    targetMask = endsWith(namesT, "_SUL_LOG");
    targets = namesT(targetMask);
    for nm = targets'
      base = extractBefore(nm, "_SUL_LOG");
      if ~ismember(nm, namesT)
        T = azvpet.features.ensureSUL_LOG(T, base);
      end
    end
  end
  if cfg.ensure.asymmetry
    T = azvpet.features.ensure_asymmetry(T);
  end
  if cfg.ensure.global_ref_pc1z
    T = azvpet.features.ensure_global_ref_pc1z(T);
  end
end
