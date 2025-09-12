function T = ensure_features(T, cfg)
  namesT = string(T.Properties.VariableNames);
  if cfg.ensure.metadata_clean
    T = azvpet.features.ensure_metadata(T);
    namesT = string(T.Properties.VariableNames);
  end
  if cfg.ensure.SUL_LOG
    for nm = namesT'
      if startsWith(nm, "Mean_") && ~endsWith(nm, "_SUL_LOG") && ~endsWith(nm, "_Left") && ~endsWith(nm, "_Right")
        T = azvpet.features.ensureSUL_LOG(T, nm);
      end
    end
    namesT = string(T.Properties.VariableNames);
  end
  if cfg.ensure.asymmetry
    T = azvpet.features.ensure_asymmetry(T);
  end
  if cfg.ensure.global_ref_pc1z
    T = azvpet.features.ensure_global_ref_pc1z(T);
  end
end
