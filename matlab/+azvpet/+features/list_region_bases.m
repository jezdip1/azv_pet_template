function baseNames = list_region_bases(T, varargin)
  % Vrátí cellstr base názvů regionů podle preferovaného prefixu.
  %   baseNames = azvpet.features.list_region_bases(T, 'Prefix',"Median")
  %
  % Volby:
  %   'Prefix'    : "Median" (default) nebo "Mean"
  %   'IncludeRx' : regulární výraz pro whitelist (volit.)
  %   'ExcludeRx' : regulární výraz pro blacklist (volit.)
  %
  % Pozn.:
  %  - Vrací názvy BEZ sufixů "_SUL" / "_SUL_LOG".
  %  - Když v T není raw sloupec, ale je jen *_SUL / *_SUL_LOG, vrátí správný base.
  %  - Duplikáty odstraní (set).
  
      p = inputParser;
      addParameter(p, 'Prefix', "Median");
      addParameter(p, 'IncludeRx', "");
      addParameter(p, 'ExcludeRx', "");
      parse(p, varargin{:});
      pref   = string(p.Results.Prefix);
      inclRx = string(p.Results.IncludeRx);
      exclRx = string(p.Results.ExcludeRx);
  
      if ~ismember(pref, ["Median","Mean"])
          error('list_region_bases:Prefix','Prefix musí být "Median" nebo "Mean".');
      end
  
      v = string(T.Properties.VariableNames);
  
      % Kandidáti s daným prefixem
      cand = v(startsWith(v, pref + "_"));
  
      % Odstripuj případné SUL sufixy na base
      % 1) *_SUL_LOG → base
      isLog = endsWith(cand, "_SUL_LOG");
      cand(isLog) = erase(cand(isLog), "_SUL_LOG");
      % 2) *_SUL → base
      isSul = endsWith(cand, "_SUL");
      cand(isSul) = erase(cand(isSul), "_SUL");
  
      % Teď jsou všechny položky formálně "base" (včetně těch, které v T nemusí mít raw variantu)
  
      % Volitelný whitelist / blacklist
      if strlength(inclRx) > 0
          cand = cand(~cellfun('isempty', regexp(cellstr(cand), inclRx, 'once')));
      end
      if strlength(exclRx) > 0
          cand = cand(cellfun('isempty', regexp(cellstr(cand), exclRx, 'once')));
      end
  
      % Zbav se zjevných ne-regionálních věcí, kdyby náhodou začínaly stejně:
      % (ponecháme jednoduchý filtr: nepouštěj metriky s jasnými technickými názvy)
      dropRx = "^(Median|Mean)_(Voxel|Noise|SUV|SUL|Global|Mask|QC|Debug|Test|Tmp|AI_)";
      keep = cellfun('isempty', regexp(cellstr(cand), dropRx, 'once'));
      cand = cand(keep);
  
      % Set/unique + vrátit jako cellstr
      baseNames = unique(cellstr(cand), 'stable');
  end
  