function [T2, nameMap, opts2] = ensure_valid_varnames(T, opts)
  % Přejmenuje proměnné v T na validní názvy a vrátí mapu + upravené opts.
  % - opts.response_list (char/string/cellstr) se přemapuje na čisté názvy
  % - opts.grouping se přemapuje na čistý název
  
      % 1) vyrob clean names
      orig = string(T.Properties.VariableNames);
      cleaned = strings(size(orig));
      for i = 1:numel(orig)
          cleaned(i) = string(matlab.lang.makeValidName(orig(i), 'ReplacementStyle','underscore'));
      end
      cleaned = string(matlab.lang.makeUniqueStrings(cellstr(cleaned), {}, namelengthmax));
  
      % 2) přejmenuj tabulku
      T2 = T;
      T2.Properties.VariableNames = cellstr(cleaned);
  
      % 3) mapa
      nameMap = table(orig(:), cleaned(:), 'VariableNames', {'orig','clean'});
  
      % 4) opts → opts2 (přemapuj responses & grouping)
      opts2 = opts;
      if isfield(opts,'response_list')
          rs = opts.response_list;
          if ischar(rs), rs = string({rs}); end
          if iscellstr(rs), rs = string(rs(:)); end
          if isstring(rs),  rs = rs(:); else, rs = string(rs); end
  
          cleanResp = strings(size(rs));
          for i = 1:numel(rs)
              j = find(nameMap.orig == string(rs(i)), 1);
              if ~isempty(j)
                  cleanResp(i) = nameMap.clean(j);
              else
                  % když už je to validní/odpovídá už čistému názvu
                  cleanResp(i) = string(matlab.lang.makeValidName(char(rs(i)), 'ReplacementStyle','underscore'));
              end
          end
          opts2.response_list = cellstr(cleanResp);
      end
  
      if isfield(opts,'grouping')
          g = string(opts.grouping);
          j = find(nameMap.orig == g, 1);
          if ~isempty(j)
              opts2.grouping = char(nameMap.clean(j));
          else
              opts2.grouping = char(matlab.lang.makeValidName(char(g),'ReplacementStyle','underscore'));
          end
      end
  end
  