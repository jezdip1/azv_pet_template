function [T, info] = ensure_metadata_vars(T, opts)
  %ENSURE_METADATA_VARS Sjednotí názvy a jednotky metadat potřebných pro SUL.
  %  - Vytvoří/naplní: InjectedDose_MBq, PatientHeight_cm, Sex, LBM_kg
  %  - Nevrací chybu, pokud některé zdroje chybí; vrací 'info' se stavem.
  %
  %  Použití:
  %    [T, info] = azvpet.features.ensure_metadata_vars(T);
  %
  %  Pozn.: LBM počítá defaultně dle 'James'. Přepni opts.lbmVersion = 'Janmahasatian'
  %         pokud chceš jiný model.
  
      if nargin < 2, opts = struct(); end
      if ~isfield(opts,'lbmVersion'), opts.lbmVersion = 'James'; end
  
      vn = T.Properties.VariableNames;
      info = struct('var','', 'status','', 'message','');
  
      function addInfo(v,s,m)
          info(end+1) = struct('var',v,'status',s,'message',m); %#ok<AGROW>
      end
      info(1)=[]; % drop prealloc
  
      % ----- InjectedDose_MBq -----
      if ismember('InjectedDose_MBq', vn)
          addInfo('InjectedDose_MBq','exists','');
      else
          if ismember('RadionuclideTotalDose', vn)
              % DICOM RadionuclideTotalDose je v Bq -> převedeme na MBq
              T.InjectedDose_MBq = T.RadionuclideTotalDose * 1e-6;
              addInfo('InjectedDose_MBq','created','Z RadionuclideTotalDose (Bq) -> MBq.');
          elseif ismember('Dose_MBq', vn)
              T.InjectedDose_MBq = T.Dose_MBq;
              addInfo('InjectedDose_MBq','alias','Z Dose_MBq.');
          elseif ismember('InjectedDose_Bq', vn)
              T.InjectedDose_MBq = T.InjectedDose_Bq * 1e-6;
              addInfo('InjectedDose_MBq','created','Z InjectedDose_Bq (Bq) -> MBq.');
          else
              addInfo('InjectedDose_MBq','missing','Nenalezen žádný z [RadionuclideTotalDose, Dose_MBq, InjectedDose_Bq].');
          end
      end
  
      vn = T.Properties.VariableNames;
  
      % ----- PatientHeight_cm -----
      if ismember('PatientHeight_cm', vn)
          addInfo('PatientHeight_cm','exists','');
      else
          if ismember('PatientSize_m', vn)
              T.PatientHeight_cm = T.PatientSize_m * 100;
              addInfo('PatientHeight_cm','created','Z PatientSize_m (m) -> cm.');
          elseif ismember('PatientHeight_mm', vn)
              T.PatientHeight_cm = T.PatientHeight_mm / 10;
              addInfo('PatientHeight_cm','created','Z PatientHeight_mm -> cm.');
          else
              addInfo('PatientHeight_cm','missing','Nenalezen žádný z [PatientSize_m, PatientHeight_mm].');
          end
      end
  
      vn = T.Properties.VariableNames;
  
      % ----- Sex (M/F) -----
      if ismember('Sex', vn)
          addInfo('Sex','exists','');
      else
          if ismember('PatientSex', vn)
              sx = string(T.PatientSex);
              sx = upper(sx);
              sx = replace(sx, "MALE", "M");
              sx = replace(sx, "FEMALE", "F");
              T.Sex = sx;
              addInfo('Sex','created','Z PatientSex.');
          else
              addInfo('Sex','missing','Nenalezen PatientSex.');
          end
      end
  
      vn = T.Properties.VariableNames;
  
      % ----- LBM_kg -----
      if ismember('LBM_kg', vn)
          addInfo('LBM_kg','exists','');
      else
          need = ismember({'PatientWeight_kg','PatientHeight_cm','Sex'}, vn);
          if all(need)
              wt = T.PatientWeight_kg;
              ht = T.PatientHeight_cm;
              sx = upper(string(T.Sex));
              isM = (sx=="M"); isF = (sx=="F");
              L = nan(height(T),1);
  
              switch lower(opts.lbmVersion)
                  case 'janmahasatian'
                      BMI = wt ./ ((ht/100).^2);
                      L(isM) = 9270*wt(isM) ./ (6680 + 216*BMI(isM));
                      L(isF) = 9270*wt(isF) ./ (8780 + 244*BMI(isF));
                  otherwise % James
                      ratio = wt ./ ht; % kg/cm
                      L(isM) = 1.10*wt(isM) - 128*(ratio(isM).^2);
                      L(isF) = 1.07*wt(isF) - 148*(ratio(isF).^2);
              end
              T.LBM_kg = L;
              addInfo('LBM_kg','created',['LBM z ', opts.lbmVersion, '.']);
          else
              addInfo('LBM_kg','missing','Chybí některé z [PatientWeight_kg, PatientHeight_cm, Sex].');
          end
      end
  end
  