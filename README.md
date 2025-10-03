# AZV PET Template

**Normative modelling of FDG-PET brain intensities** from a mixed clinical cohort (≈420 scans).  
The pipeline builds region-wise reference models (with age/sex & acquisition covariates), validates them with LOO cross-validation, calibrates predictions, and then scores new patients with prediction intervals and z-scores.

> ⚠️ Research prototype, **not** a medical device. Use under supervision of qualified clinicians.

---

## Why this exists (TL;DR)

Clinical PET data are messy: different ages, sexes, scanners, protocols… We normalize brain FDG-PETs into a common space, extract robust regional signals, and learn **expected intensity profiles** for healthy-looking brains—even when many scans come from oncology (ORL) or epilepsy. The mixed cohort still yields stable estimates thanks to sample size, robust aggregation (medians), and explicit modelling of covariates that soak up non-neurological variability.

---

## What the pipeline does

1. **Spatial normalization**  
   Each PET is affinely registered to the **MNI ICBM152** MRI template. We store fit quality (to gray-matter mask) via **MMI** and **NCC** metrics and keep them later as model covariates (`MMI_to_MNIGMS`, `NCC_to_MNIGMS`).

2. **Atlas parcellation**  
   Normalized PET is parcellated using **CerebrA** (`mni_icbm152_nlin_sym_09c_CerebrA_nifti`). Atlas parcels are **dilated** until touching neighbours to avoid “holes” and increase robustness. See the included label list (`CerebrA_LabelDetails.csv`) for region names and IDs.

3. **Intensity extraction (robust)**  
   Voxel values **< 5000** are set to **NaN** (hard threshold to suppress extra-cerebral & very low-count voxels). For every parcel we compute **median** and **mean** intensity and build one row per scan with DICOM-derived metadata (patient, injected dose, half-life, timing, scanner details …).

4. **Quantitative normalization (SUL) & variance stabilization (log)**  
   - **SUL** uses lean body mass to correct SUV bias in heavy/light patients:  
     \[
     \mathrm{SUL} = \frac{C_\text{tissue} \, [\mathrm{kBq/ml}] \cdot \mathrm{LBM}\,[\mathrm{kg}]}{D_\text{inj, decay-corrected to acquisition}\,[\mathrm{kBq}]}
     \]
     LBM is computed from patient sex/height/weight (configurable; e.g., James or Janmahasatian). Dose is decay-corrected using the radionuclide half-life and time offset from DICOM.  
   - We then use **`SUL_LOG = log(SUL)`** to reduce heteroskedasticity and improve linear model fit.

5. **Global reference (PCA)**
   From reference regions (brainstem & cerebellum; configurable) we build a subject-wise vector, run PCA across the cohort, and keep **PC1**, standardized to **`GlobalRefPC1_z`**. We persist the PCA projector to apply exactly the same transform to new scans.

6. **Statistical modelling (per region)**
   For each parcel’s **median SUL_LOG** we fit the linear model:
   ```
   AgeR1 + AgeR2 + AgeR3 + AgeR4 + Sex + BMI + lTime + lDose
 + logVoxelVol + logAcqDur_s + HasTOF + HasPSF + Subsets + FilterFWHM_mm
 + MMI_to_MNIGMS + NCC_to_MNIGMS + GlobalRefPC1_z
 + Sex:cAge + lDose:lTime
 + (1 | UNIS)
   ```
   **Definitions**
   - `AgeR1..AgeR4`: **restricted cubic spline** basis of **centered age** (`cAge = Age − mean(Age)`), 4 df (boundary knots ~ min/max, internal knots ~ quartiles).  
   - `Sex`: M/F.  
   - `BMI`: kg/m² (centered).  
   - `lTime`: log time from injection to acquisition start (min).  
   - `lDose`: log **decay-corrected** injected dose.  
   - `logVoxelVol`: log voxel volume.  
   - `logAcqDur_s`: log acquisition duration (s).  
   - `HasTOF`, `HasPSF`: scanner features flags.  
   - `Subsets`, `FilterFWHM_mm`: recon parameters.  
   - `MMI_to_MNIGMS`, `NCC_to_MNIGMS`: similarity of PET↔MNI gray-matter, used as image-quality covariates.  
   - `GlobalRefPC1_z`: PCA-based global reference covariate.  
   - `Sex:cAge`: age–sex interaction.  
   - `lDose:lTime`: dose–time interaction (kinetic confounding).  
   - `(1|UNIS)`: random intercept by **unit/scanner** (or another acquisition grouping), capturing between-unit offsets.

7. **Validation (LOO) & Prediction Intervals**
   - **Leave-One-Out CV**: refit the model N times leaving each case out; predict the held-out case.  
   - **Predictive SD**: we use LOO residuals to estimate dispersion for **prediction intervals (PI)**; CIs of the mean effect are based on model SE.  
   - **Calibration**: we decile-bin predictions and plot observed vs predicted means; an **identity line** indicates perfect calibration. If there’s a consistent slope or intercept bias, we learn a **monotone mapping** (isotonic/linear) and apply it when scoring new scans.  
   - **Outputs per region** include: mean effect vs age (by sex) with CI & PI, partial residuals, LOO Pred-vs-Obs, standardized residuals vs age, and a calibration plot.

---

## What the figures look like

> The examples below are for **Median-Hippocampus-Left-SUL_LOG** (place the PNGs under `docs/img/` to render in GitHub).

- **Age effect with CI & PI**  
  ![Age effect](docs/img/Median_Hippocampus_Left_SUL_LOG_age_effect.png)

- **Calibration (10 bins)**  
  ![Calibration](docs/img/Median_Hippocampus_Left_SUL_LOG_calibration.png)

- **Partial residuals (age term)**  
  ![Partial residuals](docs/img/Median_Hippocampus_Left_SUL_LOG_partial_resid.png)  
  ![Spline fit on partials](docs/img/Median_Hippocampus_Left_SUL_LOG_partial_resid_fitted.png)

- **Predicted vs Observed (LOO)**  
  ![Pred vs Obs](docs/img/Median_Hippocampus_Left_SUL_LOG_pred_vs_obs.png)

- **Standardized residuals vs age**  
  ![Z vs Age](docs/img/Median_Hippocampus_Left_SUL_LOG_zscore_vs_age.png)

---

## How well do regions model?

Most regions reach **R² ~0.96–0.99**; brainstem/cerebellum often exceed **0.99**, while ventricles/basal forebrain are typically lower (~0.94–0.95). For example (from `summary.json`).

---

## Using the template

### Requirements
- MATLAB R2022b+ (tested newer), with Statistics & Machine Learning Toolbox.
- (Optional) 3D Slicer for initial QC; ANTs/FSL if you re-run registrations.
- Input CSV assembled like the example **`test_pet_regions_with_mmi_JOINED_08092025.csv`** (parcel medians/means + DICOM metadata).

### Quick start

```matlab
% 1) Setup paths
>> startup

% 2) Train models on the cohort table
>> train_pipeline

% 3) Score a new patient
>> classify_from_bundle('data/to_test/patient_123.csv', 'reports/new_patients/')
```

---

## Repo layout

```
.
├── matlab/
├── data/
│   ├── raw/
│   ├── cohort/
│   └── to_test/
├── models/
├── reports/
├── docs/img/
└── CerebrA_LabelDetails.csv
```

---

## License
MIT

---
