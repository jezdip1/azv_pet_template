# AZV-PET
Clean skeleton for PET LMEM modeling in MATLAB.

## Quickstart
1) Copy local paths config:
   ```bash
   cp config/paths.local.json.example config/paths.local.json
   ```
2) Put your training table at `data/raw/pet_table.csv`.
3) In VS Code run the task **train-pipeline**.
4) Artifacts appear in `models/` and reports in `reports/latest/`.

## Requirements
- MATLAB R2023a+ with Statistics and Machine Learning Toolbox
- Linux/macOS/Windows

## Structure
See the repository tree; MATLAB code under `matlab/` with package `+azvpet`.
