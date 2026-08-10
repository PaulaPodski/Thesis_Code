# Machine Learning Classifiers to Predict Prostate Cancer Mortality from Whole Genome Bisulfite Sequencing Data 

# This repository contains the machine learning pipeline used to evaluate whether DNA methylation from whole genome bisulfite sequencing (WGBS) can predict long-term prostate cancer-specific mortality in intermediate-risk (Gleason 6/7) patients, using a 72-sample subset of the Trans-Atlantic Prostate Group (TAPG) cohort.

# Project overview
# Cohort: 72 intermediate-risk (Gleason 6/7) patients from TAPG, PBAT-WGBS, ~10 year follow-up
# Goal: identify differentially methylated CpG sites associated with PCa-specific mortality and evaluate elastic net classifiers on held-out test data
#Key methodological point: train/test partitions are processed independently through the full methylKit workflow (feature selection, coverage filtering, imputation) to avoid data leakage from an earlier full-cohort version of this pipeline
# Pipeline stages
# Stratified train/test split (70/30 and 80/20), set.seed(42)
# Independent methylKit processing per split (methRead → filterByCoverage → normalizeCoverage → unite, no reorganize())
# Differential methylation analysis on training data only: Top 20 CpG pool selection
# Missing data imputation (mean vs. KNN, k=10) on training and test beta-matrices separately
# Elastic net regression, alpha tuned per split/imputation method via CV deviance
# Manual prediction via matrix multiplication (no refitting on test data)
# Model evaluation (AUC, confusion matrices, ROC) and SHAP-based feature interpretation

## Environment Setup

This pipeline was run on QMUL's Apocrita HPC cluster. To reproduce the environment:

```bash
module load miniforge/25.3.0
conda env create -f environment.yml
conda activate methylation_env_clean
```

If you already have a conda/miniforge installation and just want to activate an existing, matching environment:

```bash
source /path/to/your/miniforge3/etc/profile.d/conda.sh
conda activate methylation_env_clean
```

All scripts require this environment to be active before running any scripts.
