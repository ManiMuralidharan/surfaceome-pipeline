# Surfaceome Pipeline v3.0

A modular R Shiny research application for cell-surface (surfaceome) target discovery: CITE-seq ADT and mass-spectrometry differential analysis, cross-modality integration, transparent target ranking, out-of-sample model validation and evaluation, with publication-ready export.

**Author:** Muralidharan Mani, Ph.D. ([ORCID 0000-0002-1194-3616](https://orcid.org/0000-0002-1194-3616)), University of Wisconsin-Madison

## Quick start

Requires R >= 4.3 (recommended) and packages:

```r
install.packages(c('shiny','bslib','ggplot2','DT','pROC','PRROC','BiocManager'))
BiocManager::install('limma')
shiny::runApp('surfaceome-pipeline')   # path to the cloned folder
```

For a reproducible environment, record your local package versions with `renv::init()` and `renv::snapshot()` after installation. A prevalidated lockfile is NOT supplied because the package environment was not available during authoring.

## Input files

1. ADT feature-by-sample CSV with `feature` in column 1 and sample IDs as remaining headers. Counts must be nonnegative integers. `.rds` containing a matrix is also supported; native Seurat objects are deliberately not automatically coerced.
2. Mass spec matrix of feature-by-sample intensities, missing values permitted if two observations per group remain for each feature. Choose whether intensities are already log2 transformed.
3. Sample metadata CSV `sample_id,group` with labels `Normal` and `Tumor` (editable in the UI).
4. Surfaceome annotation CSV `feature,surface_evidence,normal_tissue_risk`, scores in [0,1]. These data require careful external curation and provenance; no CSPA/SURFY database is bundled or automatically fetched.
5. Evaluation predictions CSV `truth,scoreA,scoreB` (`scoreB` optional). Each score must be a predicted probability in [0,1]. Keep predictions strictly out of sample for unbiased estimation.
6. ML feature table CSV `sample_id,cohort,truth,<numeric features>`; `sample_id` optional. Prediction outputs can be sent to Evaluation automatically.

See `data/` for small **SIMULATED** example files; not actual tumor studies.

## Workflow

1. Upload counts and metadata, then process ADT (CLR or logCPM, both on the log2 scale so ADT and MS fold-changes are comparable). For per-cell CITE-seq, aggregate to pseudobulk per sample first.
2. Upload mass spec intensities and metadata, normalize and fit limma.
3. Integrate on exact shared IDs. **IDs must already be harmonized**; do not merge gene symbols with UniProt IDs without curated mapping.
4. Upload external surface evidence and normal-tissue expression risk annotations, then tune weights. Only features passing the FDR gate (default BH FDR < 0.05 in at least one modality) are ranked, and concordance credit is given only for up-regulation in both ADT and MS. Tumor-vs-normal logFC is a *proxy* for specificity, not a direct measure of disease exclusivity.
5. Run out-of-sample logistic regression using a stratified random hold-out or leave-one-cohort-out. Train-only predictor centering and scaling prevents straightforward test leakage; any feature selection must be specified a priori to avoid selection leakage.
6. Evaluate probabilities. Bootstrap uses stratified within-class resampling (percentile CI), with a *paired* DeLong comparison for model A vs B. You can change the classification threshold without changing AUROC/AUPRC. Calibration plot uses empirical quantile bins.
7. Export vector PDF, RDS analysis bundle, ranking CSV and methodology TXT.

## Limitations and honest boundaries

- The HTML V2 file was used as a **visual and navigation reference**. This is an independent, native-R Shiny reimplementation, not a React-to-Shiny byte-for-byte conversion.
- CITE-seq module processes ADT matrices (not RNA, barcodes, antibody-background denoising, clustering, doublet detection or Seurat-specific assays).
- Differential limma design handles **unpaired** Tumor vs Normal samples. Add patient, batch, and experimental blocking to the design before publication where applicable. Technical replicates must not be treated as independent samples. Counts are CLR-transformed for exploratory/working analysis, not a substitute for a validated CITE-seq-specific inferential model.
- Mass-spec normalization assumes comparability after log transform and median centering; missingness is not imputed. Missing-not-at-random proteomics requires more careful approaches.
- Target scores are transparent heuristic priorities and not trained/calibrated therapeutic-success probabilities. No actual CSPA/SURFY or Human Protein Atlas data are included.
- Logistic GLM is not penalized; high-dimensional predictors or perfect class separation can cause instability. Feature selection/hyperparameter tuning must occur within training folds, ideally with nested CV. LOCO with only a few cohorts has wide uncertainty.
- Reported stratified bootstrap intervals describe conditional uncertainty for given predictions and do not account for model refitting or across-cohort sampling variation.
- Simulated-data controls are prominent. RDS export retains objects and origin flags; do not present DEMO results as discoveries.
- Data uploaded into a Shiny server may be processed on that server. Handle patient-level information in a compliant environment.
- **Research only; not clinical diagnosis, target safety confirmation, or therapeutic recommendation.**

## Testing

```r
source('tests/test_core.R')
```

Tests use synthetic matrices. They require the listed packages. A full Shiny application smoke test should include launching the app and trying all uploads, navigation, and exports.

## Project organization

- `app.R` — sole UI/server composition point
- `R/utils_data.R` — validation, normalization and limma fit
- `R/utils_metrics.R` — AUROC, AUPRC, Brier, paired DeLong, bootstrap, confusion
- `R/utils_plotting.R` — consistent ggplot figures
- `R/mod_assays.R` — reusable ADT and MS assay UI/server
- `R/mod_integrate.R` — overlap, correlation and matched differential estimates
- `R/mod_prioritize.R` — external evidence weighting/ranking
- `R/mod_ml.R` — train-only GLM, held-out/LOCO predictions
- `R/mod_eval.R` — upload, curves, calibration, intervals, thresholds
- `R/mod_export.R` — downloads and provenance

## Citation

If you use this software, please cite it using the metadata in [`CITATION.cff`](CITATION.cff) (GitHub shows a "Cite this repository" button):

> Mani M. (2026). *Surfaceome Pipeline: an R Shiny application for CITE-seq and proteomics integration and cell-surface target prioritization* (Version 3.0.0) [Computer software]. https://github.com/ManiMuralidharan/surfaceome-pipeline

Please also cite limma (Ritchie et al., 2015, NAR), pROC (Robin et al., 2011, BMC Bioinformatics) and PRROC (Grau et al., 2015, Bioinformatics).

## License

MIT; see `LICENSE`.
