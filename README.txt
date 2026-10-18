FORECASTING WITH LINEAR AUGMENTED INFORMATION COMBINATION (FLAIC)
=================================================================

Companion R code for the 2026 Master of Applied Econometrics thesis by
Branesh Prakash, Monash University.

Repository: https://github.com/branesh97/flaic-honours-thesis

OVERVIEW
--------
FLAIC augments hierarchical base forecasts with separately generated ETS
forecasts of principal components. The application uses 112 monthly Australian
domestic tourism series over January 1998–December 2019 and compares ridge
(FLAIC-S), multivariate lasso (FLAIC-L), IComb, conventional reconciliation
methods and a hierarchy-adapted FLAP benchmark.

The experiment uses 144 rolling training windows of 120 months, horizons 1–12,
four PCA centring/scaling treatments, and eight augmentation dimensions:
p = 1, 2, 3, 5, 10, 20, 40, 77. Across two regularisers this gives 64 FLAIC
configurations. Lambda is selected using the final 40 training observations
as expanding-window internal validation; predictor and response variables
are standardised with an unpenalised intercept, and selection is based on
unscaled squared forecast errors.

REPOSITORY STRUCTURE AND EXECUTION
-----------------------------------
01_data_setup.R             Load CSV, form hierarchy, save settings
02_model_functions.R        Base forecasts, reconciliation and FLAIC functions
03_run_forecasts.R          Base/MinT/IComb, principal components and FLAIC
04_run_flap.R               Hierarchy-adapted FLAP benchmark
05_evaluate.R               Raw and scaled accuracy evaluation
06_scaled_error_prial.R     Scaled-error robustness tables and figures
07a_descriptive_analysis.R  Data/component descriptive analysis
07b_analysis_workbench.R   Diagnostic figures, MCB, predictor selection, FLAP
08_tables_figures.R        Raw-error PRIAL tables
09_mcb_with_base.R         Additional MCB comparisons against ETS base forecasts

data/                      TRA-final.csv. Used by Nguyen et al. (2025).
rds/                       Saved intermediary/final R objects and result tables
figures/                   Generated plots and visualisations
coefficients/              FLAIC coefficient matrices (144 per-origin RDS files) (uploaded later)
checkpoints/               Intermediate checkpoints (not automatically resumed) (uploaded later)
supplementary/             Supplementary MCB documentation and output CSVs
thesis/                    Thesis PDF (existing folder on GitHub)
thesis-qmd/                Thesis QMD files

Requirements: R; packages hts, forecast, glmnet, future, future.apply,
tidyverse, patchwork and zoo. Exact historical versions may be recovered from
run_info$session_info in rds/flaic_03_forecasts_complete.rds, if available.

To reproduce from a fresh local R session, set the working directory to
the repository ROOT, put TRA-final.csv in data/, and run:

 Rscript 01_data_setup.R
 Rscript 03_run_forecasts.R
 Rscript 04_run_flap.R
 Rscript 05_evaluate.R
 Rscript 06_scaled_error_prial.R
 Rscript 07a_descriptive_analysis.R
 Rscript 07b_analysis_workbench.R
 Rscript 08_tables_figures.R
 Rscript 09_mcb_with_base.R

02_model_functions.R is sourced by 03 and 04. The other scripts load their
saved inputs and can be rerun after their prerequisites have been generated.
07b requires the saved per-origin files in coefficients/.

COMPUTATION AND REPRODUCIBILITY
-------------------------------
The analysis used eight parallel R workers (01_data_setup.R: workers = 8).
The main expensive run too approximately 36 hours for the FLAIC estimation stage on
eight workers, not including all other processing stages. Runtime will differ
by machine, available memory, I/O and package versions. Forecast generation
saves periodic checkpoints, but 03_run_forecasts.R does NOT automatically
resume from them. Do not rerun the 36-hour computation just to generate
summary tables; the later scripts can reuse its completed RDS outputs.

Key saved outputs:
 rds/flaic_01_data_setup.rds
 rds/flaic_03_forecasts_complete.rds
 rds/flaic_04_flap_complete.rds
 rds/flaic_05_evaluation_complete.rds
 rds/table_prial_raw_errors.rds
 supplementary/MCB_with_base_all_comparisons.csv
 supplementary/MCB_with_base_exceptions.csv
 supplementary/MCB_with_base_report.txt

The data input is not bundled. The raw CSV must contain 264 months and 308
purpose-of-travel data columns after two leading columns, grouped in sets of
four columns per region.

SUPPLEMENTARY MCB COMPARISONS
-----------------------------
09_mcb_with_base.R recomputes the series-wise Friedman/Nemenyi comparisons
at horizons 1, 6 and 12, including ETS base forecasts in each panel alongside
matching IComb and eight FLAIC specifications. Including the base forecasts
changes the Nemenyi critical distance.

The submitted thesis reports that FLAIC specifications are significantly
better ranked than base forecasts except for two cases at horizon 6: ridge
FLAIC-S with p=77 under centred-unscaled and uncentred-unscaled PCA.
The reported findings are recorded in supplementary/MCB_with_base_summary_
from_thesis.txt; numerical CSVs are NOT included until recomputed using the
completed forecast-error RDS. A lack of statistical significance does not
mean FLAIC is significantly worse than base. The hierarchy series are related
by aggregation, so the ranking blocks are not independent.

SCOPE AND LICENSING
-------------------
The scripts were prepared for repository publication but have not been rerun
in a new environment. Some R functions were adapted from the published IComb
framework; cite the original literature in the accompanying thesis. Check
adapted source-license terms before choosing a public software license.

Large intermediate RDS objects, coefficients, raw data and checkpoints are
not included in this initial code upload. GitHub normally blocks files over
100 MB; use Git LFS or external archiving if needed. Preserve an unmodified
copy of the original thesis analysis scripts as a private backup.

Thesis citation: Prakash, B. (2026). Forecasting with Linear Augmented
Information Combination. Master of Applied Econometrics thesis,
Monash University.
