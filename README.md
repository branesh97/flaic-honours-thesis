# Forecasting with Linear Augmented Information Combination (FLAIC)

Companion R scripts for the 2026 Monash Master of Applied Econometrics thesis by Branesh Prakash.

This repository implements regularised IComb augmentation with separately forecasted PCA components, together with conventional reconciliation, an adapted FLAP benchmark, empirical evaluation and supporting figures.

**Replication:** Read [README.txt](README.txt) for the input data layout, software dependencies, execution order, computation time and the limitations of the initial code release.

**Scripts:** Run stages 01, 03–09 in order from the repository root; stage 02 is sourced automatically. Full estimation takes substantial time (the FLAIC stage was approximately 36 hours using eight parallel workers).

**Supplementary results:** See [supplementary-mcb/](supplementary-mcb/) and [09_mcb_with_base.R](09_mcb_with_base.R) for the MCB comparisons including ETS base forecasts. Numerical CSVs must be regenerated from the saved evaluation objects; the accompanying text identifies findings reported in the submitted thesis.

**Data:** The TRA source CSV and large intermediate results are not provided in this initial code upload. Please check redistribution rights and file sizes before adding them.
