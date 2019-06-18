# coral-bleach-recovery
This repository contains R code for the Bayesian models accompanying Robinson, Wilson and Graham "**Abiotic and biotic controls on coral recovery 16 years after mass bleaching**", Coral Reefs 2019.

[https://link.springer.com/article/10.1007/s00338-019-01831-7](https://link.springer.com/article/10.1007/s00338-019-01831-7)

The following R packages were used to analyse data.

```
install.packages(c("tidyverse", "rethinking", "here"))
```

**data**

* `recovery_year_model.Rdata` contains model structure (*rec.year.m*) and parameter effect sizes (*rec.params*) for predicting recovery year
* `recovery_trajectory.Rdata` contains logistic model strucutre (*rec.trajectory.m*), predicted recovery years (*base*) and recovery trajectory for each reef over 100 years (*rec.trajectory*)
* `posterior_sims.Rdata` contains posterior samples for each recovery year predictor covariate (*depth, complex, init_cover, wave, herb, coral_juv, nitrogen*)
* `jacknife/` contains model structures and posterior samples for jacknife sensitivity analysis

**figures**

* `logistic_models_predictions.pdf` are model predictions for each candidate logistic model, with uncertainty intervals and plotted against observed data
* `recovery_diagnostic_jackknife.pdf` are parameter effect sizes for each jacknife subsample

**scripts**

* `1_logistic_growth_model.R` fits Bayesian logistic growth models
* `2_recovery_year_model.R` fits Bayesian linear model to predict recovery year
* `3_jacknife_analysis.R` runs jacknife sensitivity analysis on recovery year model (2)
* `scaling_function.R` is generic function for scaling covariates to mean of 0 (for continuous) or creating dummy variables (for categorical)
