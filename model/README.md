## Model

- `conceptualanalysis.*` Cribs from the first few steps of [Betancourt's workflow](https://betanalpha.github.io/assets/case_studies/principled_bayesian_workflow.html), outlining the problem and determining what domain specific knowledge can be brought to bear on the priors
- `phenology_functions.R`: Helper functions for `thermaltimemodel.R` and analysis scripts
- `thermaltimemodel.R` Thermal time models of flowering events in Stan
  - Generated Stan code in `stan_output/[sex]_[event].stan`
  - Model output in `stan_output/[sex]_[event].rds`
- `fit_model.R`: run models from the command line with fixed seeds used in publication, new random seeds, or simply download existing model fit from zenodo. Permits multiple runs of each model.
