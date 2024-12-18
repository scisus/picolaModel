## Analysis

Scripts that analyze the model results and/or write out objects to `tmp` and `output` used in later scripts.

- `1.diagnositics.R` calculate model diagnostics like Rhat and ESS
- `2.obsVSretro.R` compare observations to retrodictions
- `3.factororder.R` order factors for making good graphs
- `4.modelparameters.R` extract parameter values from thermal time model
- `5.model_check_independentdata.R` check model predictions against data in O'Reilly 1988 and Nilsson 1981
- `6.predict.R` predict thermal time for events from models 
- `7.dayofyear_translation.R` translate thermal time predictions into day of year 
- `8.doyanalysis.R` investigate patterns of flowering across sites and between provs (rank correlations, provenance differences, and variation among sites)
- `9.owens2005comp.R` calculate GDD for dates reported in Owens 2005

- `run_all_analysis.R` run all the numbered scripts
- `prediction_examples.*` examples of making predictions using the model
