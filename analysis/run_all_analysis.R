here::i_am("analysis/run_all_analysis.R")

# reads model/stan_output/*.rds
# does basic diagnostic checks
callr::rscript(here::here("analysis/1.diagnostics.R"))

# reads tmp/datlist.rds and model/stan_output/*.rds
# writes model/modells.rds
# and tmp/fretro_summary.rds
# and output/fretrocomp.rds
#
# reads inputs/forcing/dailyforc_1945_2012.csv
# and output/phenf.rds
# writes output/dretrocomp.rds
callr::rscript(here::here("analysis/2.obsVSretro.R"))

# reads output/phenf.rds (from model/thermaltimemodel.R)
# writes output/factororder.rds
callr::rscript(here::here("analysis/3.factororder.R"))

# reads output/factororder.rds
# and model/modells.rds
# and output/phenf.rds
# writes tmp/slopes.rds and tmp/intercepts.rds
# and tmp/variation.rds and output/varsummary.rds
# and output/offsets_summary and siter.rds
# and output/genotyper.rds and output/yearr.rds
callr::rscript(here::here("analysis/4.modelparameters.R"))

# reads inputs/*/*.csv and model/modells.rds
# writes tmp/indpredsummary.rds
callr::rscript(here::here("analysis/5.model_check_independentdata.R"))

# reads tmp/datlist.rds and objects/modells.rds
# and tmp/labdf.rds
# writes tmp/fpred_orch_avg.rds
# and tmp/fpred_orch_avg_summary.rds
# and tmp/fepred_allsites.rds
callr::rscript(here::here("analysis/6.predict.R"))

# reads output/factororder.rds and inputs/forcing/*.csv
# and tmp/fepred_allsites.rds
# and tmp/intercepts.rds
# writes tmp/doy_typical_all_at_PGTIS.rds
# and tmp/doy_annual_pp_sum.rds
#
# reads tmp/fpred_orch_avg.rds
# writes tmp/doy_annual_avg_pp_sum.rds
# and tmp/doy_annual_exp_sum.rds
# and output/doy_normal_subset.rds
callr::rscript(here::here("analysis/7.dayofyear_translation.R"))

# reads tmp/doy_annual_pp_sum.rds
# and tmp/doy_annual_avg_pp_sum.rds
# writes output/max_ties.rds
# and tmp/rank_correlation_wdist.rds
# and tmp/corr_model_results.rds
# and output/warmvscold.rds
callr::rscript(here::here("analysis/8.doyanalysis.R"))

# reads output/phenf.rds
# and inputs/forcing/dailyforc_1945_2012.csv
# and tmp/doy_annual_pp_sum.rds
# and tmp/intercepts.rds
# and tmp/fpred_orch_avg.rds
# and inputs/Nilsson1981/swedish-pollen-timeseries.csv
callr::rscript(here::here("analysis/9.owens2005comp.R"))
