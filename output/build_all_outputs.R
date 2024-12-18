here::i_am("output/build_all_outputs.R")

# reads output/factorroder.rds
# writes output/monthly_climate_normals.rds
callr::rscript(here::here("output/1.site_climate_change.R"))

# reads output/phenf.rds
# writes output/censdf.rds
callr::rscript(here::here("output/2.censoring.R"))

# reads tmp/datlist.rds
# writes output/replication_points.rds
callr::rscript(here::here("output/3.replication_points.R"))

# reads output/factororder.rds
# writes output/tables/sitedat.rds
#
# reads inputs/forcing/typical_year_forc.csv
# writes output/figures/siteclimplot.png
#
# reads output/monthly_climate_normals.rds
# writes output/figures/monthly_climate_normals.png
#
# reads output/replication_points.rds and output/phenf.rds
# writes output/figures/MAT.png
#
# reads tmp/phendf.rds
# writes output/figures/sampling.png
#
# reads tmp/intercepts.rds
# writes output/tables/interceptsummary.rds
#
# reads tmp/slopes.rds
# writes output/tables/slopesummary.rds
#
# reads tmp/variation.rds
# writes output/figures/varoffsets.png
#
# reads tmp/indpredsummary.rds
# writes output/figures/independent_data_comp.png
#
# reads tmp/fretro_summary.rds
# writes output/figures/obsvsretro.png and output/figures/obsvsretroconceptual.png
#
# reads tmp/fpred_orch_avg_summary.rds and tmp/phendf.rds
# writes output/figures/genpred_gdd.png
#
# reads tmp/doy_annual_avg_pp_sum.rds and tmp/phendf.rds
# writes output/figures/orchpred_doy_male.png and output/figures/orchpred_doy_female.png
# and output/figures/orchpred_doy_vertical.png
#
# reads output/warmvscold.rds
# writes output/figures/provdiffdoy.png
#
# reads tmp/rank_correlation_wdist.rds
# writes output/figures/distrankcorr.png
#
# reads tmp/corr_model_results.rds
# writes output/tables/corr_model_table.rds
#
# reads and tmp/doy_typical_all_at_PGTIS.rds
# writes output/figures/away.png
#
# reads output/doy_normal_subset.rds
# writes output/figures/normal_predictions.png
callr::rscript(here::here("output/4.graphsandtables.R"))

# reads output/tables/sitedat.rds
# and tmp/siter.rds and output/phenf.rds
# writes output/sitesdvsprovmat_rsq_range.rds
callr::rscript(here::here("output/5.effectinvestigation_site.R"))

# reads inputs/latifoliaDistribution/shapefiles/latifolia_distribution_prj.shp
# writes output/figures/siteandparentmap.png
callr::rscript(here::here("output/6.map.R"))
