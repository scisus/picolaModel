# predict flowering events
# posterior predictive includes individual observation uncertainty
# grand means ignore group specific effects
# conditional effects include group specific effects as well as the uncertainty of fixed coefficients and the uncertainty of variance parameters for groups

here::i_am("analysis/6.predict.R")

# libraries
library(purrr)
library(dplyr)
library(tidybayes)
library(tidyr)

# data

modells <- readRDS(here::here("model/modells.rds"))
alldatls <- readRDS(here::here("tmp/datlist.rds"))
sitedat <- picolaDataClimate::picola_climatebc_locs_normal_1961_1990 %>% filter(id == "site")

labdf <- readRDS("tmp/labdf.rds")

n <- 6000 # when downsampling required

siteMAT <- sitedat %>%
  filter(id == "site") %>%
  select(Site, MAT, Elevation) %>%
  mutate(MAT = round(MAT, 1)) #DUPLICATED IN DOY TRANS

# orchards ############
## generic orchard ############
## at 100 points spanning full range of provenances in the data. Use fewer MAT points with length.out for smaller object.
neworchdat_avg <- expand.grid(MAT = seq(from = range(alldatls$fbdat$MAT)[1],
                                    to = range(alldatls$fbdat$MAT)[2], length.out = 100),
                          Year = "newyear",
                          Tree = "newtree",
                          Genotype = "newgenotype",
                          Site = "neworchard",
                          event = c("begin", "end"),
                          Sex = c("FEMALE", "MALE")) %>%
  split(list(.$event, .$Sex))

# posterior prediction #########
# for each site for the full range of provenances using an average site, year, genotype, and tree (using estimated gaussian prior to generate random effects). n draws, 95% HDPI ("generic predictions) #######

fpred_orch_avg <- purrr::map2(neworchdat_avg, modells,
                          .f = function(x,y) {add_predicted_draws(newdata = x,
                                                                  object = y,
                                                                  re_formula = NULL,
                                                                  allow_new_levels = TRUE,
                                                                  sample_new_levels = "gaussian",
                                                                  ndraws = n)}) %>%
  bind_rows()
saveRDS(fpred_orch_avg, file = here::here("tmp/fpred_orch_avg.rds"))

fpred_orch_avg_summary <- fpred_orch_avg %>%
  group_by(MAT, Year, Tree, Genotype, Site, event, Sex) %>%
  median_hdci(.prediction)
saveRDS(fpred_orch_avg_summary, here::here("tmp/fpred_orch_avg_summary.rds"))



# predict the global grand means: average predicted outcome ignoring group-specific deviations in intercept or slope

# grand mean ####
# build a dataframe with one entry for each dataset - just Sex and event, no groups, and calculate the average predicted outcome ignoring group specific deviations and individual level variation.


## expectation for trees sourced from all sites ####
fepred_allsites <- purrr::map(modells, function(x) {
  add_epred_draws(newdata = select(siteMAT, Site, MAT), object = x, re_formula = NA)}) %>%
  bind_rows(.id = "model") %>%
  select(-.chain, -.iteration) %>%
  left_join(labdf)
saveRDS(fepred_allsites, file = here::here("tmp/fepred_allsites.rds"))




