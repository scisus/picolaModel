# This script fits a phenological model for lodgepole pine flowering

here::i_am("data-raw/thermaltimemodel.R")

# depends #####
library(picolaDataFlowering)
library(dplyr)
library(brms)
#library(ggplot2)
#library(tidyr)
#library(tidybayes)
#library(forcats)
#library(ggbeeswarm)
#library(lubridate)

#theme_set(theme_dark())

source(here::here('data-raw/phenology_functions.R'))


# data ####

## phenology
phendat <- picolaDataFlowering::picola_event %>%
  mutate(Tree = paste0(Orchard, Genotype, X, Y)) # create a unique Tree identifier since original data doesn't always have one

## forcing
dailyforc <- read.csv(here::here("data-raw/forcing/dailyforc_1945_2012.csv"), header=TRUE, stringsAsFactors = FALSE)

# meta
spudat <- picolaDataFlowering::picola_SPUs # read.csv("../phd/data/OrchardInfo/LodgepoleSPUs.csv", header = TRUE, stringsAsFactors = FALSE)

prov_climate <- picolaDataClimate::picola_climatebc_parent_locs_normal_1961_1990 %>%
#read.csv("../lodgepole_climate/data/climateBC/climatebc_parent_locs_Normal_1961_1990Y_v730.csv") %>%
  rename(Genotype = id1, SPZ = id2)  %>%
  select(Genotype, SPZ, MAT) %>%
  mutate(Genotype = as.character(Genotype))

## data preparation for phenology model ####
phenf <- prepare_data(phendat, clim = dailyforc, spu = spudat) %>%
  left_join(prov_climate, relationship = "many-to-many") %>%
  filter(!is.na(MAT)) # Only keep genotypes with an associated source MAT - no breeding

saveRDS(phenf, file = here::here("data/phenf.rds"))

# create 4 datasets for 4 models

fbdat <- filter_sex_event(sex = "FEMALE", event = "begin", phenf)
fedat <- filter_sex_event(sex = "FEMALE", event = "end", phenf)
mbdat <- filter_sex_event(sex = "MALE", event = "begin", phenf)
medat <- filter_sex_event(sex = "MALE", event = "end", phenf)

saveRDS(list(fbdat = fbdat, fedat = fedat, mbdat = mbdat, medat = medat), file = here::here("data/datlist.rds"))


# model ####

# This model block is faster if you run models in parallel

# initialize parameter values with the right order of magnitude
initpars <- lapply(1:6, function(id) list(sigma = 30, Intercept = 300))

# model formula
bform <- brmsformula(sum_forcing | cens(censored, upper) ~ MAT + (1|Site) + (1|Genotype) + (1|Year) + (1|Tree))

# model prior
bprior <- c(prior("gamma(3.65, 0.01)", class = "Intercept", lb = 0),
            prior("normal(0,15)", class = "sigma"),
            prior("normal(0,9)", class = "sd"),
            prior("normal(0,25)", class = "b"))

# mcmc/computation settings
niter <- 3500
ncores <- 6
nchains <- 6

# female/receptivity begin
fbfit <- brm(bform, data = fbdat,
             save_model = here::here("data-raw/stan_output/female_begin.stan"),
             file = here::here("data-raw/stan_output/female_begin.rds"),
             prior = bprior,
             init = initpars,
             iter = niter,
             cores = ncores,
             chains = nchains,
             sample_prior = TRUE,
             save_pars = save_pars(all = TRUE),
             file_refit = "always",
             seed = 475902884)

# female/receptivity end
fefit <- brm(bform, data = fedat,
             save_model = here::here("data-raw/stan_output/female_end.stan"),
             file = here::here("data-raw/stan_output/female_end.rds"),
             prior = bprior,
             init = initpars,
             iter = niter,
             cores = ncores,
             chains = nchains,
             sample_prior = TRUE,
             save_pars = save_pars(all = TRUE),
             file_refit = "always",
             seed = 614687001)

# male/pollen shed begin
mbfit <- brm(bform, data = mbdat,
             save_model = here::here("data-raw/stan_output/male_begin.stan"),
             file = here::here("data-raw/stan_output/male_begin.rds"),
             prior = bprior,
             init = initpars,
             iter = niter,
             cores = ncores,
             chains = nchains,
             sample_prior = TRUE,
             save_pars = save_pars(all = TRUE),
             file_refit = "always",
             seed = 1960616843)

# male/pollen shed end
mefit <- brm(bform, data = medat,
             save_model = here::here("data-raw/stan_output/male_end.stan"),
             file = here::here("data-raw/stan_output/male_end.rds"),
             prior = bprior,
             init = initpars,
             iter = niter,
             cores = ncores,
             chains = nchains,
             sample_prior = TRUE,
             save_pars = save_pars(all = TRUE),
             file_refit = "always",
             seed = 1938696768)
