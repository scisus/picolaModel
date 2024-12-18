# This script fits a phenological model for lodgepole pine flowering

here::i_am("model/thermaltimemodel.R")

# depends #####
library(picolaDataFlowering)
library(dplyr)
library(brms)

source(here::here('model/phenology_functions.R'))

# data ####

## phenology
phendat <- picolaDataFlowering::picola_event %>%
  mutate(Tree = paste0(Orchard, Genotype, X, Y)) # create a unique Tree identifier since original data doesn't always have one

## forcing
dailyforc <- read.csv(here::here("inputs/forcing/dailyforc_1945_2012.csv"), header=TRUE, stringsAsFactors = FALSE)

# meta
spudat <- picolaDataFlowering::picola_SPUs 
prov_climate <- picolaDataClimate::picola_climatebc_parent_locs_normal_1961_1990 %>%
  rename(Genotype = id1, SPZ = id2)  %>%
  select(Genotype, SPZ, MAT, Latitude, Longitude) %>%
  mutate(Genotype = as.character(Genotype))

## data preparation for phenology model ####
phenf <- prepare_data(phendat, clim = dailyforc, spu = spudat) %>%
  left_join(prov_climate, relationship = "many-to-many") %>%
  filter(!is.na(MAT)) # Only keep genotypes with an associated source MAT - no breeding

saveRDS(phenf, file = here::here("output/phenf.rds"))

# create 4 datasets for 4 models

fbdat <- filter_sex_event(sex = "FEMALE", event = "begin", phenf)
fedat <- filter_sex_event(sex = "FEMALE", event = "end", phenf)
mbdat <- filter_sex_event(sex = "MALE", event = "begin", phenf)
medat <- filter_sex_event(sex = "MALE", event = "end", phenf)

saveRDS(list(fbdat = fbdat, fedat = fedat, mbdat = mbdat, medat = medat), file = here::here("tmp/datlist.rds"))


# model ####

# initialize parameter values with the right order of magnitude
initpars <- lapply(1:6, function(id) list(sigma = 30, Intercept = 300))

# model formula
bform <- brmsformula(sum_forcing | cens(censored, upper) ~ MAT + (1|Site) + (1|Genotype) + (1|Year) + (1|Tree))

# model prior
bprior <- c(prior("gamma(2.7, 0.012)", class = "Intercept", lb = 0),
            prior("normal(0,15)", class = "sigma"),
            prior("normal(0,9)", class = "sd"),
            prior("normal(0,25)", class = "b"))

# mcmc/computation settings
niter <- 4500
ncores <- 6
nchains <- 6

# function that returns future::multicore when running on unix and
# non-interactively, otherwise future::sequential

default_plan <- function() {
  if (.Platform$OS.type == "unix" && !interactive()) {
    return(future::multicore)
  } else {
    return(future::sequential)
  }
}

# This function fits multiple models, by default in parallel,
# using future::future_pmap and brms. Pass future::sequential
# for the plan parameter to fit models sequentially.
fit_thermaltimemodel <- function(param_frame, plan = default_plan()) {
  old_plan <- future::plan()
  future::plan(plan)
  message("Fitting ", nrow(param_frame), " models with seeds: ", paste(param_frame$seed, collapse = ", "))
  message("Using ", ncores, " cores, ", nchains, " chains, ", niter, " iterations")
  message("Priors: ", bprior)
  adapt_delta <- 0.95
  message("Adapt delta: ", adapt_delta)
  res <- furrr::future_pmap(
    param_frame,
    function(sex, event, seed) {

      data <- filter_sex_event(sex = sex, event = event, phenf)
      name_stem <- tolower(paste0(sex, "_", event))

      message("Fitting model for ", name_stem, " with ", nrow(data), " observations")
      flush.console()

      brm(bform, data = data,
        save_model = here::here(paste0("model/stan_output/", name_stem, ".stan")),
        file = here::here(paste0("model/stan_output/", name_stem, ".rds")),
        prior = bprior,
        init = initpars,
        iter = niter,
        cores = ncores,
        chains = nchains,
        sample_prior = TRUE,
        save_pars = save_pars(all = TRUE),
        file_refit = "always",
        seed = seed,
        control = list(adapt_delta = adapt_delta))
    },
    .options = furrr::furrr_options(seed=param_frame$seed[1], stdout=TRUE)
  )
  future::plan(old_plan)
  return(res)
}
