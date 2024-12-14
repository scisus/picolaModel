# ESS and Rhatdiagnostics for stan models
#

here::i_am("data-raw/diagnostics.R")

library(dplyr)
library(purrr)
#library(rstan)

fbfit <- readRDS(here::here('data-raw/stan_output/female_begin.rds'))
fefit <- readRDS(here::here('data-raw/stan_output/female_end.rds'))
mbfit <- readRDS(here::here('data-raw/stan_output/male_begin.rds'))
mefit <- readRDS(here::here('data-raw/stan_output/male_end.rds'))

# format for Rhat and ESS calculations

sims <- list(fbsims = as.array(fbfit), fesims = as.array(fefit), mbsims = as.array(mbfit), mesims = as.array(mefit))

# calculate minimum effective sample size (bulk and tail) for a two-dimensional array whose rows are equal to the number of iterations of the Markov Chain(s) and whose columns are equal to the number of Markov Chains (preferably more than one).
miness <- function(sims) {
  bulk_ess <- apply(sims, MARGIN = 3, FUN = rstan::ess_bulk)
  tail_ess <- apply(sims, MARGIN = 3, FUN = rstan::ess_tail)

  bulk <- min(bulk_ess)
  tail <- min(tail_ess)

  return(data.frame(bulk_ess = bulk, tail_ess = tail))
}


minesses <- purrr::map(sims, miness) %>%
  dplyr::bind_rows(.id = "id")
  # Goal is for ESS to be at least 100 for each chain. e.g. if 6 chains, you're in trouble if min ESS < 600.
print(minesses)
stopifnot(all(minesses$bulk_ess > 600, minesses$tail_ess > 600))

maxrhat <- function(sims) {
  rhats <- apply(sims, MARGIN = 3, FUN = rstan::Rhat)
  maxhat <- max(rhats)
  return(maxhat)
}

# Should be less than 1.01
rhats <- purrr::map(sims, maxrhat)
print(rhats)
stopifnot(all(rhats < 1.01))

# female begin is the hardest model to fit - code below used to help diagnose problems during modeling
# fbpars <- data.frame(rstan::extract(fbfit))


# pairs(fbfit, variable = c("mu", "sigma", "sigma_cens", "mu_site", "sigma_site", "mu_year", "sigma_year", "mu_genotype", "sigma_genotype", "mu_prov", "sigma_prov"))
# pairs(fbfit, pars = c("mu", "sigma", "sigma_site", "sigma_year", "sigma_genotype", "sigma_prov"))
# pairs(fbfit, pars = c("mu", "sigma", "mu_site", "mu_year", "mu_genotype", "mu_prov"))
# nuts <- bayesplot::nuts_params(fbfit)
# draws <- as.array(fbfit)
# bayesplot::mcmc_parcoord(draws, pars = vars("mu", "sigma", starts_with("mu_"), starts_with("sigma_"), contains("alpha_site")), np=nuts, transform = function(x) {(x - mean(x)) / sd(x)})
# bayesplot::mcmc_intervals(draws, pars = vars(starts_with("mu_"), starts_with("sigma_"), contains("delta_site"), "sigma"))
# bayesplot::mcmc_intervals(draws, pars = vars( contains("year")))
# bayesplot::mcmc_areas(draws, pars = vars(contains("sigma_")))
# bayesplot::mcmc_areas(draws, pars = vars(contains("mu_")))
#female_begin_censored <- fit_model(phendat = censored, sex = "FEMALE", event = "begin")
# library(shinystan)
# launch_shinystan(fbfit)

