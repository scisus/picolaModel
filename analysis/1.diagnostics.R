# ESS and Rhatdiagnostics for stan models
#

here::i_am("analysis/1.diagnostics.R")

library(dplyr)
library(purrr)

fbfit <- readRDS(here::here('model/stan_output/female_begin.rds'))
fefit <- readRDS(here::here('model/stan_output/female_end.rds'))
mbfit <- readRDS(here::here('model/stan_output/male_begin.rds'))
mefit <- readRDS(here::here('model/stan_output/male_end.rds'))

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

fits <- list(fbsims = fbfit, fesims = fefit, mbsims = mbfit, mesims = mefit)
for (i in 1:length(fits)) {
  name <- names(fits)[i]
  fit <- fits[[i]]$fit
  seed <- rstan::get_seed(fit)
  message("------------------------------------")
  message("Checking diagnostics for ", name, " run with seed ", seed)
  print(rstan::check_hmc_diagnostics(fit))
  if (rstan::get_num_divergent(fit) > 0) {
    message("Divergent transitions detected in ", name, " run with seed ", seed)
    stop("Ending simulation early")
  }
}

# female begin is the hardest model to fit - code below used to help diagnose problems during modeling
# library(shinystan)
# launch_shinystan(fbfit)

