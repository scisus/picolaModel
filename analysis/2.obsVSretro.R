# Compare retrodictions to observations
#
## this script determines the proportion of forcing retrodictions that match observations. Observations are ranges that the event occurred in and retrodictions are from models with uncertainty intervals

# Note: this script is quite slow and resource intensive.

## Are the mean begin, end, and length retrodictions within the expected ranges? The exact flowering period is never observed because of censoring.
# Using the first and last observed flowering forcing/doy we construct the max begin and minimum end forcing/doy and minimum flowering period length in forcing and days for the data.
# Using the last observed before flowering forcing/doy and the first observed after flowering forcing/doy, we construct a minimum begin day, maximum end day, and maximum flowering period length for the data. We expect model estimates to be between the min and max ranges for the observations.

here::i_am("analysis/2.obsVSretro.R")

# libraries ####
library(dplyr)
library(purrr)
library(tidybayes)
library(tidyr)
library(ggplot2)

source(here::here('model/phenology_functions.R'))

# functions ####
# calculate the proportion of trues in a boolean vector x
calc_prop_in_int <- function(x) {sum(x)/length(x)}

# calculate the length of the phenological period in forcing units or days, depending on the df you feed in (fsim or dsim in this script). df must have columns dat_min_begin, dat_max_begin, dat_min_end, dat_max_end, retro_mean_end, retro_mean_begin, retro_sd_begin, and retro_sd
calc_len <- function(df) {
  lendf <- df %>%
    dplyr::select(-retro_min, -retro_max) %>%
    tidyr::pivot_wider(id_cols = c(Index, Sex), names_from = event, values_from = contains("_")) %>%
    dplyr::mutate(dat_range_min = dat_min_end - dat_max_begin,
                  dat_range_max = dat_max_end - dat_min_begin,
                  retro_len_mean = retro_mean_end - retro_mean_begin,
                  retro_len_sd = sqrt(retro_sd_end^2 + retro_sd_begin^2),
                  retro_len_min = retro_len_mean - retro_len_sd,
                  retro_len_max = retro_len_mean + retro_len_sd) %>%
    dplyr::group_by(Index, Sex) %>%
    dplyr::mutate(inint_mean_len = dplyr::between(x = retro_len_mean, left = dat_range_min, right = dat_range_max),
                  inint_onesd_len = any(findInterval(c(dat_range_min, dat_range_max), c(retro_len_min, retro_len_max))))

  return(lendf)
}

# build a table with the proportion of retrodiction means and 1 sigma retrodiction intervals that are contained by or overlap the event range
comp_retro2dat <- function(simdf, lendf) {
  retrocomp <- simdf %>%
    group_by(Index, Sex, event) %>%
    # are begin and end mean retrodictions in the obs interval
    mutate(inint_mean = between(x = retro_mean,
                                left = dat_min, right = dat_max)) %>%
    # do begin and end 1 sigma interval retrodictions overlap the obs interval
    mutate(inint_onesd = any(findInterval(c(dat_min,dat_max),
                                          c(retro_min,retro_max)))) %>%
    # join with length interval tests
    select(Index, Sex, event, contains("inint")) %>%
    pivot_wider(id_cols = c(Index, Sex), names_from = event, values_from = contains("inint")) %>%
    full_join(select(lendf, Index, Sex, inint_mean_len, inint_onesd_len)) %>%
    ungroup() %>%
    group_by(Sex) %>%
    summarise_at(vars(contains("inint")), calc_prop_in_int) %>%
    select(Sex, contains("begin"), contains("end"), contains("len"))
  colnames(retrocomp) <- c("Sex", "Begin", "Begin_sd", "End", "End_sd", "Length", "Length_sd")

  return(retrocomp)
}

# for day observations, what are the min or max doy for begin and end events
minmaxdat <- function(df, minormax) {
  if (minormax == "min") {
    events <- c(1,3)
    cols <- c("dat_min")
  }
  if (minormax == "max") {
    events <- c(2,4)
    cols <- c("dat_max")
  }

  stopifnot("argument minormax must be \"min\" or \"max\"" = minormax %in% c("min", "max"))

  mmdf <- df %>%
    filter(Event_Obs %in% events) %>%
    mutate(event = case_when(Event_Obs == events[1] ~ "begin",
                             Event_Obs == events[2] ~ "end")) %>%
    select(Index, Sex, event, DoY) %>%
    group_by(Index, Sex)

  colnames(mmdf)[4] <- cols

  return(mmdf)
}



## data and models ##############

alldatls <- readRDS(here::here("tmp/datlist.rds"))

modells <- list(fb = readRDS(here::here("model/stan_output/female_begin.rds")),
                fe = readRDS(here::here("model/stan_output/female_end.rds")),
                mb = readRDS(here::here("model/stan_output/male_begin.rds")),
                me = readRDS(here::here("model/stan_output/male_end.rds")))
saveRDS(modells, here::here("model/modells.rds")) #1.5GB

## modeled forcing ####
# simulate new forcing observations from the model. this is a slow step. I'm using the full model to make retrodictions, not subsampling. Makes an 12.8GB object without downsampling
fretro <- purrr::map2(alldatls, modells, function(x,y) {add_predicted_draws(newdata = x, object = y)}) %>%
  bind_rows()

fretro_summary <- fretro %>%
  group_by(Index, Year, Sex, Site, Orchard, Genotype, Tree, event, sum_forcing, censored) %>%
  median_hdci(.prediction)
saveRDS(fretro_summary, here::here("tmp/fretro_summary.rds"))

# begin & end retrodictions + data
fsim <- fretro %>%
  # summarise by observation
  group_by(Index, Sex, event, censored, sum_forcing, upper) %>%
  summarise(retro_mean = mean(.prediction), retro_sd = sd(.prediction)) %>%
  # calculate min and max of a 1 sigma interval around the mean
  mutate(retro_min = retro_mean - retro_sd, retro_max = retro_mean + retro_sd) %>%
  # identify min and max of range for observations
  mutate(dat_min = case_when(censored == "left" ~ 0,
                             censored == "right" ~ sum_forcing,
                             censored == "interval" ~ sum_forcing),
         dat_max = case_when(censored == "left" ~ sum_forcing,
                             censored == "right" ~ Inf,
                             censored == "interval" ~ upper)) %>%
  ungroup() %>%
  select(-sum_forcing, -upper, -censored)

  # length calculations
flen <- calc_len(fsim)

# fretrocomp is a table that describes how many model forcing estimates are within the observed forcing ranges
fretrocomp <- comp_retro2dat(fsim, flen)
saveRDS(fretrocomp, here::here("output/fretrocomp.rds"))

# determine the proportion of doy retrodictions that match observations ####


### historical climate data ###
dailyforc <- read.csv(here::here("inputs/forcing/dailyforc_1945_2012.csv"))# site clim with forcing
# might need to drop trench and border site

## doy data ####
phenf <- readRDS(here::here("output/phenf.rds"))

## convert forcing retrodictions to doy ####
## downsample to 1000 draws for each obs to keep sizes manageable
samp <- sample(1:max(fretro$.draw), size = 1000)
dretro <- forcing_to_doy(dailyforc, filter(fretro, .draw %in% samp), aforce = "sum_forcing", bforce = ".prediction", newdoycolname = "retro_doy") %>%
  ungroup() %>%
  select(Index, Sex, event, retro_doy)

## dat v model ####

# calculate event ranges from data and add length
ddat <- full_join(minmaxdat(phenf, minormax = "min"),
                  minmaxdat(phenf, minormax = "max"))

# calculate model estimate and 1 sigma ranges
dsim <- dretro %>%
  #summarise by observation
  group_by(Index, Sex, event) %>%
  summarise(retro_mean = mean(retro_doy), retro_sd = sd(retro_doy)) %>%
  ungroup() %>%
  #calculate min and max of a one sigma interval around the mean
  mutate(retro_min = retro_mean - retro_sd, retro_max = retro_mean + retro_sd) %>%
  left_join(ddat) %>%
  mutate(dat_min = coalesce(dat_min, 0),
         dat_max = coalesce(dat_max, Inf))


# calculate length
dlen <- calc_len(dsim)

# calculate proportion of retrodictions in the data ranges
dretrocomp <- comp_retro2dat(dsim, dlen)
saveRDS(dretrocomp, here::here("output/dretrocomp.rds"))

# Examine residuals ####
## randomized quantile residuals ####
## technically, I think these are only acceptable for the interval censored data and not the end censored data, but I think I can honestly transform end to interval with sufficiently wide interval estimates. though it might be more honest here to use prior info for extremes
fq <- fretro %>%
  mutate(dat_min = case_when(censored == "left" ~ 0,
                             censored == "right" ~ sum_forcing,
                             censored == "interval" ~ sum_forcing),
         dat_max = case_when(censored == "left" ~ sum_forcing,
                             censored == "right" ~ 800,
                             censored == "interval" ~ upper))

# Calculate randomized quantile residuals
# https://mjskay.github.io/tidybayes/articles/tidybayes-residuals.html
# Dunn & Smyth 1996
fq %>%
  summarise(
    p_lower = mean(.prediction < dat_min),
    p_upper = mean(.prediction < dat_max),
    p_residual = runif(1, p_lower, p_upper),
    z_residual = qnorm(p_residual),
    .groups = "drop_last"
  ) %>%
  ggplot(aes(x = .row, y = z_residual)) +
  geom_point(pch = 1) +
  facet_grid(Sex ~ event) +
  ggtitle("Randomized quantile residuals")

# should be uniform around 0

# qqplot - should be straight line.

fq %>%
  summarise(
    p_lower = mean(.prediction < dat_min),
    p_upper = mean(.prediction < dat_max),
    p_residual = runif(1, p_lower, p_upper),
    z_residual = qnorm(p_residual),
    .groups = "drop_last"
  ) %>%
  ggplot(aes(sample = z_residual)) +
  geom_qq() +
  geom_abline() +
  facet_grid(Sex ~ event) +
  ggtitle("qqplot")
