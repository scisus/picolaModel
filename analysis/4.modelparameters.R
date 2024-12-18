# this script reports the phenology model results

here::i_am("analysis/4.modelparameters.R")

# depends ####
library(dplyr)
library(purrr)
library(tidybayes)
library(tidyr)


source(here::here('model/phenology_functions.R'))

factororder <- readRDS(here::here("output/factororder.rds"))

# models ####
modells <- readRDS(here::here("model/modells.rds"))

# globals ####
nsamp <- 6000 # how many samples from the posterior (full posterior is big and slow)
seed <- 738

# data ####

# phenology data
phenf <- readRDS("output/phenf.rds")

# extract and summarise parameter values from the posterior
tidybayes::get_variables(modells$fb)


labdf <- data.frame(Sex = c("FEMALE", "FEMALE", "MALE", "MALE"),
                    event = c('begin', 'end', 'begin', 'end'),
                    model = c('fb', 'fe', 'mb', 'me'))
saveRDS(labdf, "tmp/labdf.rds")

# slope ####
slopes <- purrr::map(modells, gather_slope_draws) %>%
  bind_rows(.id = "model") %>%
  left_join(labdf)
saveRDS(slopes, file = "tmp/slopes.rds")

# intercept ####
intercepts <- purrr::map(modells, gather_intercept_draws) %>%
  bind_rows(.id = "model") %>%
  left_join(labdf) # label the models for plotting
saveRDS(intercepts, file = here::here("tmp/intercepts.rds")) # this is equivalent to the mean in model with no mat effect

interceptssummary <- intercepts %>%
  group_by(Sex, event) %>%
  mean_hdci(.value)

# variation ####

variation <- purrr::map(modells, gather_var_draws) %>%
  bind_rows(.id = "model") %>%
  left_join(labdf) %>%
  mutate(.variable = case_when(.variable != "sigma" ~ stringr::str_sub(.variable, 4, -12),
                               .variable == "sigma" ~ "sigma")) %>%
  mutate(.variable = factor(.variable)) %>%
  mutate(.variable = forcats::fct_relevel(.variable, "sigma", "Year", "Site",  "Genotype", "Tree"))
saveRDS(variation, file = here::here("tmp/variation.rds"))

varsummary <- variation %>%
  group_by(.variable, event, Sex) %>%
  mean_hdci(.value)

saveRDS(varsummary, here::here('output/varsummary.rds'))

# offsets ####

offsets_raw <- purrr::map(modells, gather_offset_draws) %>%
  bind_rows(.id = "model") %>%
  left_join(labdf) #%>%

# turn brms .variable names into useful names (slow)
varlevel <- offsets_raw$.variable %>% stringr::str_split_fixed("[_\\[\\,]", n=4) %>% data.frame() %>%
  select("X2", "X3")
colnames(varlevel) <- c("factor", "level")

# order factors and factor levels
yeargenotypetree <- filter(varlevel, factor %in% c("Year", "Genotype", "Tree")) %>% distinct()
yctorder <- sort(yeargenotypetree$level)

syctorder <- unique(c(factororder$site, yctorder))

offsets <- offsets_raw %>% cbind(varlevel) %>%
  ungroup() %>%
  mutate(factor = forcats::fct_relevel(factor, "Year", "Site", "Genotype", "Tree")) %>%
  mutate(level = forcats::fct_relevel(level, syctorder))

# (slow)
offsets_summary <- offsets %>%
  group_by(model, Sex, event, factor, level) %>%
  median_hdci(.value) %>%
  ungroup()
saveRDS(offsets_summary, file = here::here("output/offsets_summary.rds"))

siter <- filter(offsets, factor == "Site") %>%
  mutate(level = forcats::fct_relevel(level, factororder$site))
saveRDS(siter, file = here::here("tmp/siter.rds"))

genotyper <- filter(offsets, factor == "Genotype")
saveRDS(genotyper, file = here::here("output/genotyper.rds"))

yearr <- filter(offsets, factor == "Year") %>%
  mutate(level = forcats::fct_relevel(level, as.character(factororder$year)))
saveRDS(yearr, file = here::here("output/yearr.rds"))


# priors ####
# at some point, want to go back and compare priors
get_variables(modells$fb)[grep("prior", get_variables(modells$fb),fixed = TRUE)]
