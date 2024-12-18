here::i_am("analysis/5.model_check_independentdata.R")

library(dplyr)
library(purrr)
library(tidybayes)
library(tidyr)
library(lubridate)

# provenance and gdd data from Nilsson 1981
nilssongdd <- read.csv(here::here('inputs/Nilsson1981/swedish-pollen-timeseries.csv')) %>%
  rename(provid = provenance)
nilssonprovs <- read.csv(here::here('inputs/Nilsson1981/swedish-sites-with-MATs.csv')) %>%
  rename(Provenance = Name, provid = This.study, MAT = MAT.1961.1990) %>%
  filter(Used.Ange == "X") %>%
  select(Provenance, provid, MAT) %>%
  mutate(provid = as.character(provid), Site = "Central Sweden")
nilssondat <- left_join(nilssongdd, nilssonprovs)

# provenance and date data from o'reilly and owens 1988

ooprovsdat <- read.csv(here::here('inputs/oreilly-and-owens-1988/oreilly-and-owens-1988-table-1.csv'))
oophendat <- read.csv(here::here('inputs/oreilly-and-owens-1988/oreilly-and-owens-1988-table-2.csv'))

# format o'reilly and owens data and add forcing information
dailyforc_oo <- read.csv(here::here("inputs/forcing/dailyforc_1945_2012.csv")) %>%
  group_by(Site, Year) %>%
  mutate(index = cur_group_id()) %>% ungroup() %>%
  filter(Year == "1983", Site == 'PGTIS')

ooprovs <- ooprovsdat%>%
  rename(MAT = MAT.1961.1990, Provenance = Provenance.name.) %>%
  mutate(provid = gsub(".*\\((.*)\\).*", "\\1", Location.description)) %>%
  select(Provenance, provid, MAT)

# tree ages - 14, just reaching maturity in provenance trials at prince george

oophen <- oophendat %>%
  rename(provid = Provenance) %>%
  select(provid, starts_with("Seed"), starts_with("Pollen"), -ends_with("3.Date")) %>%
  pivot_longer(cols = contains("Date"), names_to = "event", values_to = "Date") %>%
  mutate(Sex = case_when(
    grepl("^Seed", event) ~ "FEMALE",
    grepl("^Pollen", event) ~ "MALE"),
    event = case_when(
      grepl("[26]", event) ~ "begin",
      grepl("[74]", event) ~ "end"),
    # Date to DoY
    Date_with_year = paste(Date, "1983"),
    Date_parsed = mdy(Date_with_year),
    DoY = yday(Date_parsed) ) %>%
  select(-starts_with("Date")) %>%
  left_join(ooprovs) %>%
  left_join(dailyforc_oo) %>% # add forcing
  select(event, Provenance, provid, Year, Sex, sum_forcing, MAT) %>%
  mutate(Site = "Central BC")

#model data
modells <- readRDS(here::here("model/modells.rds"))

# combine nilsson & oreilly & owens data and prep for model comparison
independentdat <- full_join(nilssondat, oophen ) %>%
  filter(!is.na(MAT), !is.na(sum_forcing)) %>%
  split(list(.$event, .$Sex))


# posterior prediction #
# ignoring random effects) #######

indpred <- purrr::map2(independentdat, modells,
                              .f = function(x,y) {add_predicted_draws(newdata = x,
                                                                      object = y,
                                                                      re_formula = NA, #NA ignore random effects, NULL include random effects
                                                                      ndraws = 6000)}) %>%
  bind_rows()

indpredsummary <- indpred %>%
  group_by(MAT, Year, Site, provid, Provenance, event, Sex, sum_forcing) %>%
  median_hdci(.prediction, .width = c(0.5, 0.95)) %>%
  mutate(ssp. = case_when(provid == 4846 ~ "contorta",
                          provid != 4846 | is.na(provid) ~ "latifolia"))

saveRDS(indpredsummary, "tmp/indpredsummary.rds")


#MAT range -4.5 to 5.4

