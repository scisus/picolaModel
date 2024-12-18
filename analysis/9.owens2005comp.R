# Compare pollen catch and cone receptivity dates in owens 2005 figures to model output ####
# Owens, J.N., Bennett, J., L’Hirondelle, S., 2005. Pollination and cone morphology affect cone and seed production in lodgepole pine seed orchards. Canadian Journal of Forest Research 35, 383–400. https://doi.org/10.1139/x04-176
#
# I'm looking to
#
# 1) make sure that the date ranges observed in Owens2005 are not totally incongruent with the dates observed in my dataset at those sites in 2000 (for sites where I also have observations)
# 2) check that model predictions overlap Owens' observations

library(dplyr)

here::i_am("analysis/9.owens2005comp.R")

# Calculate the percentage overlap between two date ranges
#
# example of use:
#
# range1 <- c(min = as.Date("2024-06-01"), max = as.Date("2024-06-12"))
# range2 <- c(min = as.Date("2024-06-05"), max = as.Date("2024-06-14"))
# overlap_pct(range1, range2)
overlap_pct <- function(range1, range2) {
  stopifnot(length(range1) == 2, length(range2) == 2)
  stopifnot(range1[1] <= range1[2], range2[1] <= range2[2])
  overlap <- max(0, min(range1[2], range2[2]) - max(range1[1], range2[1]))
  # denominator is the size of the first range
  denom <- range1[2] - range1[1]
  pct <- 100.0 * (as.double(overlap) / as.double(denom))
  message("Date-range overlap: ", overlap, ", first range: ", denom, ", overlap percentage: ", pct)
  return(pct)
}

# Check a data range in owens data overlaps with a model prediction date range by at least a certain percentage.
check_range_overlap <- function(owens_data, pp_data, site, year, sex, pct_thresh) {
  owens_begin <- owens_data %>% filter(Site == site, Year == year, Sex == sex, event == "begin") %>% pull(Date)
  owens_end <- owens_data %>% filter(Site == site, Year == year, Sex == sex, event == "end") %>% pull(Date)
  pp_begin <- pp_data %>% filter(Site == site, Year == year, Sex == sex, event == "begin") %>% pull(DoY) %>% min()
  pp_end <- pp_data %>% filter(Site == site, Year == year, Sex == sex, event == "end") %>% pull(DoY) %>% max()
  day1 <- as.Date(paste0(year, "-01-01"))
  owens_range <- c(min = as.Date(owens_begin), max = as.Date(owens_end))
  pp_range <- c(min = as.Date(pp_begin - 1, origin = day1), max = as.Date(pp_end - 1, origin = day1))
  message("Checking for ", pct_thresh, "% overlap of ranges for ", site, " ", year, " ", sex)
  message("Owens range: [", owens_range[1], ", ", owens_range[2], "], prediction range: [", pp_range[1], ", ", pp_range[2], "]")
  overlap <- overlap_pct(owens_range, pp_range)
  stopifnot(overlap >= pct_thresh)
}

phenf <- readRDS(here::here('output/phenf.rds')) # my phenology data
dailyforc <- read.csv(here::here("inputs/forcing/dailyforc_1945_2012.csv")) %>% # forcing data for obs year in owens paper (2000)
  group_by(Site, Year) %>%
  mutate(index = cur_group_id()) %>% ungroup() %>%
  filter(Year == "2000")

doy_annual_avg_pp_sum <- readRDS(here::here('tmp/doy_annual_avg_pp_sum.rds')) %>% # model predictions for year 2000
  filter(Year == 2000) %>%
  mutate(Date = as.Date(DoY - 1, origin = "1999-12-31"))

owens_data <- read.csv(here::here("inputs/owens2005/owens2005fig1thru3.csv"))

# PRT 2000 ####
# compare fig1 in owens to my observations
fig1 <- phenf %>%
  filter(Site == "PRT", Year == "2000")

# I don't have observations from PRT in 2000 in my data, but Owens observes pollen shed from 5/10 to 5/25 and receptivity from 5/15 to 5/28

dailyforc %>%
  filter(Year == 2000, Site == "PRT", Date %in% c("2000-05-10", "2000-05-28"))
# sum forcing for PRT in 2000  for that date range is 185-314 GDD

doy_annual_avg_pp_sum %>%
  filter(Site == "PRT") %>%
  arrange(Sex, Date)
# May 5 to 21 are the model predictions for receptivity, which overlaps Owens' observations by 1 week
# May 5 to 22 are the model predictions for pollen shed, which overlaps Owens' observations by 13 days
check_range_overlap(owens_data, doy_annual_avg_pp_sum, "PRT", 2000, "FEMALE", 50)
check_range_overlap(owens_data, doy_annual_avg_pp_sum, "PRT", 2000, "MALE", 50)

# PGTIS 2000 ####
# compare fig2 in owens to my observations
fig2 <- phenf %>%
  filter(Site == "PGTIS", Year == "2000", Event_Obs %in% c(2,3)) %>%
  group_by(Sex) %>%
  summarise(earliest = min(Date), latest = max(Date), minforc = min(sum_forcing), maxforc = max(sum_forcing))

# graph shows pollen on pollen monitor from 5/30 to 6/10 and receptive cones from 6/3 to 6/13 and my observations show both from June 1 to June 16 for pollen shed and through the 18th for receptivity.
#
dailyforc %>%
  filter(Year == 2000, Site == "PGTIS", Date %in% c("2000-05-30","2000-06-13"))
#sum forcing 152-253

doy_annual_avg_pp_sum %>%
  filter(Site == "PGTIS") %>%
  arrange(Date)
# May 30 to June 19
check_range_overlap(owens_data, doy_annual_avg_pp_sum, "PGTIS", 2000, "FEMALE", 50)
check_range_overlap(owens_data, doy_annual_avg_pp_sum, "PGTIS", 2000, "MALE", 50)

# Kalamalka 2000 ####
fig3 <- phenf %>%
  filter(Site == "Kalamalka", Year == 2000)

# i don't have this in my data but pollen shed is from 5/7 to 5/22 and receptivity is from 5/8 to 5/21

dailyforc %>%
  filter(Year == 2000, Site == "Kalamalka", Date %in% c("2000-05-07", "2000-05-21"))
# sum forcing 185-284

doy_annual_avg_pp_sum %>%
  filter(Site == "Kalamalka") %>%
  arrange(Date)
# May 3 to 18
check_range_overlap(owens_data, doy_annual_avg_pp_sum, "Kalamalka", 2000, "FEMALE", 50)
check_range_overlap(owens_data, doy_annual_avg_pp_sum, "Kalamalka", 2000, "MALE", 50)
