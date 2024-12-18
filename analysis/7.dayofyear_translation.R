# convert forcing to day of year

here::i_am("analysis/7.dayofyear_translation.R")

library(purrr)
library(ggplot2)
library(dplyr)
library(lubridate)
library(tidybayes)
library(tidyr)
library(forcats)

source(here::here('model/phenology_functions.R'))
focalsites <- c("Kalamalka", "KettleRiver", "PGTIS", "Trench", "Border") # warmest to coldest
factororder <- readRDS(here::here('output/factororder.rds'))

# daily forcing data ############

# daily "real" forcing
dailyforc <- read.csv(here::here("inputs/forcing/dailyforc_1945_2012.csv")) %>%
  group_by(Site, Year) %>%
  mutate(index = cur_group_id()) %>% ungroup() %>%
  mutate(Site = forcats::fct_relevel(Site, factororder$site))

# from temp mean at each site across 1945-2012
typical_year_forc <- read.csv(here::here("inputs/forcing/typical_year_forc.csv")) %>%
  mutate(Date = as.Date(Date_scale)) %>% select(-Date_scale)
# averaged over 30 year periods
normal_forc <- read.csv(here::here("inputs/forcing/normalforc_1921-2100.csv")) %>%
  group_by(Site, period, scenario) %>%  # index
  mutate(index = cur_group_id()) %>% ungroup() %>%
  mutate(Site = forcats::fct_relevel(Site, factororder$site))

fepred_allsites <- readRDS(here::here("tmp/fepred_allsites.rds"))  %>% ## expectation for trees sourced from all sites - no downsampling
  ungroup() %>%
  mutate(Site = forcats::fct_relevel(Site, factororder$site))

sitedat <- picolaDataClimate::picola_climatebc_locs_normal_1961_1990

# for each site, generate a timeseries of mean temperatures and associated accumulated forcing that reflects the general pattern of temperatures throughout the year. Do this by averaging temperatures on each day between 1945 and 2012 at each site.
# alternate source for this data would be typical_ts.csv in lodgepole_climate project

siteMAT <- sitedat %>%
  filter(id == "site") %>%
  select(Site, MAT, Elevation) %>%
  mutate(MAT = round(MAT, 1))


# HOME vs AWAY ############
# now consider doy expectations in a typical year for trees from each of the sites of interest grown at those sites of interest

## intercept only #######
intercepts <- readRDS(here::here("tmp/intercepts.rds"))
doy_typical_allsites_interceptonly <- map_dfr(split(typical_year_forc, f = list(typical_year_forc$Site), drop = TRUE),
                                              find_day_of_forcing, .id = ".id",
                                              bdf = intercepts, aforce = "sum_forcing", bforce = ".value") %>%
  rename(Site = .id, DoY = newdoycol) %>%
  ungroup() %>%
  select(-.chain, -.iteration, -.draw) %>%
  mutate(Site = forcats::fct_rev(forcats::fct_relevel(Site, factororder$site))) %>%
  left_join(select(typical_year_forc, Date, DoY) %>% distinct()) %>%
  mutate(provenance_effect = FALSE) %>%
  left_join(siteMAT )

# including MAT #######
doy_typical_allsites <- map_dfr(split(typical_year_forc, f = list(typical_year_forc$Site), drop = TRUE),
                       find_day_of_forcing, .id = ".id",
                       bdf = rename(fepred_allsites, Source = Site), aforce = "sum_forcing", bforce = ".epred") %>%
  rename(Site = .id, DoY = newdoycol) %>%
  ungroup() %>%
  select(-.row, -.draw) %>%
  mutate(Site = forcats::fct_rev(forcats::fct_relevel(Site, factororder$site))) %>%
  left_join(select(typical_year_forc, Date, DoY) %>% distinct()) %>%
  left_join(select(ungroup(intercepts), model, Sex, event) %>% distinct) %>% # add sex and event
  mutate(provenance_effect = TRUE)


# provenance effect pushes southern provenances to flower later and northern to flower earlier, shortening the overall flowering period from south to north - and increasing overlap between north and south

doy_typical_allsites_interceptonly_intermediate <- doy_typical_allsites_interceptonly %>%
  filter(Site == "PGTIS") %>%
  select(Sex, event, DoY) %>%
  rename(intercept = DoY) %>%
  group_by(Sex, event) %>%
  median_hdci(intercept) %>%
  filter(.width == 0.95) %>%
  select(-starts_with("."))

# now do all sources grown at a single site
doy_typical_all_at_PGTIS <- doy_typical_allsites %>%
  filter(Site == "PGTIS") %>%
  group_by(Site, MAT, Sex, event) %>%
  median_hdci(DoY) %>%
  filter(.width == 0.95) %>%
  select(-starts_with(".")) %>%
  merge(doy_typical_allsites_interceptonly_intermediate)
saveRDS(doy_typical_all_at_PGTIS, here::here("tmp/doy_typical_all_at_PGTIS.rds"))

 # When all sources are grown at the same Site (PGTIS), MAT effect reduces overlap


# year to year variation ####

## no site posterior prediction, 6000 draws, avg year, genotype, tree, site effects. 1945-2012. ####


# split climate data into list based on site and year
dailyforc_list <- dailyforc %>%
  arrange(Site, Year, DoY) %>%
  split(f = list(.$Site, .$Year), drop = TRUE)
fpred_orch_avg <- readRDS(here::here('tmp/fpred_orch_avg.rds')) %>%
  filter(MAT %in% c(-0.7, 6.8)) %>%
  ungroup() %>%
  select( -Year, -Site, -Tree, -Genotype, -.chain, -.iteration)

# match forcing predictions in fpred_orch_avg to doy in dailyforc_list.
doy_annual_avg_pp <- map_dfr(dailyforc_list, .f = find_day_of_forcing,
                          .id = "index",
                          bdf = fpred_orch_avg,
                          aforce = "sum_forcing",
                          bforce = ".prediction") %>%
  rename(DoY = newdoycol)

# Use strsplit to split the .id column by the period (separate toooo slow)
split_id <- strsplit(doy_annual_avg_pp$index, "\\.")

# Create new Site and Year columns by extracting the split components
doy_annual_avg_pp$Site <- sapply(split_id, `[`, 1)
doy_annual_avg_pp$Year <- as.numeric(sapply(split_id, `[`, 2))

#now summarise
doy_annual_avg_pp_sum <- doy_annual_avg_pp %>%
  select(-index) %>%
  group_by(MAT, Site, Year, event, Sex) %>%
  median_hdci(DoY, .width = c(0.50, 0.95)) %>%
  ungroup() %>%
  mutate(Year = as.numeric(Year)) %>%
  mutate(Site = forcats::fct_relevel(Site, factororder$site)) # order sites
doy_annual_avg_pp_sum$MAT_label <- paste("MAT:", doy_annual_avg_pp_sum$MAT)
saveRDS(doy_annual_avg_pp_sum, here::here("tmp/doy_annual_avg_pp_sum.rds"))

# calculate year to year variancevariance
vardoy <- doy_annual_avg_pp %>%
  rename(provMAT = MAT) %>%
  group_by(Site, event, Sex, provMAT) %>%
  summarise(standarddev = sd(DoY)) %>%
  left_join(siteMAT)

ggplot(vardoy, aes(x = MAT, y = standarddev, colour = as.factor(provMAT))) +
  geom_point(size = 2, alpha = 0.7) +
  facet_grid(event ~ Sex) +
  xlab("Site MAT") +
  ylab("Standard deviation (GDD)") +
  theme_bw()

## expectation and no random effects for y2y var and ranking correlation####
fepred_allsites_ls <- split(fepred_allsites, f = list(fepred_allsites$Site), drop = TRUE)
dailyforc_ls <- split(dailyforc, f = list(dailyforc$Site), drop = TRUE)

stopifnot("fepred_allsites_ls and dailyforc_ls names must be identical (sorted by Site)" = names(dailyforc_ls) == names(fepred_allsites_ls))

doy_annual_exp <- map2_dfr(.x = dailyforc_ls, .y = fepred_allsites_ls, .f = find_day_of_forcing_mapper, bforce = ".epred") %>%
  rename(DoY = newdoycol, index = .id) %>%
  mutate(index = as.numeric(index)) %>%
  ungroup() %>%
  left_join(select(dailyforc, index, Site, Year) %>% distinct())

### and summarize
doy_annual_exp_sum <- doy_annual_exp %>%
  group_by(MAT, Site, event, Sex, Year) %>%
  median_hdci(DoY) %>%
  ungroup() %>%
  mutate(Site = forcats::fct_relevel(Site, factororder$site)) # correct to full sites
doy_annual_exp_sum$MAT_label <- paste("MAT:", doy_annual_exp_sum$MAT)
saveRDS(doy_annual_exp_sum, here::here("tmp/doy_annual_exp_sum.rds"))

# normal periods ########
normal_forc_focal <- normal_forc %>% filter(Site %in% focalsites)
fepred_allsites_focal <- fepred_allsites %>% filter(Site %in% focalsites) %>% ungroup()

climatelist_nff <- split(normal_forc_focal, f = list(normal_forc_focal$Site), drop = TRUE)
phenlist_af <- split(fepred_allsites_focal, f = list(fepred_allsites_focal$Site), drop = TRUE)

stopifnot("climate and phenology data lists must have sites in same order" = names(climatelist_nff) == names(phenlist_af))

doy_normal <- map2_dfr(.x = climatelist_nff, .y = phenlist_af, .f = find_day_of_forcing_mapper) %>%
  rename(index = .id, DoY = newdoycol) %>%
  mutate(index = as.numeric(index)) %>%
  ungroup() %>%
  select(-.row) %>%
  left_join(select(normal_forc_focal, index, Site, period, scenario) %>% distinct()) %>%
  mutate(Site = forcats::fct_rev(forcats::fct_relevel(Site, factororder$site)))

 # graph climate change normals ####

doy_normal_subset <- doy_normal %>%
  filter(period %in% c("1951-1980", "1981-2010", "2011-2040", "2041-2070", "2071-2100"),
         Site %in% c("Kalamalka", "KettleRiver", "PGTIS", "Trench", "Border"),
         scenario %in% c("historical", "ssp245", "ssp585")) %>%
  mutate(Date = lubridate::ymd("2023-12-31") + DoY)
saveRDS(doy_normal_subset, "output/doy_normal_subset.rds")
