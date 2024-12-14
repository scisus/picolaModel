# Add forcing to climate data

here::i_am("data-raw/calc_forcing.R")

library(dplyr)
library(lubridate)

# functions #####

# calculate accumulated forcing
calc_forcing <- function(vec, threshold = 5) {
  forc <- which(vec > threshold)

  # calculate forcing
  fcol <- rep(0, length(vec))
  fcol[forc] <- vec[forc] - threshold

  return(fcol)
}

calc_accumulated_forcing <- function(df, grouping_vars, mean_temp = mean_temp, threshold = 5) {
  ndf <- df %>%
    mutate(forcing = calc_forcing(vec = {{ mean_temp }}, threshold = threshold)) %>%
    group_by(across(all_of(grouping_vars))) %>%
    mutate(sum_forcing = cumsum(forcing))

  return(ndf)

}

# globals #####
# set threshold for gdd calculation
th <- 5

# data #####

pnwnamet_adj <- picolaDataClimate::picola_PNWNAmet_adjusted %>% # read.csv("../lodgepole_climate/processed/PNWNAmet_adjusted.csv") %>%
  select(Date, Site, mean_temp_corrected) %>%
  rename(mean_temp = mean_temp_corrected) %>%
  mutate(DoY = yday(Date), Year = lubridate::year(Date))

typical_year <- picolaDataClimate::picola_typical_ts %>% # read.csv("../lodgepole_climate/processed/typical_ts.csv") %>%
  rename(mean_temp = mean_mean_temp)
normals <- picolaDataClimate::picola_normal_daily # read.csv("../lodgepole_climate/processed/normal_daily.csv")

# forcing calculation #####
dailyforc_1945_2012 <- calc_accumulated_forcing(pnwnamet_adj,
                                                grouping_vars = c("Site", "Year"),
                                                threshold = th)
typical_year_forc <- calc_accumulated_forcing(typical_year,
                                              grouping_vars = c("Site"),
                                              threshold = th)
normal_forc <- calc_accumulated_forcing(normals,
                                        grouping_vars = c("Site", "period", "scenario"),
                                        threshold = th)


# write out files #####

write.csv(dailyforc_1945_2012, here::here("data-raw/forcing/dailyforc_1945_2012.csv"), row.names = FALSE)
write.csv(typical_year_forc, here::here("data-raw/forcing/typical_year_forc.csv"), row.names = FALSE)
write.csv(normal_forc, here::here("data-raw/forcing/normalforc_1901-2100.csv"), row.names = FALSE)




