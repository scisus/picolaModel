# Order factors for plotting
#

here::i_am("analysis/3.factororder.R")

# depends
library(dplyr)

siteclim <- picolaDataClimate::picola_PNWNAmet_adjusted
provclim <- picolaDataClimate::picola_SPU_climsum
spudat <- picolaDataFlowering::picola_SPUs

phenf <- readRDS(here::here("output/phenf.rds"))

# order site, prov, and year levels by MAT
siteMAT <- siteclim %>%
  mutate(Year = lubridate::year(Date)) %>%
  group_by(Site) %>%
  summarise(MAT = mean(mean_temp_corrected)) %>%
  arrange(MAT)

sitefactororder <- siteMAT$Site

yearMAT <- siteclim %>%
  mutate(Year = lubridate::year(Date)) %>%
  right_join(data.frame(Year = as.numeric(unique(phenf$Year)))) %>%
  group_by(Year) %>%
  summarise(MAT = mean(mean_temp_corrected)) %>%
  arrange(MAT)

yearfactororder <- yearMAT$Year

factororder <- list(site = sitefactororder, year = yearfactororder)
saveRDS(factororder, here::here("output/factororder.rds"))

