# check for trends in temperature over time

here::i_am("output/1.site_climate_change.R")

library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(forcats)

factororder <- readRDS(here::here("output/factororder.rds"))

climate_normals_data <- picolaDataClimate::picola_climatebc_locs_normal_monthly

climate_normals <- climate_normals_data %>%
  filter(id == "site") %>%
  select(-starts_with("PPT"), -starts_with("Rad"), -starts_with("Tmin"), -starts_with("Tmax")) %>%
  mutate(enddate = str_extract(normal_period, "\\d{4}$")) %>%
  pivot_longer(cols = starts_with("Tave"), names_to = "Month", values_to = "Tave") %>%
  arrange(enddate) %>%
  mutate(season = ifelse(Month %in% c("Tave03", "Tave04", "Tave05", "Tave06"), "March-June", "other"),
         enddate = as.numeric(enddate),
         normal_period = fct_reorder(normal_period, enddate),
         Site = factor(Site, levels = factororder$site))

saveRDS(climate_normals, here::here("output/monthly_climate_normals.rds"))
