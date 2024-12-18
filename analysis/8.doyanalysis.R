here::i_am("analysis/8.doyanalysis.R")

library(purrr)
library(corrr)
library(geosphere)
library(broom)
library(tidyr)
library(dplyr)
library(ggplot2)

# data ###
doy_annual_exp_sum <- readRDS(here::here('tmp/doy_annual_exp_sum.rds'))

# rank patterns #####
# Are early years early across all sites? Are late years late across all sites?
#
# Patterns of flowering across seed orchards are correlated; early flowering years at one site tend to be early flowering years at other sites and the same for late flowering years.
#Flowering patterns are more highly correlated for sites that are closer together.


# rank years at each site from earliest to latest flowering
rank_years <- doy_annual_exp_sum %>%
  group_by(Sex, event, Site) %>%
  mutate(DoY_Rank = rank(DoY, ties.method = "min")) %>%
  ungroup() %>%
  select(Sex, event, Site, Year, DoY_Rank) %>%
  pivot_wider(names_from = Site, values_from = DoY_Rank) %>%
  drop_na()

# how many ties are there? (should I use spearman's or kendall's?)

# Define the site columns
site_columns <- c("Border", "Trench", "PGTIS", "KettleRiver", "Sorrento", "Tolko", "PRT", "Vernon", "Kalamalka")

# Generate unique pairs of site columns using combn
site_pairs <- as.data.frame(t(combn(site_columns, 2)), stringsAsFactors = FALSE)
colnames(site_pairs) <- c("site1", "site2")

# Function to count ties for a pair of sites
count_site_ties <- function(data, site1, site2) {
  data %>%
    group_by(Sex, event) %>%
    summarise(ties = sum(!!sym(site1) == !!sym(site2)), .groups = 'drop') %>%
    mutate(Site1 = site1, Site2 = site2)
}

# Apply the function to all unique site pairs and bind the results together
ties <- bind_rows(lapply(1:nrow(site_pairs), function(i) {
  count_site_ties(rank_years, site_pairs$site1[i], site_pairs$site2[i])
})) %>%
  arrange(desc(ties))

max_ties = round(ties$ties[1]/length(unique(rank_years$Year)), digits = 2) * 100
saveRDS(max_ties, here::here("output/max_ties.rds"))

#Calculate Kendall’s Tau correlation for each combination of Sex and event

siteorder <- readRDS(here::here('output/factororder.rds'))$site
#Calculate Kendall’s Tau correlation for each combination of Sex and event
rank_correlation <- rank_years %>%
  group_by(Sex, event) %>%
  nest() %>%
  mutate(correlation_matrix = purrr::map(data, ~ corrr::correlate(.x %>% select(-Year), method = "kendall", use = "pairwise.complete.obs"))) %>%
  select(-data) %>%
  unnest(correlation_matrix) %>%
  pivot_longer(cols = Border:Kalamalka,
               names_to = "Site2",
               values_to = "Correlation") %>%
  rename(Site1 = "term") %>%
  ungroup() %>%
  mutate(Site1 = forcats::fct_relevel(Site1, siteorder), Site2 = forcats::fct_relevel(Site2, siteorder)) %>%
  mutate(Correlation = case_when(Site1 == Site2 ~ 1,
                                 Site1 != Site2 ~ Correlation)) # assign correlation of 1 for self



# geographical isolation measure (IBD)
sitedat <- picolaDataClimate::picola_climatebc_locs_normal_1961_1990 %>%
  filter(id == "site")
sitedat <- sitedat %>%
  select(Site, Latitude, Longitude, Elevation, MAT)

# calculate distances between sites in km
distances <- sitedat %>%
  select(Longitude, Latitude) %>%
  as.matrix() %>%
  geosphere::distm(., fun = distHaversine)

distances_df <- as.data.frame(as.table(distances)) %>%
  mutate(Distance = Freq/1000, .keep = "unused") # convert from m to km
names(distances_df) <- c("Site1", "Site2", "Distance")
distances_df$Site1 <- sitedat$Site[distances_df$Site1]
distances_df$Site2 <- sitedat$Site[distances_df$Site2]

# add geographic isolation measures to rank corr df and drop duplicates
rank_correlation_wdist <- rank_correlation %>%
  left_join(distances_df, relationship = "many-to-many") %>%
  group_by(Sex, event) %>%
  distinct(Distance, .keep_all = TRUE)
saveRDS(rank_correlation_wdist, here::here('tmp/rank_correlation_wdist.rds'))

# Building models and extracting tidy summaries

corr_model_results <- rank_correlation_wdist %>%
group_by(Sex, event) %>%
  do({
    fitted_model = lm(Correlation ~ Distance, data = .)
    tidy_summary = tidy(fitted_model)
    glance_summary = glance(fitted_model)
    bind_cols(glance_summary, tidy_summary)
  }) %>%
  ungroup() %>%
  select(Sex, event, r.squared, p.value = p.value...7, term, estimate, std.error)
saveRDS(corr_model_results, here::here('tmp/corr_model_results.rds'))


## lower for warmer provs at all sites for male begin, but no consistent patterns, not even in gaps (basically 0 to huge)

# warm vs cold provs ####
# avg difference in doy for flowering between coldest and warmest prov at each site

doy_annual_avg_pp_sum <- readRDS(here::here('tmp/doy_annual_avg_pp_sum.rds')) %>%
  filter(.width == 0.95)
doy_annual_avg_pp_sum$MAT_label <- paste("Provenance MAT:", doy_annual_avg_pp_sum$MAT, "\u00B0C")

warmvscold <- doy_annual_avg_pp_sum %>%
  select(-starts_with("."), -MAT_label) %>%
  pivot_wider(names_from = MAT, names_prefix = "MAT", values_from = DoY) %>%
  mutate(diff = MAT6.8 - `MAT-0.7`) %>%
  select(-starts_with("MAT")) %>%
  group_by(Site, event, Sex) %>%
  summarise(mean_doy_diff = mean(diff), sd_diff = sd(diff))
saveRDS(warmvscold, here::here("output/warmvscold.rds"))

# variation in day of year ####
# do different sites vary more than others in flowering?
# do not remove this analysis

doy_annual_avg_pp_sum <- readRDS(here::here('tmp/doy_annual_avg_pp_sum.rds'))
doysd <- doy_annual_avg_pp_sum %>%
  group_by(MAT, Site, event, Sex) %>%
  summarise(sd = sd(DoY))

ggplot(doysd, aes(x = Site, y = sd, colour = MAT)) +
  geom_point() +
  facet_grid(Sex ~ event)

# maybe variation in start for cold provs at warm sites? no/inconclusive.
