# explore relationships between factor effect estimates and things that could potentially be influencing them that would indicate issues with the model

here::i_am("output/5.effectinvestigation_site.R")

library(dplyr)
library(broom)
library(ggdist)
library(ggplot2)
library(cols4all)

sitedat <- readRDS(here::here('output/tables/sitedat.rds')) %>%
  rename(siteMAT = MAT)
siter <- readRDS(here::here("tmp/siter.rds"))
sitersum <- siter %>%
  group_by(Sex, event, level) %>%
  summarise(median_hdci(.value), sd = sd(.value)) %>%
  rename(Site = level) %>%
  left_join(sitedat)
phenf <- readRDS(here::here("output/phenf.rds"))


# check site MAT vs site effect

ggplot(sitersum, aes(x = siteMAT, y = y)) +
  geom_point() +
  facet_grid(Sex ~ event)

## nope! phew!

# how many years of data at each site
countyrs <- phenf %>%
  select(Site, Year) %>%
  distinct() %>%
  group_by(Site) %>%
  summarise(nYears = length(Year))

# summarise site effects and join with years of data
siteryrs <-  sitersum %>%
  left_join(countyrs)

# do years of data influence site effect estimate?
ggplot(siteryrs, aes(x = nYears, y = y)) +
  geom_point() +
  facet_grid(Sex ~ event)

ggplot(siteryrs, aes(x = sd, y = y)) +
  geom_point(alpha = 0.2) +
  facet_grid(Sex ~ event)

# no

#compare provenance range to ?? site effect?

length(unique(phenf$Genotype))
length(unique(phenf$MAT))

metrange <- phenf %>%
  select(Site, MAT) %>%
  distinct() %>%
  group_by(Site) %>%
  summarise(matrange = max(MAT) - min(MAT)) %>%
  arrange(matrange)

sitersum2 <- sitersum %>%
  left_join(metrange)

# does site effect depend on range of mat provs grown at that site?
ggplot(sitersum2, aes(x = matrange, y = y)) +
  geom_point() +
  facet_grid(Sex ~ event)
## no

# does sd of site effect depend on range of mat provs grown at that site?
ggplot(sitersum2, aes(x = matrange, y = sd, colour= Site)) +
  geom_point(size = 3) +
  facet_grid(Sex ~ event) +
  scale_colour_discrete_c4a_div("dark2")
## yes, the higher the mat range covered at a site, the lower the sd of the site effect estimate

matrangeVSsd_model_results <- sitersum2 %>%
  group_by(Sex, event) %>%
  do({
    model <- lm(sd ~ matrange, data = .)
    tidy_model <- tidy(model)
    rsq <- summary(model)$r.squared
    tidy_model$rsq <- rsq
    tidy_model
  })

sitesdvsprovmat_rsq_range <- round(range(matrangeVSsd_model_results$rsq), digits = 2) *100

saveRDS(sitesdvsprovmat_rsq_range, here::here('output/sitesdvsprovmat_rsq_range.rds'))

#does site effect sd depend on nyears of observation

ggplot(siteryrs, aes(x = nYears, y = sd, colour= Site)) +
  geom_point(size = 3) +
  facet_grid(Sex ~ event) +
  scale_colour_discrete_c4a_div("dark2")

nYearsVSsd_model_results <- siteryrs %>%
  group_by(Sex, event) %>%
  do({
    model <- lm(sd ~ nYears, data = .)
    tidy_model <- tidy(model)
    rsq <- summary(model)$r.squared
    tidy_model$rsq <- rsq
    tidy_model
  })

nYearsVSsd_model_results
## eh maybe? kind of?
