# Creates objects used to build figure describing replication structure
# Identify clones represented across multiple years
# Clones represented in multiple orchards
# Clones represented by multiple trees in the same orchard

here::i_am("output/3.replication_points.R")

library(dplyr)

modeldata <- readRDS(here::here("tmp/datlist.rds")) %>%
  bind_rows()

simpledat <- modeldata %>%
  select(Tree, Genotype, Year, Site) %>% distinct()

# Clones represented across multiple years

mys <- modeldata %>%
  select(Genotype, Year) %>%
  distinct() %>%
  group_by(Genotype) %>%
  summarise(nyears = n()) %>%
  mutate(yearstf = case_when(nyears == 1 ~ FALSE,
                             nyears > 1 ~ TRUE))

# Clones represented in multiple orchards (sites)

mss <- modeldata %>%
  select(Genotype, Site) %>%
  distinct() %>%
  group_by(Genotype) %>%
  summarise(nsites = n()) %>%
  mutate(sitestf = case_when(nsites == 1 ~ FALSE,
                             nsites > 1 ~ TRUE))

# Clones represented by multiple trees in the same orchard

mtss <- modeldata %>%
  select(Tree, Genotype, Site) %>%
  distinct() %>%
  group_by(Genotype, Site) %>%
  summarise(ntrees = n()) %>%
  mutate(treestf = case_when(ntrees == 1 ~ FALSE,
                             ntrees > 1 ~ TRUE))


# join and identify no replication
# All genotypes that are replicated across sites are replicated across years.
replication_points <- mys %>%
  full_join(mss) %>%
  full_join(mtss) %>%
  group_by(Genotype, Site) %>%
  mutate(replicated = base::any(yearstf, sitestf, treestf),
         `Genotype replication` = case_when(replicated == FALSE ~ 'Not replicated',
                                 yearstf == TRUE & treestf == FALSE & sitestf == FALSE ~ 'Across years only',
                                 treestf == TRUE & yearstf == FALSE & sitestf == FALSE ~ 'Within sites only',
                                 treestf == FALSE & yearstf == TRUE & sitestf == TRUE ~ 'Across years and sites',
                                 treestf == TRUE & yearstf == TRUE & sitestf == FALSE ~ 'Across years and within sites',
                                 treestf == TRUE & yearstf == TRUE & sitestf == TRUE ~ 'Across years and sites and within sites')) %>%
  ungroup()


replication_points$`Genotype replication` <- factor(replication_points$`Genotype replication`,
                                       levels = c("Not replicated",
                                                  "Within sites only",
                                                  "Across years only",
                                                  "Across years and within sites",
                                                  "Across years and sites",
                                                  "Across years and sites and within sites"),
                                       ordered = TRUE)


# save all replication info
saveRDS(replication_points, here::here("output/replication_points.rds"))
