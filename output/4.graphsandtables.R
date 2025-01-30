# visualizations

here::i_am("output/4.graphsandtables.R")

library(dplyr)
library(ggplot2)
library(forcats)
library(ggbeeswarm)
library(tidybayes)
library(patchwork)    # combine plots with simple syntax
library(ggokabeito)   # Neat accessible color palette
library(tidyr)
library(ggrepel)
library(cols4all)
library(ggh4x)

#section headings match fig references/chunk names in paper
factororder <- readRDS(here::here("output/factororder.rds")) # site names in order from coldest to warmest MAT

# to rebuild rangemap
# source('6.map.R')

# sites ####
sitedat <- picolaDataClimate::picola_climatebc_locs_normal_1961_1990 %>%
  filter(id == "site") %>%
  select(Site, Latitude, Longitude, Elevation, MAT) %>%
  arrange(MAT)
sitedat$Type <- c("comparison", "comparison", rep("seed orchard", 7))
saveRDS(sitedat, here::here("output/tables/sitedat.rds"))

# typicalsiteclim - forcing and climate ##########
typical_year_forc <- read.csv(here::here("inputs/forcing/typical_year_forc.csv")) %>% # from temp mean at each site across 1945-2012
  mutate(Date = as.Date(Date_scale)) %>%
  select(-Date_scale) %>%
  mutate(Site = factor(Site, levels = factororder$site)) %>%
  mutate(Site_type = case_when(Site %in% c("Border", "Trench") ~ "comparison",
                               !Site %in% c("Border", "Trench") ~ "orchard")) %>%
  mutate(Site_type = fct_rev(Site_type)) %>%
  filter(DoY < 180)

meantempplot <- ggplot(typical_year_forc, aes(x = Date, y = mean_temp, color = Site, linetype = Site_type)) +
  geom_hline(yintercept = 5, colour = "dark grey") +
  geom_line() +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b") +
  ggtitle("Mean daily temperature") +
  scale_colour_viridis_d() +
  theme_bw() +
  ylab("Temperature (\u00B0C)") + xlab("") +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 8))

forcplot <- ggplot(typical_year_forc, aes(x = Date, y = forcing, color = Site, linetype = Site_type)) +
  geom_line() +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b") +
  ggtitle("Daily forcing") +
  scale_color_viridis_d() +
  theme_bw() +
  ylab("Growing Degree Days (\u00B0C)") + xlab("") +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 8))

sumforcplot <- ggplot(typical_year_forc, aes(x = Date, y = sum_forcing, color = Site, linetype = Site_type)) +
  geom_line() +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b") +
  ggtitle("Forcing accumulation") +
  scale_color_viridis_d() +
  theme_bw() +
  ylab("Growing Degree Days (\u00B0C)") + xlab("") +
  labs(linetype = "Site type", color = "Site") +
  guides(
    color = guide_legend(ncol = 1),        # Color guide with 3 columns
    linetype = guide_legend(ncol = 1)      # Linetype guide with 1 column (stacked vertically)
  ) +
  # Move the legend inside the top-left corner of the plot
  theme(
    legend.position = "right",
    legend.background = element_rect(fill = alpha("white", 0.5)),  # Semi-transparent background
    # Decrease the font size of the legend text and title
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  )
sumforcplot

siteclimplot <- meantempplot / forcplot / sumforcplot +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = 'collect')
siteclimplot
ggsave(here::here("output/figures/siteclimplot.png"), width = 5, height = 7.5)

## monthlytemps - monthly climate normals ####

monthly_climate_normals <- readRDS(here::here("output/monthly_climate_normals.rds"))
ggplot(monthly_climate_normals, aes(x = normal_period, y = Tave, group = Month, colour = season)) +
  geom_point(size = 0.5) +
  geom_line(linewidth = 0.25) +
  scale_colour_viridis_d(option = "rocket", end = 0.8, direction = -1) +
  facet_wrap("Site") +
  ggtitle("Monthly mean temperature", subtitle = "overlapping 30 year climate normal periods") +
  xlab("Climate normal period") +
  ylab("Monthly average temperature (\u00B0C)") +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust=1) )

ggsave(here::here("output/figures/monthly_climate_normals.png"), width = 6, height = 7)

# MAT - (and replication) #####
sites <- picolaDataClimate::picola_climatebc_locs_normal_1961_1990 %>%
  filter(id == "site") %>%
  select(Site, MAT) %>%
  mutate(Site = forcats::fct_rev(forcats::fct_reorder(Site, MAT)))
sites$`Site Type` <- c(rep("Seed Orchard", 7), rep("Comparison", 2))

replication_points <- readRDS(here::here("output/replication_points.rds"))

provs <- readRDS(here::here("output/phenf.rds")) %>%
  select(Tree, MAT, Site, Genotype) %>%
  distinct() %>%
  left_join(replication_points) %>%
  rename('Within Sites' = treestf, 'Across Sites' = sitestf, 'Across Years' = yearstf, Replicated = replicated) %>%
  mutate(Site = forcats::fct_relevel(Site, factororder$site))

siteplot <- ggplot(data=sites) +
  geom_point(aes(x = "Sites", y = MAT, shape = `Site Type`)) +
  geom_text_repel(aes(x = "Sites", y = MAT, label = Site), size = 2, point.padding = 0.05, min.segment.length = 0.16) +
  xlab("") +
  ylab("Mean Annual Temperature (\u00B0C)") +
  scale_y_continuous(limits = c(min(sites$MAT), max(sites$MAT))) +
  scale_shape_manual(values = c('Comparison' = 17, 'Seed Orchard' = 16)) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
  ggtitle("Site MATs") +
  guides(shape = guide_legend(nrow = 2, title.position = 'top'))


provplot <- ggplot(provs, aes(x = Site, y = MAT,
                  shape = `Genotype replication`,
                  fill = `Within Sites`,
                  color = `Across Sites`,
                  group = `Genotype replication`)) +
  geom_quasirandom(varwidth = TRUE, alpha = 0.5) +
  scale_shape_manual(values = c("Not replicated" = 4,
                                "Within sites only" = 17,
                                "Across years only" = 21,
                                "Across years and sites" = 21,
                                "Across years and within sites" = 21,
                                "Across years and sites and within sites" = 21)) +
  scale_fill_manual(values = c("TRUE" = "grey28", "FALSE" = "transparent"),
                    guide = "none") +  # Hide separate fill legend
  scale_colour_manual(values = c("TRUE" = "#DF536B", "FALSE" = 'black'),
                      guide = "none") +  # Hide separate color legend
  scale_y_continuous(limits = c(min(provs$MAT), max(provs$MAT)), position = "right") +
  ylab("Mean Annual Temperature (\u00B0C)") +
  theme(legend.position = "bottom") +
  ggtitle("Provenance MATs") +
  guides(shape = guide_legend(override.aes = list(fill = c(rep(c("transparent", "grey28"),3)),
                                                  color = c(rep("black", 4 ), rep("#DF536B", 2))),
                              title.position = 'top',
                              nrow = 3))



siteplot + provplot + patchwork::plot_layout(widths = c(1,4)) + patchwork::plot_annotation(tag_levels = 'A')

ggsave(here::here("output/figures/MAT.png"), width = 7, height = 5, scale = 1.1)

# sampling - sampling events ####
phenf <- readRDS(here::here("output/phenf.rds"))
surveydf <- phenf %>%
  select(Year, Site, Orchard, DoY) %>%
  distinct() %>%
  mutate(Site = forcats::fct_relevel(Site, factororder$site)) %>%
  group_by(Site, Year, Orchard) %>%
  mutate(sampleindex = cur_group_id()) %>%
  ungroup()

yrspersite <- surveydf %>%
  select(Site, Year) %>%
  distinct() %>%
  group_by(Site) %>%
  summarise(nyears = n()) %>%
  mutate(Site = forcats::fct_relevel(Site, factororder$site)) %>%
  arrange(as.factor(Site))

facet_labeller_site <- function(variable, value) {
  c(
    "PGTIS",
    rep("", yrspersite$nyears[1] - 1),
    "KettleRiver",
    rep("", yrspersite$nyears[2] - 1),
    "Sorrento",
    "Tolko",
    rep("", yrspersite$nyears[4] - 1),
    "PRT",
    rep("", yrspersite$nyears[5] - 1),
    "Vernon",
    rep("", yrspersite$nyears[6] - 1),
    "Kalamalka",
    rep("", yrspersite$nyears[7] - 1)
  )
}

# thanks for help from https://stackoverflow.com/questions/54178285/how-to-remove-only-some-facet-labels

ggplot(surveydf, aes(x=DoY, y=as.factor(sampleindex), colour = Site, group = as.factor(Orchard))) +
  geom_point(pch=3) +
  geom_line(alpha = 0.5) +
  facet_grid(rows=vars(Site,Year), scales="free_y",
             labeller = labeller(Site = as_labeller(facet_labeller_site))) +
  scale_color_okabe_ito() +
  scale_shape_manual(values=c(1:7)) +
  theme_bw(base_size = 14) +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        panel.border=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank()) +
  xlab("Day of Year") +
  ggtitle("Observation Dates", subtitle = "for each orchard at each site")
ggsave(here::here("output/figures/sampling.png"), width = 6, height = 8)

# parameters ####
## means ####
intercepts <- readRDS(here::here("tmp/intercepts.rds"))

interceptsummary <- intercepts %>%
  group_by(Sex, event) %>%
  mean_hdci(.value) %>%
  mutate(
    Value_CI = paste0(round(.value), " (", round(.lower), "-", round(.upper), ")")
  ) %>%
  select(-starts_with(".")) %>%
  pivot_wider(names_from = Sex, values_from = Value_CI)

saveRDS(interceptsummary, file = here::here("output/tables/interceptsummary.rds"))

### slopes ####
slopes <- readRDS(here::here("tmp/slopes.rds"))

slopesummary <- slopes %>%
  group_by(Sex, event) %>%
  mean_hdci(.value) %>%
  mutate(
    Value_CI = paste0(round(.value, 2), " (", round(.lower, 2), "-", round(.upper, 2), ")")
  ) %>%
  select(-starts_with(".")) %>%
  pivot_wider(names_from = Sex, values_from = Value_CI)
slopesummary
saveRDS(slopesummary, file = here::here("output/tables/slopesummary.rds"))

## varoffsets ####
### sd
# plot sd parameters using variation from modelparameters.R
variation <- readRDS(here::here("tmp/variation.rds")) %>% ungroup()
varplot <- ggplot(variation, aes(y = forcats::fct_rev(event), x = .value, colour = Sex, shape = event)) +
  stat_pointinterval(position = "dodge") +
  ylab("") +
  xlab("Standard deviation (GDD)") +
  facet_grid(.variable ~ ., scales = "free_y") +
  guides(shape = "none", colour = guide_legend(nrow=2)) +
  scale_colour_discrete_c4a_div(palette = "acadia") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.justification = "right") +
  geom_vline(xintercept = 0, linetype = 3, colour = 'darkgray')
varplot


### offset_medians
# plot medians of offset parameters in point clouds (like beeswarm)
offsets_summary <- readRDS(here::here("output/offsets_summary.rds")) %>%
  mutate(model = factor(model, levels = c("mb", "fb", "me", "fe")))
offsetplot <- offsets_summary %>%
  select(model, Sex, event, factor, level, .value, .point) %>% distinct() %>%
  ggplot(aes(y=forcats::fct_rev(event), x = .value, colour = Sex, shape = event, group = model)) +
  geom_quasirandom(dodge.width = 1) +
  facet_wrap("factor") +
  scale_colour_discrete_c4a_div(palette = "acadia", guide = "none") +
  ylab("") +
  xlab("Offset median (GDD)") +
  theme_bw(base_size = 10) +
  geom_vline(xintercept = 0, linetype = 3, colour = 'darkgray') +
  theme(axis.ticks.y = element_blank(), legend.position = "bottom", legend.justification = "left") +
  guides(shape = guide_legend(nrow=2)) +
  labs(shape = "Event")
offsetplot

varplot + offsetplot +
  plot_layout(widths = c(1,2.5)) +
  plot_annotation(tag_levels = 'A')

ggsave(here::here("output/figures/varoffsets.png"), width = 6.5, height = 5)

# indcomp - independent data comparison ###########

indpredsummary <- readRDS(here::here('tmp/indpredsummary.rds')) %>%
  filter(!(Site == "Central BC" & event == "end"))

## graph

# Position the contorta label
label_x <- 240  # Position the label slightly outside the plot
label_y <- 175  # Vertical position near the contorta points

# Subset the data to only include facets with "contorta" points
contorta_data <- subset(indpredsummary, ssp. == "contorta")

ggplot(indpredsummary, aes(x = sum_forcing, y = .prediction)) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper, colour = as.character(.width)),
                width = 0.1, alpha = 0.5) +
  geom_point(aes(fill = MAT, shape = Site), size = 2) +
  geom_abline(colour = "darkgrey") +
  facet_grid(Sex ~ event) +
  scale_color_manual(values = c("0.5" = "black", "0.95" = "grey"), name = "HDI") +  # For error bars
  scale_fill_viridis_c(option = "D", name = "Provenance MAT") +  # For points
  scale_shape_manual(values = c("Central BC" = 21, "Central Sweden" = 24)) +  # Custom shapes for ssp.
  theme_bw() +
  xlab("Measured forcing (GDD)") +
  ylab("Predicted forcing (GDD)") +
  # Add segments and label only in facets with contorta points
  geom_segment(data = contorta_data,
               aes(xend = label_x, yend = label_y),
               linetype = "dashed", color = "darkgrey") +
  geom_label(data = contorta_data,
             aes(x = label_x, y = label_y, label = "contorta"),
             size = 2, fill = "white", color = "black", hjust = 0) +
  theme(legend.position = "bottom") +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5),  # Correct guide for continuous fill (MAT)
         shape = guide_legend(title.position = "top", title.hjust = 0.5),  # Correct guide for shape
         colour = guide_legend(title.position = "top", title.hjust = 0.5, override.aes = list(size = 2)))  # Correct guide for colour
ggsave(here::here("output/figures/independent_data_comp.png"), width = 7, height = 5)

# retrodictions ####
#
## obsvsretro - Observed vs retrodicted ##############
fretro_summary <- readRDS(here::here("tmp/fretro_summary.rds"))
censorpal <- c4a("icefire", 3)
ggplot(fretro_summary, aes(x = sum_forcing, y = .prediction, color = censored)) +
  geom_point(alpha = .5, shape = 3) +
  facet_grid(Sex ~ event) +
  geom_abline(color = "grey20") +
  xlab("Observed accumulated forcing (GDD)") +
  ylab("Median retrodicted accumulated forcing (GDD)") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_colour_manual(values = c(censorpal[2], censorpal[3], censorpal[1])) +
  theme(legend.position = "bottom") +
  theme_bw()
ggsave(here::here("output/figures/obsvsretro.png"), width = 6, height = 5)

## compareconcept - conceptual diagram to explain retrodictions ##############
overlapex <- data.frame(data = c("observation", "model", "model", "model"),
                        label = c("observation interval", "mean in interval", "one sd overlap", "no overlap"),
                        overlap = c(NA, TRUE, TRUE, FALSE),
                        mean = c(NA, 13, 3, 25),
                        .lower = c(5, 8, -2, 20),
                        .upper = c(15, 18, 8, 30))
overlapex$label <- factor(overlapex$label, levels = unique(overlapex$label))

ggplot(overlapex, aes(x = mean, y = forcats::fct_rev(label), xmin = .lower, xmax = .upper, colour = overlap, )) +
  geom_pointinterval(size = 5, linewidth = 3) +
  ylab("") +
  facet_grid(forcats::fct_rev(data) ~ ., scales = "free_y") +
  scale_color_manual(values = c("TRUE" = "#1B9E77", "FALSE" = "#D95F02")) +
  geom_vline(xintercept = c(5,15), colour = "grey", linetype = 2) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  theme(axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank())

ggsave(here::here("output/figures/obsvsretroconceptual.png"), width = 4, height = 3)


## genpred - GDD predictions for a generic site, year, tree, etc. #################

fpred_orch_avg_summary <- readRDS(here::here("tmp/fpred_orch_avg_summary.rds"))
phenf_orchplot_prov <- readRDS(here::here("output/phenf.rds")) %>%
  filter(Event_Obs %in% c(2,3)) %>%
  mutate(Site = forcats::fct_relevel(Site, factororder$site))

ggplot(fpred_orch_avg_summary) +
  # colored ribbons for start and end
  geom_ribbon(aes(x = MAT, ymin = .lower, ymax = .upper, group = event, fill = event), alpha = 0.3) +
  geom_line(aes(x = MAT, y = .prediction, colour = event)) +
  facet_grid(. ~ Sex) +
  scale_fill_discrete_c4a_div(palette = "icefire") +
  scale_colour_discrete_c4a_div(palette = "icefire") +
  theme_bw() +
  ylab("Accumulated forcing (\u00B0C)") +
  xlab("Provenance Mean Annual Temperature (\u00B0C)") +
  theme(
    axis.text.y = element_text(size = 7),
    strip.text.y = element_text(size = 8),
    legend.position = "bottom"
  )
ggsave(here::here("output/figures/genpred_gdd.png"), width = 6, height = 4)



# DoY

doy_annual_avg_pp_sum <- readRDS(here::here("tmp/doy_annual_avg_pp_sum.rds")) %>%
  filter(.width == 0.95)
doy_annual_avg_pp_sum$MAT_label <- paste("Provenance MAT:", doy_annual_avg_pp_sum$MAT, "\u00B0C")

# phenology data for comparison
phenf_orchplot <- readRDS(here::here("output/phenf.rds")) %>%
  select(-MAT) %>%
  filter(Event_Obs %in% c(2,3)) %>%
  mutate(Year = as.numeric(Year), Site = forcats::fct_relevel(Site, factororder$site))

# orchmale - male only ####
ggplot() +
  geom_line(data = filter(phenf_orchplot, Sex == "MALE"), aes(x = Year, y = DoY, group = Year), alpha = 0.9) +
  geom_ribbon(data = filter(doy_annual_avg_pp_sum, Sex == "MALE"), aes(x = Year, ymin = .lower, ymax = .upper, group = event, fill = event), alpha = 0.3) +
  geom_line(data = filter(doy_annual_avg_pp_sum, Sex == "MALE"), aes(x = Year, y = DoY, colour = event)) +
  scale_fill_discrete_c4a_div(palette = "icefire") +
  scale_colour_discrete_c4a_div(palette = "icefire") +
  theme_bw() +
  xlab("Year") +
  ylab("Date") +
  ggtitle("Pollen shed (MALE)") +
  facet_grid(forcats::fct_rev(Site) ~ MAT_label) +
  scale_y_continuous(
    breaks = seq(1, 365, by = 14),  # Breaks every 2 weeks
    labels = format(seq(as.Date("2023-01-01"), as.Date("2023-12-31"), by = "2 weeks"), "%b %d")) +
  theme(
    axis.text.y = element_text(size = 7),
    strip.text.y = element_text(size = 8),
    legend.position = "bottom"
  )

ggsave(here::here("output/figures/orchpred_doy_male.png"), width = 6, height = 10)

# orchfemale - female only #########
ggplot() +
  geom_line(data = filter(phenf_orchplot, Sex == "FEMALE"), aes(x = Year, y = DoY, group = Year), alpha = 0.9) +
  geom_ribbon(data = filter(doy_annual_avg_pp_sum, Sex == "FEMALE"), aes(x = Year, ymin = .lower, ymax = .upper, group = event, fill = event), alpha = 0.3) +
  geom_line(data = filter(doy_annual_avg_pp_sum, Sex == "FEMALE"), aes(x = Year, y = DoY, colour = event)) +
  scale_fill_discrete_c4a_div(palette = "icefire") +
  scale_colour_discrete_c4a_div(palette = "icefire") +
  theme_bw() +
  xlab("Year") +
  ylab("Date") +
  ggtitle("Receptivity (FEMALE)") +
  facet_grid(forcats::fct_rev(Site) ~ MAT_label) +
  scale_y_continuous(
    breaks = seq(1, 365, by = 14),  # Breaks every 2 weeks
    labels = format(seq(as.Date("2023-01-01"), as.Date("2023-12-31"), by = "2 weeks"), "%b %d")) +
  theme(
    axis.text.y = element_text(size = 7),
    strip.text.y = element_text(size = 8),
    legend.position = "bottom"
  )
ggsave(here::here("output/figures/orchpred_doy_female.png"), width = 6, height = 10)

## orchpred - kalamalka and pgtis only ############

doy_annual_avg_pp_sum_kalpgtis <- doy_annual_avg_pp_sum %>%
  filter(Site %in% c("Kalamalka", "PGTIS")) %>%
  droplevels()
phenf_orchplot_kalpgtis <- phenf_orchplot %>%
  filter(Site %in% c("Kalamalka", "PGTIS")) %>%
  droplevels()
phenf_orchplot_kalpgtis_expanded <- rbind(
  transform(phenf_orchplot_kalpgtis, MAT_label = "Provenance MAT: -0.7 °C"),
  transform(phenf_orchplot_kalpgtis, MAT_label = "Provenance MAT: 6.8 °C")
)


ggplot() +
  geom_ribbon(data = doy_annual_avg_pp_sum_kalpgtis,
              aes(x = Year, ymin = .lower, ymax = .upper, group = event, fill = event), alpha = 0.3) +
  geom_line(data = doy_annual_avg_pp_sum_kalpgtis,
            aes(x = Year, y = DoY, colour = event)) +
  geom_line(data = phenf_orchplot_kalpgtis_expanded,
            aes(x = Year, y = DoY, group = Year), alpha = 0.9, linewidth = .2) +
  scale_fill_discrete_c4a_div(palette = "icefire") +
  scale_colour_discrete_c4a_div(palette = "icefire") +
  theme_bw(base_size = 7) +
  xlab("Year") +
  ylab("Date") +
  facet_nested(Site + MAT_label ~ Sex) +
  scale_y_continuous(
    breaks = seq(1, 365, by = 14),
    labels = format(seq(as.Date("2023-01-01"), as.Date("2023-12-31"), by = "2 weeks"), "%b %d")) +
  theme(
    legend.position = "bottom"
  ) +
  coord_flip()

ggsave(here::here("output/figures/orchpred_doy_vertical.png"), width = 5, height = 8)

### horizontal versions for presentation ###############

pgtisorch <- ggplot() +
    geom_line(data = filter(phenf_orchplot, Site == "PGTIS"), aes(x = Year, y = DoY, group = Year), alpha = 0.9) +
    geom_ribbon(data = filter(doy_annual_avg_pp_sum, Site == "PGTIS"), aes(x = Year, ymin = .lower, ymax = .upper, group = event, fill = event), alpha = 0.3) +
    geom_line(data = filter(doy_annual_avg_pp_sum, Site == "PGTIS"), aes(x = Year, y = DoY, colour = event)) +
    scale_fill_discrete_c4a_div(palette = "icefire") +
    scale_colour_discrete_c4a_div(palette = "icefire") +
    theme_bw(base_size = 8) +
    xlab("Year") +
    ylab("Date") +
    ggtitle("PGTIS", subtitle = "1961-1990 normal MAT: 3.9 \u00B0C") +
    facet_grid(Sex ~ MAT_label) +
    scale_y_continuous(
        breaks = seq(1, 365, by = 14),  # Breaks every 2 weeks
        labels = format(seq(as.Date("2023-01-01"), as.Date("2023-12-31"), by = "2 weeks"), "%b %d")) +
    theme(
        axis.text.y = element_text(size = 8),
        strip.text.y = element_text(size = 9),
        legend.position = "bottom"
    )

kalorch <- ggplot() +
    geom_line(data = filter(phenf_orchplot, Site == "Kalamalka"), aes(x = Year, y = DoY, group = Year), alpha = 0.9) +
    geom_ribbon(data = filter(doy_annual_avg_pp_sum, Site == "Kalamalka"), aes(x = Year, ymin = .lower, ymax = .upper, group = event, fill = event), alpha = 0.3) +
    geom_line(data = filter(doy_annual_avg_pp_sum, Site == "Kalamalka"), aes(x = Year, y = DoY, colour = event)) +
    scale_fill_discrete_c4a_div(palette = "icefire") +
    scale_colour_discrete_c4a_div(palette = "icefire") +
    theme_bw(base_size = 8) +
    xlab("Year") +
    ylab("Date") +
    ggtitle("Kalamalka", subtitle = "1961-1990 normal MAT: 8.0 \u00B0C") +
    facet_grid(Sex ~ MAT_label) +
    scale_y_continuous(
        breaks = seq(1, 365, by = 14),  # Breaks every 2 weeks
        labels = format(seq(as.Date("2023-01-01"), as.Date("2023-12-31"), by = "2 weeks"), "%b %d")) +
    theme(
        axis.text.y = element_text(size = 8),
        strip.text.y = element_text(size = 9),
        legend.position = "bottom"
    )

kalorch / pgtisorch +
    plot_annotation(tag_levels = 'A') +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

ggsave(here::here("output/figures/orchpred_doy.png"), width = 6, height = 7)
# profdiffs ####
# differences between provenances in flowering at different sites #mean difference between coldest and warmest provenance
warmvscold <- readRDS(here::here('output/warmvscold.rds'))
ggplot(warmvscold, aes(x = Site, y = mean_doy_diff, colour = Sex, shape = event)) +
  geom_point(position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = mean_doy_diff - sd_diff, ymax = mean_doy_diff + sd_diff),
                width = 0.2,   # Controls the width of the horizontal lines at the top and bottom of the error bars
                position = position_dodge(0.4)) +
  scale_colour_discrete_c4a_div(palette = "acadia") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=40, vjust=1, hjust=1)) +
  ylab("Mean difference (days)") +
  theme(legend.position = "bottom")
ggsave(here::here("output/figures/provdiffdoy.png"), width = 5.5, height = 3.25)

# distrankcorr ####
# similar patterns of flowering across sites
rank_correlation_wdist <- readRDS(here::here('tmp/rank_correlation_wdist.rds'))

ggplot(filter(rank_correlation_wdist, Correlation < 1), aes(x = Distance, y = Correlation)) +
geom_smooth(method = "lm", se = TRUE, color = "grey25", fill = "grey") +
  geom_point(alpha = 0.5, pch = 1) +
  facet_grid(Sex ~ event) +
  scale_colour_discrete_c4a_div(palette = "icefire") +
  xlab("Distance (km)") +
  ylab(expression("Correlation (Kendall's" ~ tau ~")")) +
  theme_bw() +
  theme(legend.position = "top")
ggsave(here::here("output/figures/distrankcorr.png"), width = 5, height = 4)

corr_model_results <- readRDS(here::here('tmp/corr_model_results.rds'))

# corrmodeltable - distance (slope) only ####
corr_model_table <- corr_model_results %>%
  filter(term == "Distance") %>%
  select(-term, -p.value) %>%
  mutate(across(c(r.squared:std.error), \(x) signif(x, digits = 2)))
saveRDS(corr_model_table, here::here('output/tables/corr_model_table.rds'))

## homeaway - Home vs away #######

#Change in flowering day of year expectation with MAT effect, typical year, trees grown at PGTIS
doy_typical_all_at_PGTIS <- readRDS(here::here("tmp/doy_typical_all_at_PGTIS.rds")) %>%
  pivot_longer(cols = c(intercept, DoY), names_to = "proveffect", values_to = "DoY")

pgtis_intercepts <- doy_typical_all_at_PGTIS %>%
  filter(proveffect == "intercept") %>%
  select(Sex, event, DoY) %>%
  distinct()

awayplot <- ggplot(doy_typical_all_at_PGTIS, aes(x = DoY, y=MAT, colour = proveffect, group = MAT)) +
  geom_point() +
  geom_vline(data = pgtis_intercepts, aes(xintercept = DoY), colour = "darkgrey", linetype = 3) +
  geom_line(colour = 'darkgrey') +
  facet_grid(Sex ~ event) +
  labs(x = "Day of Year",
       y = "Provenance MAT (\u00B0C)",
       colour = "") +  # Removes the title of the colour legend
  scale_colour_manual(values = c("#1B9E77", "darkgrey"),
                      labels = c("With provenance effect", "No provenance effect")) +
  theme_bw() +
  theme(legend.position = "bottom")
awayplot
ggsave(here::here("output/figures/away.png"), width = 5, height = 4)

# When all sources are grown at the same Site (PGTIS), MAT effect reduces overlap, increases differences between provenances

# normalphenology - climate change prediction ####
doy_normal_subset <- readRDS(here::here("output/doy_normal_subset.rds"))
doy_normal_subset$MATlabel <- paste(doy_normal_subset$Site,
                                      " (",
                                      doy_normal_subset$MAT,
                                      "\u00B0C",
                                      ")",
                                      sep = "")


historicalonly <- doy_normal_subset %>% filter(scenario == "historical") %>%
  rename(type = scenario) %>%
  merge(data.frame(scenario = c("ssp245", "ssp585")))
doy_normal_plotting <- doy_normal_subset %>%
  filter(scenario != "historical") %>%
  mutate(type = "future") %>%
  full_join(historicalonly) %>%
  arrange(desc(MAT)) %>%
  mutate(MATlabel = factor(MATlabel, levels = unique(MATlabel), ordered = TRUE))

#convert DoY to date labels
doy_to_date <- function(doy) {
  # Assuming January 1st is Day 1
  date_labels <- format(as.Date(doy - 1, origin = "2024-01-01"), "%b %d")
  return(date_labels)
}

ggplot(filter(doy_normal_plotting, event == "begin"), aes(x = period, y = DoY, colour = Sex, shape = event)) +
    stat_pointinterval(position = position_dodge(width = 1), alpha = 0.5, .width = c(0.50, 0.95)) +
    stat_pointinterval(data = filter(doy_normal_plotting, event == "end"), position = position_dodge(width = 1), alpha = 0.5, .width = c(0.50, 0.95)) +
    facet_grid(scenario ~ MATlabel, labeller = labeller(scenario = c(ssp245 = "SSP2 4.5 W/m²", ssp585 = "SSP5 8.5 W/m²"))) +
    scale_colour_discrete_c4a_div(palette = "acadia") +
    scale_y_continuous(labels = doy_to_date, breaks = seq(100, 200, by = 14)) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.title.x = element_text(vjust = -1)) +
    ylab(NULL) +
    xlab("Normal period") +
    guides(shape = guide_legend(override.aes = list(size = 3)))

ggsave(here::here("output/figures/normal_predictions.png"), width = 9, height = 6, units = "in")




