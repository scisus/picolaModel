library(patchwork)

pointmap

growingseason_sumforcplot <- ggplot(typical_year_forc, aes(x = Date, y = sum_forcing, color = Site, linetype = Site_type)) +
    geom_line() +
    scale_x_date(limits = as.Date(c("2024-03-15", max(typical_year_forc$Date))), date_breaks = "1 month", date_labels =  "%b") +
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

growingseason_sumforcplot

provs <- readRDS(here::here("output/phenf.rds")) %>%
    select(Tree, MAT, Site, Genotype) %>%
    distinct() %>%
    left_join(replication_points) %>%
    rename('Within Sites' = treestf, 'Across Sites' = sitestf, 'Across Years' = yearstf, Replicated = replicated) %>%
    mutate(Site_plot = factor(Site, levels = site_levels))


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


siteprovplot <- ggplot(provs, aes(x = Site, y = MAT)) +
    geom_quasirandom(varwidth = TRUE, alpha = 0.5, shape = 3) +
    scale_y_continuous(limits = c(min(provs$MAT), max(provs$MAT)), position = "right") +
    ylab("Mean Annual Temperature (\u00B0C)") +
    theme(legend.position = "bottom") +
    ggtitle("Provenance MATs")
siteprovplot

site_levels <- factororder$site %>%
    if_else(. %in% c("Border", "Trench"), "Comparison", .) %>%
    unique()

site_mat_lines <- sites %>%
    mutate(
        Site_plot = if_else(`Site Type` == "Seed Orchard", Site, "Comparison"),
        Site_plot = factor(Site_plot, levels = site_levels)
    )

comparison_labels <- site_mat_lines %>%
    filter(`Site Type` == "Comparison")

siteprovplot <- ggplot(provs, aes(x = Site_plot, y = MAT)) +
    geom_quasirandom(
        varwidth = TRUE,
        alpha = 0.7,
        shape = 3
    ) +
    geom_errorbar(
        data = site_mat_lines,
        aes(x = Site_plot, ymin = MAT, ymax = MAT),
        inherit.aes = FALSE,
        width = 0.65,
        linewidth = 1.1,
        colour = "black"
    ) +
    geom_text(
        data = comparison_labels,
        aes(
            x = Site_plot,
            y = MAT,
            label = Site
        ),
        inherit.aes = FALSE,
       # nudge_x = 0.18,
        nudge_y = 0.25,
        hjust = 0.5,
        vjust = 0,
        size = 3
    ) +
    scale_x_discrete(
        limits = site_levels,
        expand = expansion(add = c(0.4, 1.1))
    ) +
    scale_y_continuous(
        limits = range(c(provs$MAT, sites$MAT), na.rm = TRUE)
    ) +
    coord_cartesian(clip = "off") +
    ylab("Mean annual temperature (°C)") +
    xlab("Site") +
    theme_bw() +
    theme(
        legend.position = "bottom",
        plot.margin = margin(5, 35, 5, 5)
    )
siteprovplot
