here::i_am("output/6.map.R")

library(ggplot2)
theme_set(theme_bw())
library(sf)
library(cols4all)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(viridis)
library(dplyr)
library(ggspatial)
library(terra)
library(elevatr)
library(ggrepel)

# NAD83 / Canada Atlas Lambert 3978
# WGS 84 -- WGS84 - World Geodetic System 1984, used in GPS 4326

# data ########
# country outlines
basedat <- ne_states(country = c("canada", "united states of america"), returnclass = "sf") %>%
    st_transform(basedat, crs = 3978)

# terrain

## Create a SpatialPointsDataFrame with the extent coordinates
extent_points <- st_as_sf(data.frame(x = c(-145, -65), y = c(45,70)),
                          coords = c("x", "y"),
                          crs = 4326)
extent_points_zoom <- st_as_sf(data.frame(x = c(-119.5, -119), y = c(50.1,50.4)),
                               coords = c("x", "y"),
                               crs = 4326)

## Download & format the DEM

get_elevations <- function(extent, zoom, projection = st_as_text(st_crs(basedat))){
  elevations_dat <- get_elev_raster(extent, z = zoom, clip = "bbox", verbose = FALSE)
  elevations_dat@data@names <- "elevation"


  # Create a SpatRaster from the elevations_dat
  elevations <- rast(elevations_dat) %>%
    terra::project(projection)
  elevations_df <- as.data.frame(elevations, xy = TRUE)
  elevations_df$elevation[elevations_df$elevation<0] <- NA # ignore oceans

  return(elevations_df)

}

elevations <- get_elevations(extent = extent_points, zoom = 4)

# lodgepole pine distribution shapefile
pcontorta <- st_read(here::here("inputs/latifoliaDistribution/shapefiles/latifolia_distribution_prj.shp")) %>%
    st_make_valid()

genotypes <- select(picolaDataFlowering::picola_event, Genotype) %>% distinct()

parents <- picolaDataFlowering::picola_parent_locs %>%
    mutate(Genotype = Clone) %>%
    select(Genotype, lat, lon) %>%
    filter(Genotype %in% genotypes$Genotype) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326, agr = "constant")

# orchard locations
sitedat <- picolaDataFlowering::picola_site_coord_elev

# 4 sites are very close together. create different dfs for mapping labels at different map scales
sitezoomout <- sitedat %>%
  filter(!Site %in% c("Vernon","Tolko","PRT", "Kalamalka")) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, agr = "constant")

sitezoomin <- sitedat %>%
  filter(Site %in% c("Vernon","Tolko","PRT", "Kalamalka")) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, agr = "constant")

sites <- sitedat %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, agr = "constant")

# bboxes ##########

## bounding box for distribution
bboxdist <- st_bbox(pcontorta)

# parent bounding box
bboxparents <- st_bbox(parents) %>%
    st_as_sfc() %>%
    st_transform(crs = 3978) %>%
    st_bbox()

# sites bounding boxes
bboxsites <- st_bbox(sites) %>%
    st_as_sfc() %>%
    st_transform(crs = 3978) %>%
    st_bbox()

bboxsitezoom <- sites %>%
    filter(Site %in% c("Kalamalka", "Vernon", "PRT", "Tolko")) %>%
    st_bbox() %>%
    st_as_sfc() %>%
    st_transform(crs = 3978) %>%
    st_bbox()

# maps ###########


basemap <- ggplot(data = basedat) +
  geom_raster(data = elevations, aes(x = x, y = y, fill = elevation)) +
  scale_fill_gradientn(colours = grey.colors(20, end = 0.8), na.value = NA) +
  geom_sf(data = pcontorta, alpha = 0.2, fill = "darkolivegreen3") +
  geom_sf(fill = NA) +
  annotation_north_arrow(location = "bl", which_north = "true", # location set to "tl"
                         pad_x = unit(0.05, "in"), pad_y = unit(0.05, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "br", width_hint = 0.22) + # added scale bar
  coord_sf(xlim = c(bboxparents$xmin - 4e5, bboxparents$xmax + 1e5),
           ylim = c(bboxsites$ymin - 1e5, bboxsites$ymax + 3e5)) +
  theme(legend.position = "none") +
  ylab("") + xlab("")

print(basemap)

pointmap <- basemap +
  geom_sf(data = sites, size = 2, aes(shape = orchard, color = orchard)) +
  scale_shape_manual(values = c('FALSE' = 17, 'TRUE' = 16)) +
  scale_colour_brewer(type = "qual", palette = "PuOr") +
  geom_sf(data = parents, shape = 3, alpha = 0.8) +
  coord_sf(xlim = c(bboxparents$xmin - 3e5, bboxparents$xmax + 2e5),
           ylim = c(bboxsites$ymin - 1e5, bboxsites$ymax + 3e5)) +
  geom_label_repel(data = sitezoomout, aes(label = Site, geometry = geometry),
                   stat = "sf_coordinates", nudge_x = 5e5) +
  geom_label_repel(data = sitezoomin, aes(label = Site, geometry = geometry),
                   stat = "sf_coordinates", nudge_x = -5e5, nudge_y = 5e4)

print(pointmap)

# Save the plot
ggsave(filename = here::here("output/figures/siteandparentmap.png"), plot = pointmap, width = 7, height = 7, dpi = 300, units = "in")
