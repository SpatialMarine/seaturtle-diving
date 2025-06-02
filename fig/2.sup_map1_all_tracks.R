# Suplemmentary Figure - All track seaturtle distribution


library(dplyr)
library(data.table)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)

library(tidyterra)
library(raster)
library(marmap)
library(terra)

library(plotly)
library(sf)
library(ggshadow)
library(ggforce)
library(giscoR)

source("fun/reading_tools.R")
source("setup.R")


# -----------------------------------------------------------------------------
# Mediterranean and Wes-Africa zoomed figures ---------------------------------

# Load SSM results
# Set path to tracking data
trackpath <- "C:/Users/J. Menéndez Blázquez/SML_Dropbox/SML Dropbox/gitdata/seaturtle-diving/output/tracking/locdata/L2_loc"

# search files
track_files <- tibble(
  file = list.files(trackpath, pattern = "loc.csv", recursive = TRUE, full.names = TRUE))

# batch import
data <- rbindlist(lapply(track_files$file, fread), fill=TRUE)

# rename variables
data <- data %>%
  rename(OrganismID = organismID, 
         trip = tripID)

# import landmask and coastline
world <- ne_countries(scale = "medium", returnclass = "sf")
coastline <- gisco_get_coastallines(year = "2016", epsg = "4326", resolution = "03")
# transform coastline in the same crs that landmask
coastline <- st_transform(coastline, crs = st_crs(world))


# load bathymetry GEBCO2020
b1 <- rast(paste0(input_dir,"/gis/GEBCO_2020_REDUCE_bathymetry.tif"))
# use only values == or < 0 as a bathymetry
b1[b1 > 0] <- NA #0 data as NA

# bathymetry for mediterrenanean sea
b2 <- rast(paste0(input_dir,"/gis/GEBCO_2020_Mediterranean_bathymetry.tif"))
# use only values == or < 0 as a bathymetry
b2[b2 > 0] <- NA #0 data as NA

b <- mosaic(b1, b2, fun = mean)



# color ramp for bathymetry
cols <- colorRampPalette(rev(c('#ecf9ff','#BFEFFF','#97C8EB','#4682B4','#264e76','#162e46')))(100)
cols <- adjustcolor(cols, alpha.f = 0.75) 

# limit represent in the plot for Mediterranean area
xlim <- c(-25.5, 35)
ylim <- c(13.5, 45.2)


# Obtain values of batyhymetry in the plot area xlim and ylim
visible_data <- crop(b, extent(c(xlim, ylim)))
visible_range <- range(values(visible_data), na.rm = TRUE)

# Plot map ----------------------------
p <- ggplot() +
  
  # add bathymetry
  tidyterra::geom_spatraster(data = b) +
  # bathymettry color ramp
  scale_fill_gradientn(colors = cols,
                       name = "Depth (m)",
                       limits = visible_range, # limits of values in the represented area
                       guide = guide_colorbar(frame.colour = "grey5", ticks.colour = "grey5"),
                       na.value = "#FFFFFF") +
  
  # add tracks 
  geom_path(data = data,
            aes(x = longitude, y = latitude, group = interaction(OrganismID, trip)),
            size = 1.7, alpha = 0.25, colour = "tan1") +
  geom_path(data = data,
            aes(x = longitude, y = latitude, group = interaction(OrganismID, trip)),
            size = 1.2, alpha = 0.35, colour = "darkorange") +
  geom_path(data = data,
            aes(x = longitude, y = latitude, group = interaction(OrganismID, trip)),
            size = 0.2, alpha = 0.95, colour = "darkorange4") +
  
  # add land and coastline
  geom_sf(data = world, fill="grey98", colour = "grey85", size = .03) +
  geom_sf(data = coastline, fill = "transparent", colour = "grey60", size = 3, alpha = 0) +

  # spatial bounds
  coord_sf(xlim = xlim, ylim = ylim, expand=T) +

  # x y labels
  xlab("") +  ylab("") +
  # theme
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.ticks = element_line(size = 0.75),
        axis.ticks.length = unit(6, "pt"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.2),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.justification = "center",
        legend.key.width = unit(15, "pt"),
        legend.key.height = unit(20, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8))

p


# save plot
p_png <- paste0(output_dir, "/fig/ai/sup_map1/sup_map.png")
p_svg <- paste0(output_dir, "/fig/ai/sup_map1/sup_map.svg")
ggsave(p_png, p, width=20, height=12, units="cm", dpi=450, bg="white")
ggsave(p_svg, p, width=20, height=12, units="cm", dpi=450, bg="white")




