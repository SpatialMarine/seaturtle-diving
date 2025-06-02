


####### create map of states for the differenrt Sea-turtles ########  ------------------------------------------

library(marmap)

library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(tidyterra)
library(raster)
library(terra)

library(plotly)
library(sf)
library(ggshadow)
library(ggforce)
library(giscoR)
library(ggtext)

source("setup.R")

# load data:
depths_bind <-read.csv("C:/Users/J. Menéndez Blázquez/SML_Dropbox/SML Dropbox/gitdata/seaturtle-diving/output/jess/plots/dive_summariesLabStates060224_FIN1304.csv")

# prepare unique turtle data states
newdf <- depths_bind %>% filter(depths_bind$maximumDepth > 10)
rm(depths_bind)

colnames(newdf)[colnames(newdf) == "latitude"] ="Latitude"
colnames(newdf)[colnames(newdf) == "longitude"] ="Longitude"

ids <- unique(newdf$organismID)

# Prepare spatial layers -------------------------------------------------------
# import landmask and coastline
world <- ne_countries(scale = 10, returnclass = "sf")
coastline <- gisco_get_coastallines(year = "2016", epsg = "4326", resolution = "03")
# transform coastline in the same crs that landmask
coastline <- st_transform(coastline, crs = st_crs(world))


# load bathymetry GEBCO2020
b <- rast(paste0(input_dir,"/gis/GEBCO_2020_Mediterranean_bathymetry.tif"))
# use only values == or < 0 as a bathymetry
b[b > 0] <- NA #0 data as NA

# color ramp for bathymetry
cols <- colorRampPalette(rev(c('#ecf9ff','#BFEFFF','#97C8EB','#4682B4','#264e76','#162e46')))(100)
cols <- adjustcolor(cols, alpha.f = 0.65) 

# limit represent in the plot
xlim <- c(-5.2, 15.4)
ylim <- c(35.2, 44.2)

# Obtain values of batyhymetry in the plot area xlim and ylim
visible_data <- crop(b, extent(c(xlim, ylim)))
visible_range <- range(values(visible_data), na.rm = TRUE)

visible_range <- visible_range + (visible_range * 0.10) # add 5% of the due the expand = T limite in the plot
# or use expand limit = F and addjust ylim and xlim


# create a progress bar for loop in to differet ptt or IDs
pb <- txtProgressBar(min = 0,
                     max = length(ids),
                     style = 3,
                     width = length(ids), # avoid differents prints
                     char = "=") 


t <- Sys.time()
for (i in seq_along(ids)) {
  
  # extact id
  id <- ids[i]
  # filter by id
  filt_df <-newdf[newdf$organismID == id,]
  
  # filter first and last track position from SDMs
  locs <- read.csv(paste0(output_dir,"/tracking/locdata/L2_loc/",id,"_L2_loc.csv"))
  
  start <- head(locs, 1)
  end  <- tail(locs, 1)
  
  
  # Plot map ----------------------------
  p <- ggplot() +
    # add bathymetry
    tidyterra::geom_spatraster(data = b) +
    # color ramp
    scale_fill_gradientn(colors = cols,
                         name = "Depth (m)",
                         limits = visible_range, # limits of values in the represented area
                         guide = guide_colorbar(frame.colour = "grey5", ticks.colour = "grey5"),
                         na.value = "#FFFFFF") +
    # sea-turtle track
    geom_path(data = filt_df,
              aes(x = Longitude, y = Latitude),
              size = 0.75, colour = "grey15") +
    # states values
    geom_point(aes(x=Longitude, y=Latitude, col=as.factor(states4)),
               data=filt_df,
               size = 2.55,
               color = "grey25",
               alpha = 0.1) + 
    geom_point(aes(x=Longitude, y=Latitude, col=as.factor(states4)),
               data=filt_df,
               size = 2.35,
               color = "grey10",
               alpha = 0.8) + 
    geom_point(aes(x=Longitude, y=Latitude, col=as.factor(states4)), #in aes, shape=as.factor(class)
               data=filt_df,
               size = 1.6) + 
    scale_color_manual(name = "State", labels = c("1", "2", "3","4"), values = c("#ca7dcc","#E69F00","#56B4E9", "coral")) +

    # add land and coastline
    geom_sf(data = world, fill="grey98", colour = "grey85", size = .03) +
    geom_sf(data = coastline, fill = "transparent", colour = "grey60", size = 3, alpha = 0) +
    
    # first track position registered
    geom_point(data = start, aes(x = longitude, y = latitude),
               color = "grey20",
               shape = 19,
               alpha = 0.4,
               size = 5) +
    geom_point(data = start, aes(x = longitude, y = latitude),
               fill = "#B3EE3A",
               color = "#000000",
               shape = 21,
               size = 4) +
    
    # last track position registered
    geom_point(data = end, aes(x = longitude, y = latitude),
               color = "grey20",
               shape = 17,
               alpha = 0.4,
               size = 5.1) +
    geom_point(data = end, aes(x = longitude, y = latitude),
               fill = "#FF6347",
               color = "#000000",
               shape = 24,
               size = 4) +
    
    # add id annotate
    # annotate("text", x = -5, y = 44.05, 
    #         label = paste("ptt:", bold(id)),
    #         parse = TRUE, size = 5, color = "black", fill = "white") +
    
    geom_label(aes(x = -2.3, y = 43.6, label = paste("ptt:", id)), 
               size = 8, color = "black", fill = "white", fontface = "bold") +
  
    # spatial bounds
    # coord_sf(xlim = xl, ylim = yl, expand= TRUE) +
    coord_sf(xlim = xlim, ylim = ylim, expand=T) +
    
    # x y labels
    xlab("") +  ylab("") +
    # theme
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_line(size = 0.75),
          axis.ticks.length = unit(6, "pt"),
          panel.border = element_rect(color = "black", fill = NA, size = 1.2),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none")
          
  # export map
  p_png <- paste0(output_dir,"/fig/ai/sup_map4/ptt/",id,"_sup_map4.png")
  # p_svg <- paste0(output_dir,"/fig/ai/sup_map4/ptt/",id,"_sup_map4.svg")
  ggsave(p_png, p, width=20, height=12, units="cm", dpi=350, bg="white")
  #ggsave(p_svg, p, width=20, height=12, units="cm", dpi=350, bg="white")
  
  print(paste0("id: ",id," completed"))
  
  # progress bar update --------
  setTxtProgressBar(pb, i)
}

Sys.time() - t
close(pb) #close progres bar












# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Plot id with bast extension
# -----------------------------------------------------------------------------

# select ids of interest

#------------------------------------------------------------------------------
# Seaturtuele which reached east-med -----------------------------------
id <- "151934"
filt_df <-newdf[newdf$organismID == id,]

# limit represent in the plot for this turtle (mediterranean area)
xlim <- c(-4, 35)
ylim <- c(30.5, 45.5)


# Seaturtule id for west-africa -----------------------------------
id <- "200045"
filt_df <-newdf[newdf$organismID == id,]

# limit represent in the plot for West-Africa
xlim <- c(-25.5, 11)
ylim <- c(13.5, 41)


# prepare spatial layers for plot Med and West-Africa --------------------------
# import landmask and coastline
world <- ne_countries(scale = 10, returnclass = "sf")

# load bathymetry GEBCO2020
b1 <- rast(paste0(input_dir,"/gis/GEBCO_2020_REDUCE_bathymetry.tif"))
# use only values == or < 0 as a bathymetry
b1[b1 > 0] <- NA #0 data as NA

# bathymetry for mediterrenanean sea
b2 <- rast(paste0(input_dir,"/gis/GEBCO_2020_Mediterranean_bathymetry.tif"))
# use only values == or < 0 as a bathymetry
b2[b2 > 0] <- NA #0 data as NA

b <- mosaic(b1, b2, fun = mean)

# check mosaic
plot(b)

# -------------------------------- Plot the ID of interest ---------------------

# color ramp for bathymetry
cols <- colorRampPalette(rev(c('#ecf9ff','#BFEFFF','#97C8EB','#4682B4','#264e76','#162e46')))(100)
cols <- adjustcolor(cols, alpha.f = 0.7) 

# Obtain values of batyhymetry in the plot area xlim and ylim
visible_data <- crop(b, extent(c(xlim, ylim)))
visible_range <- range(values(visible_data), na.rm = TRUE)

visible_range <- visible_range + (visible_range * 0.10) # add 5% of the due the expand = T limite in the plot
# or use expand limit = F and addjust ylim and xlim

# filter first and last track position from SDMs
locs <- read.csv(paste0(output_dir,"/tracking/locdata/L2_loc/",id,"_L2_loc.csv"))

start <- head(locs, 1)
end  <- tail(locs, 1)


# Plot map ----------------------------
p <- ggplot() +
  
  # add bathymetry
  tidyterra::geom_spatraster(data = b) +
  # color ramp
  scale_fill_gradientn(colors = cols,
                       name = "Depth (m)",
                       limits = visible_range, # limits of values in the represented area
                       guide = guide_colorbar(frame.colour = "grey5", ticks.colour = "grey5"),
                       na.value = "#FFFFFF") +
  
  # sea-turtle track 
  # geom_sf(data = track_line, colour = "grey15", size = 1) +
  
  geom_path(data = filt_df,
            aes(x = Longitude, y = Latitude),
            size = 0.75, colour = "grey15") +
  
  # states values
  geom_point(aes(x=Longitude, y=Latitude, col=as.factor(states4)),
             data=filt_df,
             size = 2.55,
             color = "grey25",
             alpha = 0.1) + 
  geom_point(aes(x=Longitude, y=Latitude, col=as.factor(states4)),
             data=filt_df,
             size = 2.35,
             color = "grey10",
             alpha = 0.8) + 
  geom_point(aes(x=Longitude, y=Latitude, col=as.factor(states4)), #in aes, shape=as.factor(class)
             data=filt_df,
             size = 1.6) + 
  scale_color_manual(name = "State", labels = c("1", "2", "3","4"), values = c("#ca7dcc","#E69F00","#56B4E9", "coral")) +
  
  # first and last position registered
  geom_point(data = start, aes(x = longitude, y = latitude),
             color = "grey20",
             shape = 19,
             alpha = 0.4,
             size = 5) +
  geom_point(data = start, aes(x = longitude, y = latitude),
             fill = "#B3EE3A",
             color = "#000000",
             shape = 21,
             size = 4) +
  
  # first and last position registered
  geom_point(data = end, aes(x = longitude, y = latitude),
             color = "grey20",
             shape = 17,
             alpha = 0.4,
             size = 5.1) +
  geom_point(data = end, aes(x = longitude, y = latitude),
             fill = "#FF6347",
             color = "#000000",
             shape = 24,
             size = 4) +
  
  # add land and coastline
  geom_sf(data = world, fill="grey98", colour = "grey85", size = .03) +
  geom_sf(data = coastline, fill = "transparent", colour = "grey60", size = 3, alpha = 0) +
  
  # med = 1.8, 43.9
  #from atlantci x = 5 y = 15)
  
  geom_label(aes(x = 1.8, y = 43.9, label = paste("ptt:", id)), 
             size = 6, color = "black", fill = "white", fontface = "bold") +
  
  # spatial bounds
  # coord_sf(xlim = xl, ylim = yl, expand= TRUE) +
  coord_sf(xlim = xlim, ylim = ylim, expand = T) +
  # x y labels
  xlab("") +  ylab("") +
  
  # theme
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.ticks = element_line(size = 0.7),
        axis.ticks.length = unit(6, "pt"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.2),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        legend.justification = "center",
        legend.key.width = unit(15, "pt"),
        legend.key.height = unit(20, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8))

p


# export map
p_png <- paste0(output_dir,"/fig/ai/sup_map4/ptt/",id,"_sup_map4.png")
p_svg <- paste0(output_dir,"/fig/ai/sup_map4/ptt/",id,"_sup_map4.svg")
ggsave(p_png, p, width=20, height=12, units="cm", dpi=400, bg="white")
ggsave(p_svg, p, width=20, height=12, units="cm", dpi=400, bg="white")
