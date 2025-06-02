

# seaturtle-diving 
# updating the figures and script 
# by J.Menéndez-Blázquez (August, 2024)



# load packages
rm(list=ls())

library(depmixS4)
library(ggplot2)
library(tidyr)
library(dplyr)
library(caret)
library(egg)
library(grid)
library(ggdist)


# load data:
depths_bind <- read.csv("C:/Users/J. Menéndez Blázquez/SML_Dropbox/SML Dropbox/gitdata/seaturtle-diving/output/jess/plots/dive_summariesLabStates060224_FIN1304.csv")

toplot<-subset(depths_bind, select = c(organismID,diveID, states3, states4,states5))


#randomly select n dives/ state for plotting
set.seed(9)
toplot2 <- toplot %>% group_by(states4) %>% slice_sample(n=200)

#read raw TDR data
tot_dives<-read.csv("C:/Users/J. Menéndez Blázquez/SML_Dropbox/SML Dropbox/gitdata/seaturtle-diving/output/jess/plots/allDiveDepths.csv")

#label each dataframe concatenated organsimID and dive number
tot_dives$newid<-paste0(tot_dives$organismID,"_",tot_dives$diveno)

toplot2$newid<-paste0(toplot2$organismID,"_",toplot2$diveID)


#merge dataframes to include state information
mergedf<-merge(toplot2, tot_dives, by= "newid")

#subset randomly selected dives
subDives<-subset(mergedf, newid %in% toplot2$newid)


subDives<-subDives[order(subDives$X.1),]
subDives$date <- as.POSIXct(subDives$date,format="%Y-%m-%d %H:%M:%S",tz="CET")


#make sure TDR readings are in correct order
subDives2<-subDives %>%
  group_by((newid)) %>%
  mutate(numbering= row_number())


df2 <- subDives2[order(subDives2$X.1),]

subDives2<-as.data.frame(subDives2)


#flag start of each dive
subDives2<-subDives2 %>%
  mutate(t = case_when(numbering == 1 ~ '0',
                       TRUE ~ ''))

##find difference between each TDR reading in secs to plot
subDives3<-subDives2 %>%
  group_by(newid) %>%
  mutate(idletime = difftime(date, date[1], units = "sec"))
as.data.frame(subDives3)

subDives3$idletime<-as.numeric(subDives3$idletime)


#use this function to cleanup graphs

cleanup=theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background=element_blank(),
              axis.line=element_line(color = "black"))

# names for the plot
names <- c(
  `1` = "State 1",
  `2` = "State 2",
  `3` = "State 3",
  `4` = "State 4"
)

# ------------------------------------------------------------------------------
# Fig 2. -----------------------------------------------------------------------
#plot dives for each states in facets # ----------------------------------------


p <- ggplot(subDives3, aes(x = idletime, y =-depth, color = as.factor(subDives3$states4)), alpha = 0.25) +
        # states lines
        geom_line(aes(group=as.factor(diveID)), alpha = 0.25, linewidth = 0.7) +
        facet_grid(~ states4,labeller = as_labeller(names)) +
        # axys
        xlab("Time (Seconds)") +
        ylab("Depth (m)") +
        scale_colour_manual(name = "State", values=c("#ca7dcc","#E69F00","#56B4E9", "coral")) +
        # theme
        theme_bw() +
        theme(axis.text.y = element_text(size = 10),
              axis.title.y = element_text(margin = margin(r = 8)),
              axis.text.x = element_text(size = 10),
              axis.title.x = element_text(margin = margin(t = 8)),
              axis.ticks = element_line(size = 0.7),
              axis.ticks.length = unit(7.5, "pt"),
              strip.text.x = element_text(size = 12, face="bold"),
              strip.background = element_blank(),
              panel.spacing = unit(2, "lines"),
              panel.border = element_rect(color = "black", fill = NA, size = 1.2),
              panel.background = element_blank(),
              panel.grid = element_blank(),
              legend.position = "none")

p        

# export fig 3a
p_png <- paste0(output_dir,"/fig/fig2.png")
p_svg <- paste0(output_dir,"/fig/fig2.svg")
ggsave(p_png, p, width=30, height=10, units="cm", dpi=300, bg="white")
ggsave(p_svg, p, width=30, height=10, units="cm", dpi=300, bg="white")






# --------------------------------------------------------------------------------
# Figure 3a ----------------------------------------------------------------------
# Link this figure with fig 3b   -------------------------------------------------
#now plot timeseries of dives coloured for each state, like fig. 4 in tenneson et al., 


table(depths_bind$organismID)

depths_bind_plot<- depths_bind%>% filter(organismID=="235396") %>%
  select(organismID,diveID,states4)

tot_dives_plot<- tot_dives%>% filter(organismID=="235396") %>%
  select(organismID,diveno,date,depth)

#label each dataframe concatenated organsimID and dive number
tot_dives_plot$newid<-paste0(tot_dives_plot$organismID,"_",tot_dives_plot$diveno)

depths_bind_plot$newid<-paste0(depths_bind_plot$organismID,"_",depths_bind_plot$diveID)
#merge dataframes to include state information
divegram<-merge(depths_bind_plot, tot_dives_plot, by= "newid")

divegram$date <- as.POSIXct(divegram$date,format="%Y-%m-%d %H:%M:%S",tz="CET")
divegram<-divegram[order(divegram$date),]

divegram[which(divegram$states4=="1"),]$states4 <-"s1"
divegram[which(divegram$states4=="2"),]$states4 <-"s2"
divegram[which(divegram$states4=="3"),]$states4 <-"s3"
divegram[which(divegram$states4=="4"),]$states4 <-"s4"


data_wide <- spread(divegram, states4,depth)
data_wide

# plot
p <- ggplot(data_wide, aes(x = date)) +
  #Subset data for black line, i.e., shared line
  geom_line(aes(y = -s1), data = data_wide,
            color = "#ca7dcc", alpha = 0.4) + 
  geom_line(aes(y = -s2), data = data_wide,
            color = "#E69F00", alpha = 0.4) + 
  geom_line(aes(y = -s3), data = data_wide,
            color = "#56B4E9", alpha = 0.4) + 
  geom_line(aes(y = -s4), data = data_wide,
            color = "coral", alpha = 0.4) +
  scale_colour_manual(name = "State", values=c("#ca7dcc","#E69F00","#56B4E9", "coral")) +
  # axys
  ylim(-divegram$depth[which.max(divegram$depth)], -3)+
  xlab("Date") +
  ylab("Depth (m)")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%m/%Y")+
  # theme
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.ticks = element_line(size = 0.7),
        axis.ticks.length = unit(7.5, "pt"),
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

# export fig 3a
p_png <- paste0(output_dir,"/fig/3a_raw.png")
p_svg <- paste0(output_dir,"/fig/3a_raw.svg")
ggsave(p_png, p, width=20, height=14, units="cm", dpi=300, bg="white")
ggsave(p_svg, p, width=20, height=14, units="cm", dpi=300, bg="white")
# ------------------------------------------------------------------------------





### time budget plots ##
# ------------------------------------------------------------------------------
# Fig 5a -----------------------------------------------------------------------
# boxplot - Proportion of total dives per state   ------------------------------

##to match teneson papwer
p <- ggplot(freq.df2,# Draw barplot with grouping & stacking
             aes(x = state ,
                 y = prop,
                 fill = state)) + 
        ggdist::stat_halfeye(alpha = 0.5, 
                             adjust = 0.7,  # adjust bandwidth
                             justification = 0, # move to the right
                             .width = 0,   # remove the slub interval
                             point_colour = NA) +
        geom_boxplot(width = 0.25, size = 0.55, color = "#000000") +
        geom_jitter(width = 0.075, size = 1) +
        labs(fill = "State")+
        xlab("State")+
        ylab("Proportion of Dives")+
        scale_fill_manual(values=c("#ca7dcc", "#E69F00","#56B4E9","coral"))+
        # theme
        theme_bw() +
        theme(axis.text.y = element_text(size = 10),
              axis.title.y = element_text(size = 10, margin = margin(r = 8)),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks = element_line(size = 0.75),
              axis.ticks.x = element_blank(),
              axis.ticks.length = unit(6, "pt"),
              panel.border = element_rect(color = "black", fill = NA, size = 1),
              panel.background = element_blank(),
              panel.grid = element_blank(),
              legend.position = "none")


p

# export fig 5a
p_png <- paste0(output_dir,"/fig/5a.png")
p_svg <- paste0(output_dir,"/fig/5a.svg")
ggsave(p_png, p, width=21, height=12, units="cm", dpi=300, bg="white")
ggsave(p_svg, p, width=21, height=12, units="cm", dpi=300, bg="white")



# Fig. 5b         --------------------------------------------------------------
##### plot % time per dive ###

depths_bind$timeState<-depths_bind$diveTime+depths_bind$surface_int

freq.df<-as.data.frame(table(depths_bind$organismID, depths_bind$states4))

freq.df <- freq.df %>% 
  rename("id"="Var1",
         "state"= "Var2")

freq.df2<-as.data.frame(freq.df %>%
                          group_by(id) %>% #id
                          mutate(tot_dives =sum(Freq)) %>% 
                          ungroup)

freq.df2<-as.data.frame(freq.df2 %>%
                          group_by(id) %>% 
                          mutate(dives_not = tot_dives-Freq))


freq.df2<-as.data.frame(freq.df2 %>%
                          group_by(id) %>% #id
                          mutate(prop =  Freq/sum(Freq)) %>% 
                          ungroup)

freq.df2 <- freq.df2 %>% 
  rename("dives_performed"="Freq")

p <- ggplot(freq.df2,# Draw barplot with grouping & stacking
       aes(x =id ,
           y = prop,
           fill = state)) + 
      geom_bar(stat = "identity",
               position = "stack",
               alpha = 0.85) + 
      labs(fill = "State")+
      xlab("Organism ID")+
      ylab("Proportion of Dives")+
      scale_fill_manual(values=c("#ca7dcc", "#E69F00","#56B4E9","coral"))+
      # theme
      theme_bw() +
      theme(axis.text.y = element_text(size = 10),# angle = 90, hjust = 0.5),
            axis.title.y = element_text(size = 10, margin = margin(r = 8)),
            axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5),
            axis.title.x = element_text(margin = margin(t = 8)),
            axis.ticks = element_line(size = 0.75),
            axis.ticks.length = unit(6, "pt"),
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none")

p

# export fig 5b
p_png <- paste0(output_dir,"/fig/5b.png")
p_svg <- paste0(output_dir,"/fig/5b.svg")
ggsave(p_png, p, width=21, height=12, units="cm", dpi=300, bg="white")
ggsave(p_svg, p, width=21, height=12, units="cm", dpi=300, bg="white")









####### create map of states ########  ------------------------------------------

#' Good idea showing the same in a map. This map would need some improvements:
#'   - Increase font size from axis labels
#' - Land map and coastline are very coarse. Need to replace by higher resolution.
#' - Depth contours and bathymetry seem to put some noise in terms of colors. It may need to change or even remove.
#' - Turtle track: I would add the turtle track using a line, and then the dives on top. Also, the start and end of track can be represented with different symbol.
#' - There is enough space on the top left side (Spain mainland) to add a world map with box of the zoomed area. So people not familiar with study location can identify easily.
#' 
#' @Javi, I suggest that you could help with this figure as you already have R code to do most of these suggestions.

library(marmap)

library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(tidyterra)
# library(raster)
library(terra)

library(plotly)
library(sf)
library(ggshadow)
library(ggforce)
library(giscoR)

# source("setup.R")


# load data:
depths_bind <- read.csv("C:/Users/J. Menéndez Blázquez/SML_Dropbox/SML Dropbox/gitdata/seaturtle-diving/output/jess/plots/dive_summariesLabStates060224_FIN1304.csv")

library(dplyr)
# prepare unique turtle data states
newdf <- depths_bind %>% filter(depths_bind$maximumDepth>10)

colnames(newdf)[colnames(newdf) == "latitude"] ="Latitude"
colnames(newdf)[colnames(newdf) == "longitude"] ="Longitude"

#newdf2<- newdf %>% filter (Longitude>-5)
id <- "235396"
filt_df <- newdf[newdf$organismID == id,]


# prepare spatial layers --------------------------
# import landmask and coastline
world <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf")

# load bathymetry GEBCO2020
b <- rast(paste0(input_dir,"/gis/GEBCO_2020_Mediterranean_bathymetry.tif"))
plot(b)

# use only values == or < 0 as a bathymetry
# b[b > 0] <- NA #0 data as NA

# color ramp for bathymetry
cols <- colorRampPalette(rev(c('#ecf9ff','#BFEFFF','#97C8EB','#4682B4','#264e76','#162e46')))(100)
cols <- adjustcolor(cols, alpha.f = 0.75) 

# limit represent in the plot
xlim <- c(-5, 5)
ylim <- c(34.9, 40)

# Obtain values of batyhymetry in the plot area xlim and ylim
visible_data <- crop(b, raster::extent(c(xlim, ylim)))
visible_range <- range(values(visible_data), na.rm = TRUE)

# sea-turtle track
track_line <- st_read(paste0(input_dir,"/gis/235396_L2_loc_track.gpkg"))

locs <- read.csv(paste0(output_dir,"/tracking/locdata/L2_loc/",id,"_L2_loc.csv"))
# filter first and last track position
start <- head(locs, 1)
end  <- tail(locs, 1)


# Note** Don't load raster() package some issues with terra package.
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
  geom_sf(data = track_line, colour = "grey15", size = 1) +
  
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
  geom_sf(data = world, fill="grey98", colour = "grey35", size = 1) +
  
  # spatial bounds
  # coord_sf(xlim = xl, ylim = yl, expand= TRUE) +
  coord_sf(xlim = xlim, ylim = ylim, expand=T) +
  # x y labels
  xlab("") +  ylab("") +
  # theme
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.ticks = element_line(size = 0.75),
        axis.ticks.length = unit(7.5, "pt"),
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
p_png <- paste0(output_dir,"/fig/3b_raw.png")
p_svg <- paste0(output_dir,"/fig/3b_raw.svg")
ggsave(p_png, p, width=20, height=14, units="cm", dpi=300, bg="white")
ggsave(p_svg, p, width=20, height=14, units="cm", dpi=300, bg="white")



# ------------------------------------------------------ SCINTIFIC REPORTS
# revision fig panel == Figure 3B but split into differenrt maps panel by divign state
# run previous code for seaturtle location, bathymetry, etc.


filt_df$states4 <- as.factor(filt_df$states4)

# Plot panel maps for this figure  ----------------------------
p <- ggplot() +
  
  # p <- ggplot() +
  
  # add bathymetry
  tidyterra::geom_spatraster(data = b) +
  scale_fill_gradientn(colors = cols,
                       name = "Depth (m)",
                       limits = visible_range,
                       guide = guide_colorbar(frame.colour = "grey5", ticks.colour = "grey5"),
                       na.value = "#FFFFFF") +
  
  # sea-turtle track 
  geom_sf(data = track_line, colour = "grey15", size = 1) +
  
  # states values (puntos en capas superpuestas para efecto visual)
  geom_point(aes(x = Longitude, y = Latitude),
             data = filt_df,
             size = 2.50,
             color = "grey25",
             alpha = 0.1) + 
  geom_point(aes(x = Longitude, y = Latitude),
             data = filt_df,
             size = 2.30,
             color = "grey10",
             alpha = 0.8) + 
  geom_point(aes(x = Longitude, y = Latitude, col = states4),
             data = filt_df,
             size = 1.5) + 
  scale_color_manual(name = "State", 
                     labels = c("1", "2", "3", "4"), 
                     values = c("#ca7dcc", "#E69F00", "#56B4E9", "coral")) +
  
  # first position
  geom_point(data = start, aes(x = longitude, y = latitude),
             color = "grey20", shape = 19, alpha = 0.4, size = 5.3) +
  geom_point(data = start, aes(x = longitude, y = latitude),
             fill = "#B3EE3A", color = "#000000", shape = 21, size = 4.2) +
  
  # last position
  geom_point(data = end, aes(x = longitude, y = latitude),
             color = "grey20", shape = 17, alpha = 0.4, size = 5.3) +
  geom_point(data = end, aes(x = longitude, y = latitude),
             fill = "#FF6347", color = "#000000", shape = 24, size = 4.2) +
  
  # land and coastline
  geom_sf(data = world, fill = "grey98", colour = "grey35", size = 1) +
  
  # spatial bounds
  coord_sf(xlim = xlim, ylim = ylim, expand = TRUE) +
  
  # labels y tema
  xlab("") + ylab("") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9.5),
        axis.text.x = element_text(size = 9.5),
        axis.ticks = element_line(size = 0.70),
        axis.ticks.length = unit(7.0, "pt"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.1),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        legend.justification = "center",
        legend.key.width = unit(15, "pt"),
        legend.key.height = unit(20, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        strip.text = element_blank()) +
  
  # facetea por estado
  facet_wrap(~ states4)

p



# export map
p_png <- paste0(output_dir,"/fig/3b_raw_facet.png")
p_svg <- paste0(output_dir,"/fig/3b_raw_facet.svg")
ggsave(p_png, p, width=30, height=20, units="cm", dpi = 350, bg="white")
# ggsave(p_svg, p, width=20, height=14, units="cm", dpi=300, bg="white")














# ------------------------------------------------------------------------------
# Figure 4   ------------------------------------------------------------------
### plot move persistance histograms #####

sum<-subset(depths_bind, select = c(maximumDepth, bottomTime,diveTime, surface_int, mpm, states4))
sum$states4<-as.factor(sum$states4)

facet_names <- list(
  '1'="State 1",
  '2'="State 2",
  '3'="State 3",
  '4'="State 4"
)


name_labeller <- function(variable,value){
  return(facet_names[value])
}

p <- sum %>%
        ggplot(aes(x=mpm, fill=states4)) +
        geom_histogram( color="grey20", position = 'identity') +
        scale_fill_manual(values=c("#ca7dcc","#E69F00","#56B4E9", "coral"),guide="none") +
        #labs(fill="")+
        theme(legend.position = "none")+
        facet_wrap(~states4, labeller=name_labeller) +
        theme_classic()+
        xlab("Move Persistence") +
        ylab("Count")+
        # theme 
        theme_article()+
        theme_bw() +
        theme(axis.text.y = element_text(size = 10),
              axis.text.x = element_text(size = 10),
              axis.ticks = element_line(size = 0.75),
              axis.ticks.length = unit(7, "pt"),
              axis.title.y = element_text(margin = margin(r = 8), size = 11, face = "bold"),
              axis.title.x = element_text(margin = margin(t = 8), size = 11, face = "bold"),
              panel.border = element_rect(color = "black", fill = NA, size = 1.2),
              panel.background = element_blank(),
              panel.grid = element_blank(),
              strip.text.x = element_text(size = 12, face="bold"),
              strip.background = element_blank())

p



# export map
p_png <- paste0(output_dir,"/fig/fig4.png")
p_svg <- paste0(output_dir,"/fig/fig4.svg")
ggsave(p_png, p, width=23, height=17, units="cm", dpi=300, bg="white")
ggsave(p_svg, p, width=23, height=17, units="cm", dpi=300, bg="white")










# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Supplementary Figure 2 -------------------------------------------------------


dive_mpm<-read.csv("C:/Users/J. Menéndez Blázquez/SML_Dropbox/SML Dropbox/gitdata/seaturtle-diving/output/jess/plots/dive_summariesLabStates060224_FIN1304.csv")


sum<-subset(dive_mpm, select = c(maximumDepth, bottomTime, diveTime, surface_int, mpm, states4))
sum$states4<-as.factor(sum$states4)

facet_names <- list(
  '1'="State 1",
  '2'="State 2",
  '3'="State 3",
  '4'="State 4"
)


name_labeller <- function(variable,value){
  return(facet_names[value])
}

library(GGally)

sum %>%
  ggplot(aes(x=mpm, fill=states4)) +
  geom_histogram( color="#e9ecef", position = 'identity') +
  scale_fill_manual(values=c("#ca7dcc","#E69F00","#56B4E9", "coral"),guide="none") +
  #labs(fill="")+
  theme(legend.position = "none")+
  facet_wrap(~states4, labeller=name_labeller)+
  theme_classic()+
  xlab("Move Persistence") +
  ylab("Count")+
  theme_article()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size=14,face="bold"))

names(sum)[names(sum) == 'maximumDepth'] <- 'Maximum Depth'
names(sum)[names(sum) == 'bottomTime'] <- 'Bottom Time'
names(sum)[names(sum) == 'surface_int'] <- 'Surface Interval'
names(sum)[names(sum) == 'diveTime'] <- 'Dive Time'

# pair plot ----------------------------------------------
my_cols <- c("#ca7dcc","#E69F00","#56B4E9", "coral")  

par(cex.axis = 1.5,  # Aumenta el tamaño de las etiquetas de los ejes
    cex.lab = 1.5)  

pairs(sum[,1:4], pch = 19,  cex = 1,
      col = my_cols[sum$states4],
      lower.panel=NULL)











######################################################################################
#validation graphs for supplimentary using oceanografics data

################################################################################
#now compare one second data and resampled 5 min data



# Suplementary material 1 ------------------------------------------------------

output_data <- paste0("C:/Users/xcarrj/Desktop/laptop/tort/")
fiveMin<- read.csv(paste0(output_data, "dive_summarys_fiveMin.csv"))

oneSec<-read.csv(paste0(output_data, "dive_summarys_oneSec.csv"))


#add sample interval for combining
oneSec$sample<-"1 sec"
fiveMin$sample<-"300 sec"
library(dplyr)

oneSec<-oneSec %>% select (diveStartTime, bottomTime,maximumDepth, diveTime,surface_int,diveEndTime,organismID,concatID,sample)
fiveMin<-fiveMin %>% select (diveStartTime, bottomTime,maximumDepth, diveTime,surface_int,diveEndTime,organismID,concatID,sample)

#filter dives with max depth greater than 10m
one_filt<-oneSec %>% filter (maximumDepth>10)
five_filt<-fiveMin %>% filter (maximumDepth>10)

#how manys dives for each sampling interval?
nrow(one_filt)
nrow(five_filt)

#bind dives together
allSamples<-rbind(oneSec,fiveMin)

#frequency table
table(allSamples$sample,allSamples$maximumDepth)

#bin maximum depths 
tot_bin<-allSamples %>% mutate(new_bin = cut(maximumDepth, breaks=c(10, 20, 30,40,50,60,70,80)))

#make a frequency table of bins
freqTab<-table(tot_bin$new_bin, tot_bin$sample)

df<-as.data.frame(freqTab)

ggplot(df, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="dodge", stat="identity")+
  theme_classic()+
  scale_x_discrete(name ="Maximum Depth (M)", 
                   labels=c("(0-10]","(10-20]","(20-30]", "(30-40]", "(40-50]", "(50-60]", "(60-70]", "(70-80]"))+
  scale_fill_discrete(name = "Sampling Interval", labels = c("1 Sec", "300 Sec"))+
  ylab("Frequency")+
  theme_article()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size=14,face="bold"),legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14))

