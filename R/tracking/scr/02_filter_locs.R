#-------------------------------------------------------------------------------------
# filter_locs    Filter locations
#-------------------------------------------------------------------------------------
# This script processes animal tracking data following a common approach between
# different species.
#
# Main steps are:
# - Selection of tracks given a defined criteria
# - Filter location data: Near-duplicate positions, filter, angle and point on land
# - Different processing for PTT/GPS (eg. time gaps, LC class)



#---------------------------------------------------------------
# 2. Select data
#---------------------------------------------------------------

# import metadata
# basic data for further processing: codeName, instrumentType, organismID
organism_meta <- read.csv(paste0(output_dir, "/tracking/metadata/metadataL0.csv"))

# import landmask
world <- ne_countries(scale = "medium", returnclass = "sf")


#---------------------------------------------------------------
# 3. Filter data for each organismID from metadata
#---------------------------------------------------------------

for (i in 1:nrow(organism_meta)){

  print(paste("Processing tag", i))
  
  # get basic information
  indiv <- organism_meta$organismID[i]
  sname <- organism_meta$scientificName[i]
  cname <- organism_meta$codeName[i]
  instr <- organism_meta$instrumentType[i]

  # set filter parameters according to taxonomic group
  filt_vmax <- seaturtle_filt_vmax
  
  # import L0 location data
  input_data <- paste0(output_dir, "/tracking/locdata/L0_loc/")
  infile <- sprintf(paste0(input_data, "%s_L0_loc.csv"), indiv)
  data <- read.csv(infile)
  data$time <- parse_date_time(data$time, "Ymd HMS") # parse time
  
  
  # Trim tracks into segments 
  # Tracks with data gaps in excess of a given time are broken up for separate modeling
  data$tripID <- timedif.segment(data$time, thrs = trip_time_gap)
  data$tripID <- paste(indiv, str_pad(data$trip, 3, pad = "0"), sep="_")

  ### Select trips according to multiple criteria
  
  # summarize data per trip
  trips <- summarizeTrips(data = data, id = "organismID", trip = "tripID", date ="time", lon = "longitude", lat = "latitude")

  # filter trips
  trips <- filter(trips,
                  duration_h >= sel_min_dur,
                  n_loc >= sel_min_loc,
                  distance_km >= sel_min_dist, 
                  !id %in% sel_exclude)

  ### Process data per trip
  
  # select trips
  selected_trips <- unique(trips$trip)
  if(length(selected_trips) == 0) next
  
  # create empty list to append data
  trip_data <- list()
  
  for(j in 1:length(selected_trips)){
    
    ## subset data per trip
    sdata <- dplyr::filter(data, tripID == selected_trips[j])
    
    ## Process Argos data
    if(instr == "Argos"){
      
      # Set params
      filt_step_time <- argos_filt_step_time
      filt_step_dist <- argos_filt_step_dist
      
      # Remove near-duplicate positions
      sdata <- filter_dup(sdata, step.time = filt_step_time, step.dist = filt_step_dist)
      
      # Filter out Z location classess
      sdata <- filter(sdata, argosLC != "Z")
      
      # Filter positions by speed and angle for PTT
      sdata$argosfilter <- sdafilter(lat = sdata$latitude,
                                    lon = sdata$longitude,
                                    dtime = sdata$time,
                                    lc = sdata$argosLC,
                                    vmax = filt_vmax, # in m/s
                                    ang = filt_ang, # No spikes are removed if ang=-1
                                    distlim = filt_distlim)
      
      # select position that are not removed by the filter
      sdata <- dplyr::filter(sdata, argosfilter == "not")
    }

    
    ## append data
    trip_data[[j]] <- sdata
  }
  
  # combine trip data
  data <- rbindlist(trip_data, fill=TRUE)
  if(nrow(data) == 0) next
  
  # write level 1 results
  output_data <- paste0(output_dir, "/tracking/locdata/L1_loc/")
  if (!dir.exists(output_data)) dir.create(output_data, recursive = TRUE)
  outfile <- sprintf(paste0(output_data, "%s_L1_loc.csv"), indiv)
  write.csv(data, outfile, row.names = F)
  
  ### Plot map 
  ### define lon/lat bounds
  # range from data
  xl <- extendrange(data$longitude)
  yl <- extendrange(data$latitude)
  # get centroid
  zoom_to <- c(mean(xl), mean(yl))  # center of the range
  # define zoom level
  lon_span <- xl[2]-xl[1]
  lat_span <- yl[2]-yl[1]
  zoom_lon <- floor(log2(360/lon_span))
  zoom_lat <- floor(log2(180/lat_span))
  zoom_level <- min(zoom_lon, zoom_lat)
  # define span
  lon_span <- 360 / 2^zoom_level
  lat_span <- 180 / 2^zoom_level
  # define boundaries
  lon_bounds <- c(zoom_to[1] - lon_span / 2, zoom_to[1] + lon_span / 2)
  lat_bounds <- c(zoom_to[2] - lat_span / 2, zoom_to[2] + lat_span / 2)
  
  
  # Plot
  p <- ggplot() +
    # add land
    geom_sf(fill="grey80", colour="grey80", linewidth=0.3, data=world)+
    # add tracks
    geom_path(data = data, 
              aes(x = longitude, y = latitude, group=organismID, color=organismID),
              size=0.8) +
    # spatial bounds
    # define spatial extent
    #p <- p +
    coord_sf(xlim = lon_bounds, ylim = lat_bounds, expand = T, ndiscr = 1000) +
    # labs
    labs(title = sname, subtitle=indiv, x ="", y = "") +
    # theme
    theme_article() +
    theme(legend.position = "none")
  
  # Export map
  outfile <- sprintf(paste0(output_data, "%s_L1_loc.png"), indiv)
  ggsave(outfile, p, width=20, height=10, units="cm", dpi=150, bg="white")

}



#---------------------------
# Update metadata
#---------------------------

# Import L0 metadata
metadata <- read.csv(paste0(output_dir, "/tracking/metadata/metadataL0.csv"))

# import all location data
data_dir <- paste0(output_dir, "/tracking/locdata/")
txt_files <- list.files(data_dir, full.names = TRUE, recursive = TRUE, pattern = "L1_loc.csv")
data <- rbindlist(lapply(txt_files, fread), fill=TRUE)

# Find deployment from metadata that were not processed
sel <- which(metadata$organismID %in% data$organismID)
selected <- metadata[sel,]
discarded <- metadata[-sel,]

# Export metadata
write.csv(selected, paste0(output_dir, "/tracking/metadata/metadataL1.csv"), row.names=F)
write.csv(discarded, paste0(output_dir, "/tracking/metadata/discardedL1.csv"), row.names=F)


print("Filtering ready")
