#--------------------------------------------------------------------------------
# process_ttdr.R
#--------------------------------------------------------------------------------
# Process pressure and temperature data from WC tags configured using time seies
#
# TTDR processing includes the following steps:
# 1. Data standardization
# 2. Zero offset correction. + Identfication of ascent and descent phases
# 3. Estimation of depth and temperature errors
# 4. Temperature QC (regional test)
# 5. Interpolate locations from SSM
#
# Then, other scripts can be followed:
# 2. SST product and QC to detected insolation events.
# 3. Diving analysis
# 4. Derive MLD product
# 5. Compare SST and MLD with other products
# Custom functions are found in "scr/fun_ttdr"


#---------------------------------------------------------------
# 1. Set data repository
#---------------------------------------------------------------
location_data <-  paste0(output_dir, "/tracking/locdata/L1_loc")
output_data <- paste0(output_dir, "/ttdr")
if (!dir.exists(output_data)) dir.create(output_data, recursive = TRUE)




#---------------------------------------------------------------
# 2. Process metadata
#---------------------------------------------------------------

# import metadata
metadata <- read.csv(paste0(output_dir, "/tracking/metadata/metadataL1.csv"))





#---------------------------------------------------------------
# 3. Process TTDR
#---------------------------------------------------------------

#cl <- makeCluster(10)
#registerDoParallel(cl)

## Process delayed mode (location data)
for (i in 1:nrow(metadata)){
#foreach(i=1:nrow(db), .packages=c("dplyr", "stringr", "lubridate", "diveMove", "move", "foieGras")) %dopar% {
  
  print(paste("Processing tag", i, "of", nrow(db)))
  
  ## get tag information
  id <- metadata$organismID[i]
  
  # Import data
  loc_file <- list.files(input_dir, recursive=TRUE, full.names=TRUE, pattern = sprintf("%s-Series.csv", id))
  data <- read.csv(loc_file)
  
  #-------------------------------
  # Step 1. Data standardization
  #-------------------------------
  
  # Standardize data
  # Need to keep the file for diveMove
  ttdr <- wc2ttdr(data, locale = locale, date_deploy = metadata$deploymentDateTime[i], tfreq = "5 min")  
  ttdr_file <- paste0(output_data, "/", id, "_L0_ttdr.csv")
  write.csv(ttdr, ttdr_file, row.names = FALSE)
  
  
  #-------------------------------
  # Step 2. Zero offset correction
  #-------------------------------
  
  # remove duplicates
  # priority to registers with [depth & temp] > [depth] > [temp]
  if (any(duplicated(ttdr$time))){
    ttdr$idx <- 1:nrow(ttdr)
    dups <- ttdr %>% group_by(time) %>% filter(n()>1) %>% filter(is.na(depth))
    ttdr <- ttdr[-dups$idx,]
    ttdr <- dplyr::select(ttdr, -idx)
  }
  
  ##Create a TDR class##
  #TDR is the simplest class of objects used to represent Time-Depth recorders data in diveMove.
  tdrdata <- createTDR(time = ttdr$time, depth = ttdr$depth,
                       dtime = 300,  # sampling interval (in seconds)
                       file = ttdr_file)  # path to the file
  
  ## Calibrate with ZOC using filter method.
  #The method consists of recursively smoothing and filtering the input time series 
  #using moving quantiles.It uses a sequence of window widths and quantiles, and starts
  #by filtering the time series using the first window width and quantile in the specified
  #sequences 
  dcalib <-calibrateDepth(tdrdata,
                          wet.thr = 3610,  # (seconds) At-sea phases shorter than this threshold will be considered as trivial wet.Delete periods of wet activity that are too short to be compared with other wet periods.
                          dive.thr = 3,    # (meters) threshold depth below which an underwater phase should be considered a dive.
                          zoc.method ="filter",  # see Luque and Fried (2011)
                          k = c(12, 240),  # (60 and 1200 minutes) Vector of moving window width integers to be applied sequentially.
                          probs = c(0.5, 0.05),   # Vector of quantiles to extract at each step indicated by k (so it must be as long as k)
                          depth.bounds = c(-5, 15),  # minimum and maximum depth to bound the search for the surface.
                          na.rm = TRUE)
  
  # save dcalib object for further analysis
  saveRDS(dcalib, paste0(output_data, "/", id, "_dcalib.rds"))
  
  ## Incorporate adjusted depth and offset
  ttdr <- cbind(ttdr, depthAdjusted = dcalib@tdr@depth, depthOffset = ttdr$depth - dcalib@tdr@depth)
  
  
  #-------------------------------
  # Step 3. QC for depth
  #-------------------------------
  
  # Calculate depth error
  d_error <- depth_error(depth = ttdr$depthAdjusted, drange = ttdr$depthRange)
  ttdr$depthUpperError <- round(d_error$upper.error, 2)
  ttdr$depthLowerError <- round(d_error$lower.error, 2)
  
  
  #-------------------------------
  # Step 4. QC for temperature
  #-------------------------------
  
  # Calculate temperature error
  ttdr$temperatureError <- temp_error(ttdr$temperatureRange)
  
  # Temperature regional range test (Mediterranean Argo)
  # 1: good data; 4: bad data
  ttdr$temperatureQC1 <- trange_test(ttdr$temperature, ttdr$temperatureError, tmin=10, tmax=40)
  
  
  #-------------------------------
  # Step 5. Derive diving metrics
  #-------------------------------

  ## dive.activity: data frame with details about all dive and postdive periods found"dive.id",
  #"dive.activity", and "postdive.id"
  #L:dry, W:wet, U:underwater, D:diving, Z:brief wet 
  dive_act <- getDAct(dcalib)
  
  #dive.phases: This identifies each reading with a particular dive phase. Thus, each reading belongs to one
  #of descent (D), descent/bottom (DB), bottom (B), bottom/ascent (BA), and ascent (A) phases.
  dive_phases <- getDPhaseLab (dcalib)
  
  ## Merge diving data with Time Series
  ttdr <- cbind(ttdr, dive_act, dive_phases)
  
  
  #-------------------------------
  # Step 6. QC for SST
  #-------------------------------
  
  #-------------------------------
  # 5.1. Generate SST product
  #-------------------------------
  ### Generate SST product
  ### Select first records for each postdive.
  ### We asume that for these records there is no insolation effect.
  ### This function remove first 3 dives and select data
  sst <- getSST(ttdr)  
  
  ### Interpolate SST to TimeSeries
  ttdr$sst <- approx(x = sst$time, y = sst$temperature, xout = ttdr$time, method="linear", rule=2)$y
  ttdr$sst <- round(ttdr$sst, digits=1)
  ttdr$sst_qc <- 8  # interpolated data
  ttdr$sst_qc[ttdr$time %in% sst$time] <- 1  # good data, no interpolated
  
  ### SST product
  #L1_sst <- select(data, ptt, date, sst, sst_qc)
  
  
  #-------------------------------
  # Temperature above SST test
  #-------------------------------
  ### This test identifies temperature records that are above the estimated SST, given a temperature threshold
  ### We define a temperature threshold of 2 degrees
  ### We do not select points on the surface because there may be no depth available when having temperature data
  ### 4: bad data; 1: good data
  ttdr$aboveSST <- aboveSST(temp = ttdr$temperature, temp.er = ttdr$temperatureError,
                                sst = ttdr$sst, temp.thr = 1)
  

  
  #-------------------------------
  # Step 7. Add locations
  #-------------------------------
  loc_file <- paste0(location_data, "/", id, "_L1_loc.csv")
  loc_filt <- read.csv(loc_file)
  loc_filt$time <- parse_date_time(loc_filt$time, "Ymd HMS") # parse time
  
  # Calculate 95% confidence interval of location error
  # SSM provides standard errors in km
  # Calculate SE average between lon and lat
  # Estimate 95% confidence interval
  #loc_filt$xy_error <- 1.96*1000*((loc_filt$x.se +  loc_filt$y.se)/2)  # return in m
  

  
  ### Interpolate -----------------------------------------
  
  # No extrapolation is performed. So check that time stamps in series are within time domain of loc
  ttdr <- dplyr::filter(ttdr, time >= min(loc_filt$time) & time <= max(loc_filt$time))
  
  # convert location data to move class
  loc_filt_move <- move(x=loc_filt$longitude, y=loc_filt$latitude, time=loc_filt$time, data=loc_filt,
                        proj=CRS("+proj=longlat +ellps=WGS84"), animal=loc_filt$organismID)
  
  # linear interpolation
  lint <- interpolateTime(loc_filt_move, time = ttdr$time, spaceMethod = "greatcircle")
  
  # incorporate coordinates to TTDR data.frame
  ttdr$longitude <- lint@coords[,1]
  ttdr$latitude <- lint@coords[,2]
  
  # Interpolate horizontal error to TTDR data
  #ttdr$xy_error <- akima::aspline(x=(as.numeric(loc_filt$date)), y=(loc_filt$xy_error), xout=(as.numeric(ttdr$date)))$y
  
  
  #-------------------------------
  # Step 8. Add day/night
  #-------------------------------
  ttdr$daynight <- daynight(lon = ttdr$longitude, lat = ttdr$latitude, time = ttdr$time)
  
  
  #-------------------------------
  # Export data
  #-------------------------------
  ttdr_file <- paste0(output_data, "/", id, "_L1_ttdr.csv")
  write.csv(ttdr, ttdr_file, row.names = FALSE)
}
