#----------------------------------------------------------------------------
# process_dives
#----------------------------------------------------------------------------




#---------------------------------------------------------------
# 1. Set data repository
#---------------------------------------------------------------
input_data <- paste0(output_dir, "/tracking/", sp_code, "/TTDR")
output_data <- paste0(output_dir, "/tracking/", sp_code, "/dives")
if (!dir.exists(output_data)) dir.create(output_data, recursive = TRUE)



#---------------------------------------------------------------
# 2. Import data
#---------------------------------------------------------------

# Import all ttdr. readTrack
ttdr_files <- list.files(input_data, full.names=TRUE, pattern = "_L1_ttdr.csv")



#---------------------------------------------------------------
# 3. Process individual files
#---------------------------------------------------------------


cl <- makeCluster(6)
registerDoParallel(cl)

foreach(i=1:length(ttdr_files), .packages=c("dplyr", "stringr", "lubridate", "diveMove")) %dopar% {
  
  # import ttdr
  data <- readTrack(ttdr_files[i])
  id <- data$id[i]
  
  # import dcalib
  dcalib <- readRDS(paste0(input_data, "/", id, "_dcalib.rds"))
  
  # dive summary
  dive_sum <- diveSummary(dcalib) %>%
    filter(begdesc >= min(data$date), endasc <= max(data$date))
  
  # calculate vertical speeds (m/s)
  dive_sum$desc.speed <- dive_sum$descdist/dive_sum$desctim
  dive_sum$asc.speed <- dive_sum$ascdist/dive_sum$asctim
  
  dive_sum$depth_upper_error <- NA
  dive_sum$depth_lower_error <- NA
  for (j in 1:nrow(dive_sum)){
    
    # subset TDR data during dive
    sdata <- data %>%
      filter(date >= dive_sum$begdesc[j], date <= dive_sum$endasc[j])
    
    # find max depth and get errors
    sel <- which(sdata$depth_adj == max(sdata$depth_adj, na.rm=T))
    dive_sum$depth_upper_error[j] <- sdata$depth_upper_error[sel]
    dive_sum$depth_lower_error[j] <- sdata$depth_lower_error[sel]
  }
  
  
  # Add latitute and longitude from TTDR dataset
  dive_sum <- dive_sum %>% inner_join(dplyr::select(data, date, lon, lat), by = c("begdesc" = "date"))
  
  # select variables
  df <- dplyr::select(dive_sum, dive.id, lon, lat, begdesc, endasc, divetim,
                      pdd, pdd_qc, maxdep, depth_upper_error, depth_lower_error, dtype, xy_error)
  
  # add id
  df <- cbind(id, df)
  
  # Export data
  outfile <- paste0(output_data, "/", id, "_dive.csv")
  write.csv(df, outfile, row.names = FALSE)
}

stopCluster(cl)




