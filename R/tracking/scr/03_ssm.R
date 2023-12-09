#-------------------------------------------------------------------------------------
# 03_regularize_ssm     Interpolate tracks into regular time steps using foieGras
#-------------------------------------------------------------------------------------
#
# Main steps are:
# - Regularize tracks

#devtools::install_github("ianjonsen/foieGras")  # github version fixes problem with GPS
library("foieGras")


#---------------------------------------------------------------
# Prepare cluster
#---------------------------------------------------------------
cl <- makeCluster(cores)
registerDoParallel(cl)


#---------------------------------------------------------------
# 2. Select data
#---------------------------------------------------------------

# import metadata
# basic data for further processing: codeName, instrumentType, organismID
metadata <- read.csv(paste0(output_dir, "/tracking/metadata/metadataL1.csv"))

# summarize metadata per organismID
organism_meta <- metadata %>%
  group_by(organismID, codeName, instrumentType) %>%
  summarize()




#---------------------------------------------------------------
# 3. Filter data for each organismID from metadata
#---------------------------------------------------------------

#for (i in 1:nrow(organism_meta)){
foreach(i=1:nrow(organism_meta), .packages=c("dplyr", "ggplot2", "foieGras", "stringr", "lubridate", "animalsensor")) %dopar% {  

  print(i)
  
  # get basic information
  indiv <- organism_meta$organismID[i]
  cname <- organism_meta$codeName[i]
  
  # import L1 location data
  input_data <- paste0(output_dir, "/tracking/locdata/", cname, "/L1_loc/")
  infile <- sprintf(paste0(input_data, "%s_L1_loc.csv"), indiv)
  data <- read.csv(infile)
  data$time <- parse_date_time(data$time, "Ymd HMS") # parse time

  # summarize data per trip
  trips <- summarizeTrips(data)

  # filter trips
  trips <- filter(trips,
                  duration_h >= sel_min_dur,
                  n_loc >= sel_min_loc,
                  distance_km >= sel_min_dist, 
                  !id %in% sel_exclude)
  
  # subset data
  # filter by id and selected trips
  data <- filter(data, tripID %in% trips$trip)
  
  ###### State-Space Model
  
  # convert to foieGras format
  indata <- data %>%
    rename(id = tripID,
           date = time,
           lc = argosLC,
           lon = longitude,
           lat = latitude,
           smaj = argosSemiMajor,
           smin = argosSemiMinor,
           eor = argosOrientation) %>%
    dplyr::select(id, date, lc, lon, lat, smaj, smin, eor)
  
  # filter location class data with NA values
  # very few cases, but creates an error in fit_ssm
  indata <- dplyr::filter(indata, !is.na(lc))
  
  # fit SSM
  # we turn sdafilter off because we previously filtered data
  # we run the model with multiple trips at once
  fit <- fit_ssm(indata, model = "crw", time.step = reg_time_step,
                 control = ssm_control(verbose = 0), spdf = FALSE,
                 map = list(psi = factor(NA)))

  # fit Time-varying move persistence
  # Jonsen et al. 2018. https://doi.org/10.1002/ecy.2566
  fmp <- fit_mpm(fit, what = "predicted", model = "jmpm", control = mpm_control(verbose = 0))
  
  #plot(fmp, pages = 1, ncol = 3, pal = "Cividis", rev = TRUE)
  #m <- fmap(fit, fmp, what = "predicted", pal = "Cividis")
  
  # get fitted behavioural state
  # Time-varying move persistence
  datafmp <- data.frame(grab(fmp, what = "fitted", as_sf = FALSE))
  
  # get fitted locations
  data <- data.frame(grab(fit, what = "predicted", as_sf = FALSE))
  
  # combine and arrange data
  data <- data %>%
    # join datasets
    left_join(datafmp, by=c("id", "date")) %>%
    # rename and prepare data for further steps
    rename(tripID = id, time = date, longitude = lon, latitude = lat) %>%
    arrange(time) %>%
    mutate(organismID = indiv) %>%
    dplyr::select(organismID, everything())
  
  # write level 2 results
  output_data <- paste0(output_dir, "/tracking/locdata/", cname, "/L2_loc/")
  if (!dir.exists(output_data)) dir.create(output_data, recursive = TRUE)
  outfile <- sprintf(paste0(output_data, "%s_L2_loc.csv"), indiv)
  write.csv(data, outfile, row.names = F)
  
  # export convergence status
  convergence <- data.frame(organismID = indiv, tripID = fit$id, converged_fit = fit$converged, converged_fmp = fmp$converged)
  outfile <- sprintf(paste0(output_data, "%s_L2_convergence.csv"), indiv)
  write.csv(convergence, outfile, row.names = FALSE)
}


stopCluster(cl) # Stop cluster
print("Regularization ready")

#---------------------------------------------------------------
# 4. Summarize processed data
#---------------------------------------------------------------

summarizeId <- function(data, id = "organismID", date ="time", lon = "longitude", lat = "latitude"){

  df <- data %>%
    # rename varaibles
    dplyr::rename(id = `id`, date = `date`, lon = `lon`, lat = `lat`) %>%
    # order by date
    dplyr::arrange(date) %>%
    # group by id and trip
    dplyr::group_by(id) %>%
    dplyr::summarize(date_deploy = first(date),
                     lon_deploy = first(lon),
                     lat_deploy = first(lat),
                     date_last = last(date),
                     time_interval_h = median(as.numeric(difftime(tail(date, -1), head(date, -1), units="hours"))),
                     distance_km = sum(geosphere::distGeo(p1 = cbind(lon, lat)), na.rm=TRUE)/1000,  # segment distance
                     n_loc = n()) %>%  # get first and last observations
    dplyr::mutate(duration_d = round(difftime(date_last, date_deploy, units="days"))) %>%  # calculate duration of the track
    # rename to original variables
    dplyr::rename(`id` = id)

  return(df)
}

# Import L0 metadata
metadata <- read.csv(paste0(output_dir, "/tracking/metadata/metadataL1.csv"))

# import all location data
data_dir <- paste0(output_dir, "/tracking/locdata/")
txt_files <- list.files(data_dir, full.names = TRUE, recursive = TRUE, pattern = "L2_loc.csv")
data <- rbindlist(lapply(txt_files, fread), fill=TRUE)

# summarize data per trip
tripsId <- summarizeId(data)

# combine data
selected <- metadata %>%
  left_join(tripsId, by = c("organismID" = "id"))

# Export metadata
write.csv(selected, paste0(output_dir, "/tracking/metadata/metadataL2.csv"), row.names=F)
