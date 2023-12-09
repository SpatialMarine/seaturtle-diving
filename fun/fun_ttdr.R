#------------------------------------------------------------------------------------
# fun_ttdr.R    Suite of function for processing TTDR data
#------------------------------------------------------------------------------------
# aboveSST     Temperature above SST test
# daynight          Return if location with time stamp is during day/nigh
# depth_error.R     Function to calculate depth error from Time-Series data on SPLASH tags
# diveSummary       Classify dive types and QC for postdive
# getSST            Generate SST product
# temp_error     Function to calculate depth error from Time-Series data on SPLASH tags
# trange_test     QC test using global temperature range
# wc2ttdr  Processes Wildlife Computers Time-Series csv files




#--------------------------------------------------------------------------------
# aboveSST     Temperature above SST test
#--------------------------------------------------------------------------------
aboveSST <- function(temp, temp.er, sst, temp.thr = 4){
  # In order to be more consevative, we use the lower band of the temperature data (temperature - temperature error)
  
  tdif <- (temp - temp.er) - sst
  ifelse(tdif >= temp.thr, 4, 1)  # 4: bad data; 1: good data
}
#--------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# daynight          Return if location with time stamp is during day/nigh
#-------------------------------------------------------------------------------
daynight <- function(lon, lat, time){
  # given a vectors of lon, lat and time, returns "day" or "night"
  
  library(maptools)
  
  # internal function to classify day/night
  daynight <- function (time, sunrise, sunset){
    ifelse(time >= sunrise & time <= sunset, "day", "night")
  }
  
  ### Incorporate sun and lunar metrics
  ### Calculate: sunrise and sunset times; moon phase
  sunrise <- sunriset(crds=cbind(lon, lat), dateTime=time, direction=c("sunrise"), POSIXct.out=TRUE)$time
  sunset <- sunriset(crds=cbind(lon, lat), dateTime=time, direction=c("sunset"), POSIXct.out=TRUE)$time
  
  ## Derive day/night
  daynight <- mapply(daynight, time=time, sunrise=sunrise, sunset=sunset)
  return(daynight)
}
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------
# depth_error.R     Function to calculate depth error from Time-Series data on SPLASH tags
#-------------------------------------------------------------------------------------------
depth_error <- function(depth, drange, percent=0.01, res=0.5){
  # Arguments:
  # depth           Depth reading (meters, positive value)
  # drange          Resolution of the reading
  # percent         Percentage of the reading
  # res             Sensor resolution
  #
  # Value:
  # Upper and lower depth error (in meters, positive value)
  #
  # Information about the sensor:
  # range: -40 to +1000m
  # resolution: 0.5 m
  # accuracy: +-1% of the reading. Also: 1% of the reading +/-2 depth sensor resolutions
  #
  # Decription:
  # This function only applies for Time-Series data from Wildlife Computers. 
  # Dynamic range of the dive readings during the summary period determine
  # the Time-Series Data resolution (Kevin Lay, Wildlife Computers).
  #
  # Usage:
  # depth <- 60
  # drange <- 6.25
  # e <- depth_error(depth=depth, drange)
  
  # Calculate upper and lower resolution bands
  upper.band <- depth + drange
  lower.band <- depth - drange
  
  # Calculate sensor accuracy for upper and lower bands
  upper.accuracy <- upper.band * percent + 2 * res
  lower.accuracy <- lower.band * percent + 2 * res
  
  # Calculate total depth error (resolution + accuracy)
  upper.error <- upper.band + upper.accuracy
  lower.error <- lower.band - lower.accuracy
  # if(lower.error < 0) lower.error <- 0
  if (any(lower.error < 0)) lower.error[lower.error < 0] <- 0
  
  return(data.frame(upper.error=upper.error, lower.error=lower.error))
}
#-------------------------------------------------------------------------------------------



#------------------------------------------------------------------
# diveSummary         Classify dive types and QC for postdive
#------------------------------------------------------------------
# Description: Classify dives into one of three simple dive shapes (Square, V, U)
diveSummary <- function(dcalib){
  
  library(akima)
  library(dplyr)
  
  ### Data frame with summary statistics for each dive
  diveDF <- diveStats(dcalib, depth.deriv=TRUE)
  
  # add dive id
  diveDF$dive.id <- c(1:nrow(diveDF))
  
  # correction of dive duration (divetime) considering the difference of 300s between depth 
  #at begdesc and depth at the same time in data
  diveDF$begdesc <- diveDF$begdesc - 300
  diveDF$divetim <- diveDF$divetim + 300
  diveDF$endasc <- diveDF$begdesc + diveDF$divetim  # time of beggining of descent + dive time (s)
  
  # Filtered dives by time
  diveDF <- filter(diveDF, divetim > 300)  # filter out dives <300s duration
  data <- data.frame(date = dcalib@tdr@time, depth = dcalib@tdr@depth)
  
  ####################
  # Dive Type Analysis
  ####################
  diveDF$botdepth <- diveDF$maxdep*0.8  # bottom depth (criteria by WC)
  
  ## function to calculate bottom time duration data
  
  # for dive i
  # select beg.time (begdesc) of dive i => t1
  # select end.time (endasc) of dive i => t2
  # select TDR data from t1 to t2
  diveDF$botttim <- NA
  diveDF$dtype <- NULL
  for (i  in 1:(nrow(diveDF))) {
    
    # get time period of a dive
    t1 <- diveDF$begdesc[i] # begining of descent of dive i
    t2 <- diveDF$endasc[i]  # end ascent of dive i
    
    # select TDR data within a dive 
    
    fdata <- filter(data, date >= t1 & date <= t2)
    d <- fdata$depth#depths
    t <- fdata$date# time (in seconds)
    
    #xout <- seq(t1,t2,60)# time sequence to predict
    xout <- seq(t1,t2,by="15 sec")
    
    t <- as.numeric(t-t1)
    xout <- as.numeric(xout-t1)
    
    pred <- aspline(t, d, xout=xout, n=2)  # prediction
    bottom <- which(pred$y>diveDF$botdepth[i])# get depths below bottom depth
    botstart <- xout[bottom[1]]  # first time where depth is > than bottom depth
    botend <- xout[bottom[length(bottom)]] # last time where depth is > than bottom depth
    diveDF$botttim[i] <- botend-botstart # duration of the bottom time
    
    # dive type
    if (diveDF$botttim[i] > diveDF$divetim[i]*0.5) diveDF$dtype[i] <- "S"
    if (diveDF$botttim[i] <= diveDF$divetim[i]*0.2) diveDF$dtype[i] <- "V"
    if (diveDF$botttim[i] > diveDF$divetim[i]*0.2 & diveDF$botttim[i] <= diveDF$divetim[i]*0.5) diveDF$dtype[i] <- "U"
  }
  
  
  ####################
  # Postdive duration
  ####################
  
  # Calculate postdive duration
  diveDF$pdd <- NULL
  for (i  in 1:(nrow(diveDF)-1)){
    endasc <- diveDF$endasc[i]  # end ascent of dive i
    begdesc <- diveDF$begdesc[i+1]  # begining of descent of dive i+1
    diveDF$pdd[i] <- difftime(begdesc, endasc, units="secs")
  }
  
  # Quality control of Postdive Duration
  # for dive i
  # select end.time of dive i => t1
  # select start.time of dive i+1 => t2
  # select TDR data from t1 to t2
  # is NA data in TDR subset?
  #    yes: bad PDD
  #    no: good PDD
  # return a data frame (dive id, ppd_qc)
  
  diveDF$pdd_qc <- NULL
  for (i  in 1:(nrow(diveDF)-1)){
    
    # get time period of postdive duration
    t1 <- diveDF$endasc[i]  # end ascent of dive i
    t2 <- diveDF$begdesc[i+1]  # begining of descent of dive i+1
    
    # select TDR data within such period
    fdata <- filter(data, date >= t1 & date <= t2)
    nacheck <- any(is.na(fdata$depth))
    
    # define quality of postive data
    # (1) good, (4) bad
    ifelse(nacheck == TRUE, diveDF$pdd_qc[i] <- 4, diveDF$pdd_qc[i] <- 1) 
  }
  
  
  # Reorder with dive.id as first column
  diveDF <- diveDF %>%
    dplyr::select(dive.id, dtype, pdd, pdd_qc, everything())
  
  return(diveDF) 
}
#------------------------------------------------------------------




#--------------------------------------------------------------------------------
# getSST     Generate SST product
#--------------------------------------------------------------------------------
getSST <- function(data){
  # This function fixed several problems found on the initial test.
  # 1) I found high temperature values for some points. This is mainly because I initially consider
  #    all dives and postdives. I have sorted this out by:
  #    a) filtering dives >10 minutes long (ie. at least 2 registers per dive).
  #    b) filtering out dives <5 m deep (use lower error estimate). shallow dives could be affected)
  # 2) There were low temperature values for two reasons:
  #    a) Some postdives present without depth data. I filtered those out
  #    b) First register of postdive very similar to previous underwater register. I selected the second register per postdive. 
  
  # ttdr.df       Data.frame generated in advance
  
  
  ## First, select dives with more than 10 minutes of duration, with at least one point at >5m
  dive10 <- data %>%
    group_by(dive.id) %>%
    summarize(n=n(),
              maxdep=max(depthLowerError)) %>%
    filter(n>2, maxdep > 5)
  
  ## Indetify second time stamp of a postdive
  sst <- data %>%
    group_by(postdive.id) %>%
    arrange(time) %>%
    slice(2) %>%  # take the second register from the postdive
    ungroup
  
  ## Filter data
  sst <- filter(sst, postdive.id > 3,  # remove first 3 dives
                postdive.id %in% dive10$dive.id,
                temperatureQC1 == 1,  # select temperature data within range
                !is.na(depth)) # select data with depth information
  
  ## Select fields
  sst <- dplyr::select(sst, organismID, time, depthAdjusted, depthUpperError,
                       depthLowerError, temperature, temperatureError, temperatureQC1)
  return(sst)
}
#--------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------
# temp_error     Function to calculate depth error from Time-Series data on SPLASH tags
#-------------------------------------------------------------------------------------------
temp_error <- function(trange, accuracy = 0.1){
  # Arguments:
  # trange          Resolution of the reading
  # percent         Accuracy of the sensor
  #
  # Value:
  # Temperature error (in degrees)
  #
  # Information about the sensor:
  # range: -40 to 60 degrees celsius
  # resolution: 0.05 degrees
  # accuracy: +- 0.1 degrees
  #
  # Decription:
  # This function only applies for Time-Series data from Wildlife Computers. 
  # Dynamic range of the readings during the summary period determine
  # the Time-Series Data resolution (Kevin Lay, Wildlife Computers).
  #
  # Usage:
  # trange <- 0.85
  # e <- temp_error(trange)
  
  # Calculate total depth error (resolution + accuracy)
  e <- trange + accuracy
  return(e)
}
#-------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# trange_test     QC test using global temperature range
#--------------------------------------------------------------------------------
trange_test <- function(temp, temp.error=0, tmin=10, tmax=40){
  # Arguments
  # t         temperature
  # terror    temperature error. Default is 0 (no error)
  # tmin      minimum temperature range
  # tmax      maximum temperature range
  #
  # Values
  # QC flags: 0, 1, 2, 3, 4, 5
  # Returns NA when temperature data is not available
  #
  # Note
  # temperature range (tmin and tmax) have been calculated for Western Mediterranean...
  # Check what trange variable from Wildlife Computers means
  #
  # Refs:
  # - Argo regional values for Mediterranean Sea (10-40): http://www.argodatamgt.org/content/download/341/2650/file/argo-quality-control-manual-V2.7.pdf
  # - https://www.ukargo.net/data/quality_control/
  #
  # Usage
  # temp <- 18.8
  # temp.error <- 0.95
  # trange_test(temp, temp.error)
  ifelse(temp + temp.error >= tmin & temp - temp.error <= tmax, 1, 4)  # 1: good data; 4: bad data
}
#--------------------------------------------------------------------------------


#---------------------------------------------------------------------
# wc2ttdr  Processes Wildlife Computers Time-Series csv files
#---------------------------------------------------------------------
wc2ttdr <- function(data, locale = "English", date_deploy=NULL, tfreq = "5 min"){
  
  require(dplyr)
  require(lubridate)
  
  # Convert to POSIXct
  data$time <- paste(data$Day, data$Time)
  data$time <- parse_date_time(data$time, c("dmY HMS", "Ymd HMS"), locale=locale, tz="UTC")
  
  # Select data collected from the date of deployment
  if (!is.null(date_deploy)) data <- filter(data, time >= date_deploy)

  ### Rename columns
  names(data)[names(data)=="Ptt"] <- "organismID"
  names(data)[names(data)=="Depth"] <- "depth"
  names(data)[names(data)=="DRange"] <- "depthRange"
  names(data)[names(data)=="Temperature"] <- "temperature"
  names(data)[names(data)=="TRange"] <- "temperatureRange"
  
  ## create a time series for the whole period with NAs
  id <- data$organismID[1]  # Get animal id
  min.time <- min(data$time)  # min time stamp
  max.time <- max(data$time)  # max time stamp
  ts.seq <- seq(from = min.time, to = max.time, by = tfreq)  # create complete TS at regular period of time
  ts.seq <- data.frame(time = ts.seq)  # convert to data.frame
  data <- merge(ts.seq, data, by="time", all.x = TRUE)  # merge data with complete sequence
  data$organismID[is.na(data$organismID)] <- id
  
  # Reorder column names
  data <- dplyr::select(data, organismID, time, depth, depthRange, temperature, temperatureRange)
  return(data)
}
#---------------------------------------------------------------------



