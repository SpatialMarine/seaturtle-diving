#-----------------------------------------------------------------------
# main_tracking.R
#-----------------------------------------------------------------------

source("setup.R")


#---------------------------------------------------------------
# 1. Set parameters for processing
#---------------------------------------------------------------

### Trip definition
trip_time_gap <- 7 * 24  # (used if trip_type == time) Tracks with data gaps in excess of [seg_time_gap] hours were broken up for separate modeling

### Track/trip selection
sel_min_loc <- 20  # minimum number of locations per trip
sel_min_dur <- 10 * 24 # minimum duration of trip, in hours
sel_min_dist <- 15 # minimum distance of tracks, in km
sel_exclude <- NULL # custom selection of tags based on exploration of data

### Track filtering

# Params for Argos filter
argos_filt_step_time <- 2/60  # time difference to consider duplicated positions, in hours
argos_filt_step_dist <- 0/1000  # spatial distance to consider duplicated poisitions, in km
filt_ang <- c(15, 25) # value of the angle using in sdafilter, no spikes are removed if ang=-1
filt_distlim <- c(2500, 5000) # value of the limite distance using in sdafilter, no spikes are removed if ang=-1

# Params for all
filt_land <- FALSE  # remove locations on land

# Params by species
seaturtle_filt_vmax <- 2  # value of the maximum of velocity, in m/s


### Track regularization
reg_time_step <- 2  # time step to interpolate positions, in hours




#---------------------------------------------------------------
# Processing location data
#---------------------------------------------------------------

#-----------------------------------------------
# Step 1. Pre-process data and standardize data (L0)
#-----------------------------------------------
# Import all tracking datasets
# Transforms different data sources into a common format
# convert tracking data into L0
# create a metadata file with all deployments
# Output data can be found: /output/tracking/locdata (and metadata)
# - /locdata: tracking data files are stored per organismID
# - /metadata: metadata records represent deployments
source("R/tracking/scr/01_preproc.R")


#-----------------------------------------------
# Step 2. Filter data (L1)
#-----------------------------------------------
source("R/tracking/scr/02_filter_locs.R")


#-----------------------------------------------
# Step 3. Regularize data
#-----------------------------------------------
# SSM for Argos tags
cores <- 10
source("analysis/tracking/scr/03_ssm.R")
