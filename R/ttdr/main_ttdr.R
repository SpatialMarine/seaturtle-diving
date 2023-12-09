#-----------------------------------------------------------------------
# main_tracking.R
#-----------------------------------------------------------------------

source("setup.R")
source("fun/fun_ttdr.R")

#---------------------------------------------------------------
# 1. Set parameters for processing
#---------------------------------------------------------------

# TTDR data
tfreq <- 5 * 60  # time interval from TTDR data, in seconds



#---------------------------------------------------------------
# Processing TTDR data
#---------------------------------------------------------------


### TDR data
source("R/ttdr/scr/01_process_ttdr.R")  # process TDR



source("R/ttdr/scr/02_process_dives.R")  # convert to dive


