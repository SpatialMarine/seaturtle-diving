#------------------------------------------------------------------------------
# setup.R         Project set up
#------------------------------------------------------------------------------

# set parameters
source("config.R")

# Load required packages
pacman::p_load("data.table", "tidyr", "dplyr", "lubridate", "stringr", "readxl", "tools",# data manipulation
               "ggplot2", "egg", "pals", "ggpubr", # plots
               "diveMove", "foieGras", # diving and movement
               "rnaturalearthdata", "rnaturalearth",  # base data
               "foreach", "doParallel",  # parallel computing
               "raster", "sf", "animalsensor", "move", "argosfilter", "geosphere")  # spatial

# Create data paths
input_dir <- paste(main_dir, "input", sep="/")
if (!dir.exists(input_dir)) dir.create(input_dir, recursive = TRUE)

output_dir <- paste(main_dir, "output", sep="/")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
