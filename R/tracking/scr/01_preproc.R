#-----------------------------------------------------------------------------------------
# 01_preproc.R        Pre-process loggerhead tracking data
#-----------------------------------------------------------------------------------------
# This script pre-processes tracking data. The main goal is to standardize
# the tracking dataset


#---------------------
# Set parameters
#---------------------

## define dataset metadata
#datasetID <- "icss"
#ownerInstitutionAbbrev <- "ICSS"

# set input data folder
deployment_metadata <- sprintf("%s/tracking/TODB_2023-12-09.xlsx", input_dir)

# set output metadata folder
output_metadata <- sprintf("%s/tracking/metadata/", output_dir)
if (!dir.exists(output_metadata)) dir.create(output_metadata, recursive = TRUE)


#---------------------
# Deployment data
#---------------------

# import metadata
metadata <- read_excel(deployment_metadata, sheet = "metadata")

# prepare deployment metadata
metadata <- metadata %>%
  # rename variables and change to character
  mutate(codeName = getSpeciesCode(scientificName),
         organismID = ptt) %>%
  # recode variables
  dplyr::mutate(organismSex = recode(organismSex, Male = "male", Female = "female")) %>%
  # filter animals to process
  dplyr::filter(wcTimeSeries == "y", codeName == "CARCAR")

# Export metadata
write.csv(metadata, paste0(output_dir, "/tracking/metadata/metadataL0.csv"), row.names=F)

#---------------------------------------------------------------
# 3. Process tracking data
#---------------------------------------------------------------

## Process delayed mode (location data)
for (i in 1:nrow(metadata)){

  print(paste("Processing tag", i, "of", nrow(metadata)))
  
  ## get tag information
  id <- metadata$ptt[i]
  codeName <- metadata$codeName[i]

  # If data from Wildlife Computers -------------------------
  # One folder per individual

    # Import Argos data
    loc_file <- list.files(input_dir, recursive=TRUE, full.names=TRUE, pattern = sprintf("%s-Locations.csv", id))
    data <- read.csv(loc_file)
    # Standardize data  
    dataL0 <- animalsensor::wcLoc2L0(data, locale = locale, date_deploy = metadata$deploymentDateTime[i]) 
  
    # set output data folder
    output_data <- sprintf("%s//tracking/locdata/L0_loc/", output_dir)
    if (!dir.exists(output_data)) dir.create(output_data, recursive = TRUE)
    
    # export data
    outfile <- sprintf("%s/%s_L0_loc.csv", output_data, id)
    write.csv(dataL0, outfile, row.names = F)
}

print("Pre-processing ready")

