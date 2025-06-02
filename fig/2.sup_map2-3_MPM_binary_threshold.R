

# Process locs data to obtain g value threshold for id.


# load data:
depths_bind <-read.csv("C:/Users/J. Menéndez Blázquez/SML_Dropbox/SML Dropbox/gitdata/seaturtle-diving/output/jess/plots/dive_summariesLabStates060224_FIN1304.csv")

# prepare unique turtle data states
newdf <- depths_bind %>% filter(depths_bind$maximumDepth > 10)
rm(depths_bind)

colnames(newdf)[colnames(newdf) == "latitude"] ="Latitude"
colnames(newdf)[colnames(newdf) == "longitude"] ="Longitude"

ids <- unique(newdf$organismID)

# empty dataframe to store new threshold value
df <- data.frame(id = as.character(),
                 g_threshold = as.integer())


for (i in seq_along(ids)) {
  # extract ID of interest
  id <- ids[i]
  # read L2 loc files
  locs <- read.csv(paste0(output_dir,"/tracking/locdata/L2_loc/",id,"_L2_loc.csv"))
  
  g_threshold <- median(locs$g)

  # results threshold
  r <- data.frame(id = id,
                  g_threshold = g_threshold)
  
  # add to df results
  df <- rbind(df,r)
}

write.csv(df, paste0(output_dir,"/tracking/mpm/L2_loc_mpm_threshold.csv") )
