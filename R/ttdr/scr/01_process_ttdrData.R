#Script to process low resolution dives Part 1. (based on D. March turtle dataset). Jessica Harvey-Carroll (carrolljessh@gmail.com/ jessica.carroll@bioenv.gu.se). 2023#

# TTDR processing low resolution consists of the following steps:
# 1. Data standardization
# 2. Zero offset correction. 
# 3. Identfication of ascent and descent phases
# 4. Predict times turtle crossed given threshold based on trajectory
# 5. Relabel dive phases 
# 6.  Prepare data for manipulation 
# 7. Remove dives if NA's fall within dive intervals
# 8. Remove dives if gap in transmission 
# 9. Split dives if turtle can cross 5m threshold based on ascent speed
#10. Relabel new processed data
#11. if transmission stops during the last dive then remove
#12. Split 'hangouts' if possible for turtle the reach surface and return to next reading

#set paths
input_data <- paste0(input_dir, "/tracking/")
output_data <- paste0(output_dir, "/TTDR")
metadata <- read.csv(paste0(input_dir, "/tracking/TODB_2023-12-09_diveAnalysis.csv"))


for (i in 1:nrow(metadata)){
################################################################################
#------------------------------------------------------------------------------#
# Step 1. Process meta-data #
#------------------------------------------------------------------------------#
################################################################################
organismID <- metadata$ptt[i]

print(paste("Processing", organismID))
input_ttdrData <- paste0(input_data,organismID,"/",organismID, "-Series.csv" )

data <- read.csv(input_ttdrData)
# TTDR data, WCT 5 minute bins so:
tfreq <- 5 * 60  # time interval from TTDR data, in seconds

depl <- parse_date_time(metadata$deploymentDateTime[i], c("dmY"), tz="UTC")

#wc2ttdr function Processes Wildlife Computers Time-Series csv files
ttdr <- wc2ttdr(data, date_deploy = depl, tfreq = "5 min")  

#find and save NA depth readings
depthNA <- ttdr %>%
  filter(is.na(depth))

ttdr <- ttdr %>%
  filter(!is.na(depth))

ttdr<- ttdr %>% filter(row_number() <= last(which(!is.na(depth))))
################################################################################
#------------------------------------------------------------------------------#
# Step 2. Zero offset correction #
#------------------------------------------------------------------------------#
################################################################################
# remove duplicates
# priority to registers with [depth & temp] > [depth] > [temp]
if (any(duplicated(ttdr$date))){
  ttdr$idx <- 1:nrow(ttdr)
  dups <- ttdr %>% group_by(date) %>% filter(n()>1) %>% filter(is.na(depth))
  ttdr <- ttdr[-dups$idx,]
  ttdr <- dplyr::select(ttdr, -idx)
}

##Create a TDR class##
#TDR is the simplest class of objects used to represent Time-Depth recorders data in diveMove.
tdrdata <- createTDR(time = ttdr$date, depth = ttdr$depth,
                     dtime = 300,  # sampling interval (in seconds)
                     file = "ttdr.csv")  # path to the file

#change this number for desired dive threshold
customDive_threshold=8.00


## Calibrate with ZOC using filter method.
#The method consists of recursively smoothing and filtering the input time series 
#using moving quantiles.It uses a sequence of window widths and quantiles, and starts
#by filtering the time series using the first window width and quantile in the specified sequences
dcalib <-calibrateDepth(tdrdata,
                        wet.thr = 3610,  # (seconds) At-sea phases shorter than this threshold will be considered as trivial wet.Delete periods of wet activity that are too short to be compared with other wet periods.
                        dive.thr = customDive_threshold,    # (meters) threshold depth below which an underwater phase should be considered a dive.
                        zoc.method ="filter",  # see Luque and Fried (2011)
                        k = c(12, 240),  # (60 and 1200 minutes) Vector of moving window width integers to be applied sequentially.
                        probs = c(0.5, 0.05),   # Vector of quantiles to extract at each step indicated by k (so it must be as long as k)
                        depth.bounds = c(-5, 15),  # minimum and maximum depth to bound the search for the surface.
                        na.rm = TRUE)

## Incorporate adjusted depth and offset. Now use this as raw data.
ttdr <- cbind(ttdr, depth_adj = dcalib@tdr@depth, depth_offset = ttdr$depth - dcalib@tdr@depth)

#calculate error ranges for device
## Calculate depth error
d_error <- depth_error(depth = ttdr$depth_adj, drange = ttdr$drange)
ttdr$depthUpperError <- round(d_error$upper.error, 2)
ttdr$depthLowerError <- round(d_error$lower.error, 2)

# Calculate temperature error
ttdr$temperatureError <- temp_error(ttdr$trange)

# Temperature regional range test (Mediterranean Argo)
# 1: good data; 4: bad data
ttdr$temperatureQC1 <- trange_test(ttdr$temperature, ttdr$temperatureError, tmin=10, tmax=40)

depthTempError_file <- paste0(output_data, "depthError", organismID, ".csv")
write.csv(ttdr, depthTempError_file)

################################################################################
#------------------------------------------------------------------------------#
# Step 3. Identify phases #
#------------------------------------------------------------------------------#
################################################################################
#specify columns for code
ttdr<-subset(ttdr, select=c(date,depth,drange, temperature,trange,depth_adj,depth_offset))

#identify dive and surface phases using a given threshold ##
ttdr_thresh<-ttdr %>%
  mutate(phase = case_when(depth_adj >= customDive_threshold ~ 'dive',
                           TRUE ~ 'surface'))


#first reading of dive below threshold is labelled as decent. Last reading below threshold is marked as ascent.
ttdr_threshLab<-ttdr_thresh %>%
  mutate(thresh = case_when(lag(phase) == "surface" & phase == "dive" ~ 'descent',
                            lag(phase) == "dive" & phase == "surface" ~ 'ascent',
                            lag(phase) == "dive" & phase == "dive" ~ 'dive',
                            lag(phase) == "surface" & phase == "surface" ~ 'surface',
                            TRUE ~ ''))

#identify descents to count dive numbers#
desc<-ttdr_threshLab %>% filter(thresh == "descent")
desc$diveno <- 1:nrow(desc) 

#add numbered descents to labelled ttdr_thresh (ttdr_threshLab)
ttdr_numDesc<-merge(ttdr_threshLab, desc, by = 'date', all.x= TRUE)

#remove empty duplicate rows from y(labelled descent) dataframe from merging dataframes
ttdr_numDesc <- subset(ttdr_numDesc, select=-c(depth.y,drange.y,temperature.y,trange.y,depth_adj.y,depth_offset.y,phase.y,thresh.y))

#repeate dive numbers
ttdr_numDescRep<-na.locf(ttdr_numDesc,na.rm = FALSE)

#remove repeates for surface intervals
ttdr_clean <- within(ttdr_numDescRep, diveno[thresh.x == 'surface'] <- NA)

################################################################################
#------------------------------------------------------------------------------#
# Step 4. Predict times turtle crossed given threshold based on trajectory #
#------------------------------------------------------------------------------#
################################################################################
#calculate the point turtles pass descent/ ascent thresholds for each dive. Add in artificial times.if predicted time same as next time then omit

##descent/acent.func # calculate the point turtles pass descent thresholds for each dive. Add in artificial times. If predicted time same as next time then omit

ttdr_desc<-descent.func(ttdr_clean)
ttdr_descAsc<-ascent.func(ttdr_desc)

#descent and ascent times now in columns, need to combine added threshold times into single date-time column.
#merge columnds together
data_ascDesc<-ttdr_descAsc

orig<-subset(data_ascDesc, select=c(date, depth.x,drange.x,temperature.x,trange.x,depth_adj.x,depth_offset.x,phase.x, thresh.x,diveno))

data_descTime<-subset(data_ascDesc, select=c(desc.thresh.tim,desc.thresh.dep,descenttim))#descent time (descenttim), depth and seconds (for reference)

data_ascTime<-subset(data_ascDesc, select=c(asc.thresh.tim, asc.thresh.dep, ascenttim)) #ascent time (ascenttim), depth and seconds (for reference)

# rename third column as date and 2nd column as depth for descent and ascent dataframes. 6th column in original dataset is the adjusted depth from ZOC and renamed as depth
colnames(data_ascTime)[colnames(data_ascTime) == "ascenttim"] ="date"
colnames(data_ascTime)[colnames(data_ascTime) == "asc.thresh.dep"] ="depth"

colnames(data_descTime)[colnames(data_descTime) == "desc.thresh.dep"] ="depth"
colnames(data_descTime)[colnames(data_descTime) == "descenttim"] ="date"

colnames(orig)[colnames(orig) == "depth_adj.x"] ="depth"

#add for merge
data_ascTime$phase.x<- paste0("dive")
data_ascTime$thresh.x<-paste0("ascent")

data_descTime$phase.x<- paste0("dive")
data_descTime$thresh.x<-paste0("descent")

#merge original and ascent
dataMerge<-merge(x = orig, y = data_ascTime, by=c("date", "depth", "phase.x", "thresh.x"), all = TRUE)

dataMerge<-dataMerge %>% drop_na(date)

#merge new data and descent
dataMerge2<-merge(x = dataMerge, y = data_descTime, by=c("date", "depth","phase.x", "thresh.x"), all = TRUE)

colnames(dataMerge2)[colnames(dataMerge2) == "depth"] ="depth_adj"

dataMerge2<-dataMerge2 %>% drop_na(date)

dataMergeLocf<-na.locf(dataMerge2,na.rm = FALSE)


################################################################################
#------------------------------------------------------------------------------#
# Step 5. Relabel dive phases #
#------------------------------------------------------------------------------#
################################################################################
#relabel phases after adding in ascent and descent predictions#
ascDescPhase<-dataMergeLocf %>%
  mutate(new_thresh = case_when(lag(phase.x) == "surface" & phase.x == "dive" ~ 'dive',
                                lag(phase.x) == "dive" & phase.x == "surface" ~ 'surface',
                                lag(phase.x) == "dive" & phase.x == "dive" ~ 'dive',
                                lag(phase.x) == "surface" & phase.x == "surface" ~ 'surface',
                                TRUE ~ ''))
ascDescPhase2<-ascDescPhase %>%
  mutate(new_thresh_2 = case_when(lag(phase.x) == "surface" & phase.x == "dive" ~ 'descent',
                                  lag(phase.x) == "dive" & phase.x == "surface" ~ 'ascent',
                                  lag(phase.x) == "dive" & phase.x == "dive" ~ 'dive',
                                  lag(phase.x) == "surface" & phase.x == "surface" ~ 'surface',
                                  TRUE ~ ''))

#relabel dive numbers again
ascDescFilt<-ascDescPhase2 %>% filter(new_thresh_2 == "descent")
ascDescFilt$diveno_2 <- 1:nrow(ascDescFilt) 

#merge dive numbers with data
ascDescNums<-merge(ascDescPhase2, ascDescFilt, by = 'date', all.x= TRUE)

#remove empty duplicate columns for second dataframe. Keep diveno_2
ascDescNums<-subset(ascDescNums, select=-c(depth_adj.y, phase.x.y,thresh.x.y,depth.x.y, drange.x.y,temperature.x.y,trange.x.y, depth_offset.x.y,diveno.y,asc.thresh.tim.y,desc.thresh.tim.y, new_thresh.y,new_thresh_2.y))

#carry foward dive numbers
ascDescNumsLocf<-na.locf(ascDescNums,na.rm = FALSE)

#make sure surface is NA for dive numbers
ascDescFin <- within(ascDescNumsLocf, diveno_2[new_thresh.x == 'surface'] <- NA)

#new dataframe using date, adjusted depth, new_threshold and dive number
ascDescFin_short<-subset(ascDescFin, select = c(date,depth_adj.x,new_thresh.x,diveno_2))

colnames(ascDescFin_short)[colnames(ascDescFin_short) == "depth_adj.x"] ="depth"
colnames(ascDescFin_short)[colnames(ascDescFin_short) == "new_thresh.x"] ="phase"
colnames(ascDescFin_short)[colnames(ascDescFin_short) == "diveno_2"] ="diveno"

#remove any duplicated rows#
ascDescFin_short<-ascDescFin_short[!duplicated(ascDescFin_short), ]
ascDescFin_short

#write CSV of preprocessed data!!
#write.csv(ascDescFin_short, "C:/Users/user/Desktop/tort/finalrun_7mthresh/preprocessed_dives_151935.csv")
ascDescFin_short_file <- paste0(output_data, "preprocessed_dives_", organismID, ".csv")
write.csv(ascDescFin_short, ascDescFin_short_file)

################################################################################
#------------------------------------------------------------------------------#
# Step 6. Prepare data for predictions #
#------------------------------------------------------------------------------#
################################################################################
ascDescFin_short<- read.csv(ascDescFin_short_file)

##grep to add 00's to midnight
ascDescFin_short$date[grep("[0-9]{4}-[0-9]{2}-[0-9]{2}$",ascDescFin_short$date)] <- paste(
  ascDescFin_short$date[grep("[0-9]{4}-[0-9]{2}-[0-9]{2}$",ascDescFin_short$date)],"00:00:00")

#make sure datetime is posixct
ascDescFin_short$date <- as.POSIXct(ascDescFin_short$date,format="%Y-%m-%d %H:%M:%S",tz="CET")

#new df for dives
dataDive<-ascDescFin_short %>%
  filter(phase == "dive")

#new df for surface
dataSurf<-ascDescFin_short %>%
  filter(phase == "surface")

#split dives by dive number
data_sep<-split(dataDive, dataDive$diveno)

#find total number of dives
lastDive<-tail((names(data_sep)), n=1)
divNum<-as.numeric(lastDive[1])

#character string of dive names
names<-sprintf("dive%d", 1:divNum)
split_names <- as.character(names)

#make each dive a seperate dataframe now
for (i in 1:length(data_sep)) {     
  assign(split_names[i], data_sep[[i]])
}


################################################################################
#------------------------------------------------------------------------------#
# Step 7. Remove dives if NA's fall within dive intervals
#------------------------------------------------------------------------------#
################################################################################
# newint function to get start and fin time of each dive. Apply to find start and fin times of all dives
for(i in  1:length(names)){
  data =  eval(parse(text = paste0( names[i])))
  assign(paste(names[i],"interval",sep = ""), t<-newint(data))
  print(paste(names[i],"interval",sep = ""))
}

#merge all start and fin times into single dataframe
interval_list = mget(ls(pattern = "interval"))
interval_listbind= do.call(what = rbind, args = interval_list)
interval_listbind<-interval_listbind[order(interval_listbind$diveno),]

#make sure all in same timezones
interval_listbind$start<-force_tz(interval_listbind$start, "UTC")
interval_listbind$fin<-force_tz(interval_listbind$fin, "UTC")

#depthNA is times with missing depths from first step
depthNA$date<-force_tz(depthNA$date, "UTC")

#find dives where depth NA's present
df1<-  interval_listbind
df2<- depthNA

setDT(df1); setDT(df2)
# initialise new column with "N"
df1[, dateoccurs := "N"]
#update join#
hasNA<-  df1[df2, dateoccurs := "Y", on = .(start <= date, fin >= date)][]
hasNA<-as.data.frame(hasNA)
torem<-  hasNA %>% filter (dateoccurs== "Y")

#string of dive id's to remove
rem<-torem$diveno

tostay<-  hasNA %>% filter (dateoccurs== "N")

#string of dives to stay
stay<-tostay$diveno

#character string of dives to stay (with no depth NA's within)
names1<-sprintf("dive%d", stay)

################################################################################
#------------------------------------------------------------------------------#
# Step 8. Remove dives if gap in transmission #
#------------------------------------------------------------------------------#
################################################################################
dfgap<-ttdr
dfgap$gap <- c(NA, with(dfgap, date[-1] - date[-nrow(dfgap)]))

# now, how often was the gap more than some amount of time?
gap_threshold <- 6 #minutes #300 secs
dfgap$over_thresh <-as.integer(as.logical( dfgap$gap > gap_threshold)) #new column highlights if gap was over previously given value
dfgap

#find entries where gap over given time occurs
gapFlag<-which(dfgap$over_thresh ==1)
gapFlagLast<-c(gapFlag-1,gapFlag) #last recording before gap

#times where last reading before gap
gapFlaga<-as.data.frame(dfgap[gapFlagLast,1])

names(gapFlaga)[names(gapFlaga) == "dfgap[gapFlagLast, 1]"] <- "date"

## flag for dive to remove because large gap of failed recordings
df1<-  interval_listbind
df2<- gapFlaga
setDT(df1); setDT(df2)
# initialise new column with "N"
df1[, dateoccurs := "N"]
#update join#
hasGap<-  df1[df2, dateoccurs := "Y", on = .(start <= date, fin >= date)][]
hasGap<-as.data.frame(hasGap)
torem2<-  hasGap %>% filter (dateoccurs== "Y")
rem2<-torem2$diveno

tostay2<-  hasGap %>% filter (dateoccurs== "N")
stay2<-tostay2$diveno

names2<-sprintf("dive%d", stay2)

###dives to carry forward, combining dives to keep from names1 (NA depth) and names2 (failed recordings)
common_values<-intersect(names1,names2)
names<-common_values

################################################################################
#------------------------------------------------------------------------------#
# Step 9. Split dives if turtle can cross 5m threshold based on ascent speed#
#------------------------------------------------------------------------------#
################################################################################

# function to split dive if peaks detected and possible for turtle to surface based on trajectory
#screen all dives, apply when necessary
for(i in  1:length(names)){
  data =  eval(parse(text = paste0(names[i])))
 # print(data)
  m<-split.dive(data)
  w<-is.null(m)
  if (w== FALSE){
    t<-assign(paste(names[i],"_split",sep = ""),m)
  }
  else (t<-assign(paste(names[i],"_notsplit",sep = ""),data))
}

#gather and bind split/ not split dives
split_list = mget(ls(pattern = "_split"))
notsplit_list= mget(ls(pattern = "_notsplit"))

bygsplit_data = do.call(what = rbind, args = split_list)
bygnotsplit_data = do.call(what = rbind, args = notsplit_list)

#merge all dives
splitBind <- rbind(bygsplit_data, bygnotsplit_data)

#merge together surface and dives
diveSurfBind <- rbind(splitBind, dataSurf)
#order by date
diveSurfBind_ord<-diveSurfBind[order(diveSurfBind$date),]

diveSurfBind_ord<-diveSurfBind_ord[!duplicated(diveSurfBind_ord), ]

diveSurfBind_ord_file <- paste0(output_data, "split_dive_", organismID, ".csv")
write.csv(diveSurfBind_ord, diveSurfBind_ord_file)


################################################################################
#------------------------------------------------------------------------------#
# Step 10. Relabel new processed data#
#------------------------------------------------------------------------------#
################################################################################

diveSurfBind_ord<- read.csv(diveSurfBind_ord_file)

##grep to add 00's to midnight
diveSurfBind_ord$date[grep("[0-9]{4}-[0-9]{2}-[0-9]{2}$",diveSurfBind_ord$date)] <- paste(
  diveSurfBind_ord$date[grep("[0-9]{4}-[0-9]{2}-[0-9]{2}$",diveSurfBind_ord$date)],"00:00:00")

diveSurfBind_ord$date <- as.POSIXct(diveSurfBind_ord$date,format="%Y-%m-%d %H:%M:%S",tz="CET") ## Your dates need to be in 

diveSurfBind_ord<- diveSurfBind_ord %>% drop_na(date)
#colnames(diveSurfBind_ord)

diveSurfBind_ord<-subset(diveSurfBind_ord, select = c(date,depth))

#leave datetime and depth
diveSurfBind_ord

#again rounded to customDive_threshold for threshold

# function to label surface/ dive phases
data_labsFinal<-phases(diveSurfBind_ord)

# functions to calculate depth difference (m), time difference (s) and speed difference (m/s) from previous entry

data_labsFinal<-depth.diff(data_labsFinal)
data_labsFinal<-time.diff(data_labsFinal)  
data_labsFinal<-speed(data_labsFinal)

#identify and label changes in direction, negative = decrease in depth
data_labsDirec<-data_labsFinal %>%
  mutate(direction = case_when(diff > 0 ~ 'positive',
                               diff < 0 ~ 'negative',
                               diff == 0 ~ 'no change',
                               TRUE ~ ''))

#now need to loop through each dive, get summary stats and put in new df
#split dives by number
data_sep<-split(data_labsDirec, data_labsDirec$diveno)

#find total number of dives
lastDive<-tail((names(data_sep)), n=1)
divNum<-as.numeric(lastDive[1])

#print names of all dives
names<-sprintf("dive%d", 1:divNum) 

split_names <- as.character(names)

#make each dive a seperate dataframe
for (i in 1:length(data_sep)) {
  assign(split_names[i], data_sep[[i]])
}

################################################################################
#------------------------------------------------------------------------------#
# Step 11. if transmission stops during the last dive then remove#
#------------------------------------------------------------------------------#
################################################################################

#if transmission stops during the last dive then remove 'incomplete' dive
last.dive.name<-tail(names, n=1)
last.dive<-eval(parse(text = paste0(last.dive.name)))
if(tail(last.dive$depth,n=1)>=customDive_threshold){
 names<-names[1:length(names)-1]
}else{
  names<-names
}
#run dive.phaselab function to label phases of dive
#apply to all dives
for(i in  1:length(names)){
  data =  eval(parse(text = paste0( names[i])))
  assign(paste(names[i],"phase",sep = "_"), t<-dive.phaselab(data))
  print(paste(names[i],"phase",sep = "_"))
}

################################################################################
#------------------------------------------------------------------------------#
# Step 12. Split 'hangouts' if possible for turtle the reach surface and return to next reading#
#------------------------------------------------------------------------------#
################################################################################
#split.hangouts function to split 'hangouts' if possible for turtle to reach surface based on trajectory. Note- Hangouts indicate ascent, followed by time spent stationary at depth and then descent

final<- paste(names, "phase", sep = "_")

#screen all dives, apply when necessary
for(i in  1:length(final)){
  print(final[i])
  data =  eval(parse(text = paste0(final[i])))
 #print(final[i])
  m<-split.hangouts(data)
  w<-is.null(m)
  if (w== FALSE){
    t<-assign(paste(final[i],"_splithangs",sep = ""),m)#, m<-split.dive(data))
    print(paste(final[i],"_splithangs",sep = ""))
  }
}

if (length(mget(ls(pattern = "_splithangs")))>0){

#gather and bind split/ not split dives
splithang_list = mget(ls(pattern = "_splithangs"))
#notsplithang_list= mget(ls(pattern = "_notsplithangs"))
bygsplithang_data = do.call(what = rbind, args = splithang_list)

#merge split hang dives with previous split dives 
#make sure time and depth only, remove bogus columns

newhangs<-subset(bygsplithang_data, select = c(date,depth))
dives<-subset(splitBind, select = c(date,depth))#splitBind[ -c(3:5) ]

#combine total dives
splitBind2 <- rbind(newhangs, dives)

#now add surface values
#remove bogus columns
dataSurf2<-subset(dataSurf, select = c(date,depth)) #dataSurf[ -c(1,4,5) ]

#merge together for final dataset for processing
diveSurfBind2 <- rbind(splitBind2, dataSurf2)

diveSurfBind2.1<-diveSurfBind2[order(diveSurfBind2$date),]
final_dives<-diveSurfBind2.1[!duplicated(diveSurfBind2.1), ]
}else{
 dives<-subset(splitBind, select = c(date,depth))
  splitBind2 <- dives
  
  #now add surface values
  #remove bogus columns
 dataSurf2<-subset(dataSurf, select = c(date,depth)) #dataSurf[ -c(1,4,5) ]
  
  #merge together for final dataset for processing
  diveSurfBind2 <- rbind(splitBind2, dataSurf2)
  
  diveSurfBind2.1<-diveSurfBind2[order(diveSurfBind2$date),]
 final_dives<-diveSurfBind2.1[!duplicated(diveSurfBind2.1), ]
}
  
final_dives_file <- paste0(output_data, "final_forprocess_", organismID, ".csv")
write.csv(final_dives, final_dives_file)

#clean environment from previous loop iteration
rm(list=setdiff(ls(), c("output_data","input_data","metadata", "i")))
}
################################## Finish loops#################################
beep()
