#Clean code to process low resolution dives (based on D. March turtle dataset). Jessica Harvey-Carroll (carrolljessh@gmail.com/ jessica.carroll@bioenv.gu.se). 2023#

rm(list=ls())
require(dplyr)
require(lubridate)
require(diveMove)
library(zoo)
library(changepoint)
library(akima)
library(data.table)
library(beepr)
library(pracma)
library(changepoint)
library(tidyr)

setwd("C:/Users/user/Desktop/tort")

# TTDR data, WCT 5 minute bins so:
tfreq <- 5 * 60  # time interval from TTDR data, in seconds
data<- read.csv("C:/Users/user/Desktop/tort/input/tracking/CAR/wc/DONE/181762/181762-Series.csv")

#meta data

meta<-read.csv("TODB_2022-07-23.csv")

organismID<- meta %>%filter (ptt== "181762")

depl <- parse_date_time(organismID$deploymentDateTime, c("dmY"), tz="UTC")

###for turtle 200043 input is different
# depl <- parse_date_time(format(as.POSIXct(organismID$deploymentDateTime,format='%d/%m/%Y %H:%M'),format='%d/%m/%Y'), c("dmY"), tz="UTC")

#-------------------------------------------------------------------------------
# wc2ttdr  Processes Wildlife Computers Time-Series csv files
#-------------------------------------------------------------------------------
wc2ttdr <- function(data, locale = "English", date_deploy=NULL, tfreq = "5 min"){
  
  
  # Convert to POSIXct
  data$date <- paste(data$Day, data$Time)
  data$date <- parse_date_time(data$date, "dmY HMS", tz="UTC")
  
  # Select data collected from the date of deployment
  if (!is.null(date_deploy)) data <- filter(data, date >= date_deploy)
  
  ### Rename columns
  names(data)[names(data)=="Ptt"] <- "id"
  names(data)[names(data)=="Depth"] <- "depth"
  names(data)[names(data)=="DRange"] <- "drange"
  names(data)[names(data)=="Temperature"] <- "temperature"
  names(data)[names(data)=="TRange"] <- "trange"
  
  ## create a time series for the whole period with NAs
  # id <- data$id[1]  # Get animal id
  min.time <- min(data$date)  # min time stamp
  max.time <- max(data$date)  # max time stamp
  ts.seq <- seq(from = min.time, to = max.time, by = tfreq)  # create complete TS at regular period of time
  ts.seq <- data.frame(date = ts.seq)  # convert to data.frame
  data <- merge(ts.seq, data, by="date", all.x = TRUE)  # merge data with complete sequence
  #data$id[is.na(data$id)] <- id
  
  # Reorder column names
  data <- dplyr::select(data, date, depth, drange, temperature, trange)
  return(data)
}
#-------------------------------------------------------------------------------

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
customDive_threshold=7


## Calibrate with ZOC using filter method.
#The method consists of recursively smoothing and filtering the input time series 
#using moving quantiles.It uses a sequence of window widths and quantiles, and starts
#by filtering the time series using the first window width and quantile in the specified
#sequences 
#NOTE organismID 200043 = 7m, 200045= 6m and 151935= 6m dive threshold.
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


#-------------------------------------------------------------------------------
#descent.func # calculate the point turtles pass descent thresholds for each dive. Add in artificial times. If predicted time same as next time then omit
#-------------------------------------------------------------------------------
descent.func<-function(dataframe){
  
  rwz = which(dataframe[,9] == "descent")
  desc.thresh.tim = rep(NA, nrow(dataframe))
  desc.thresh.dep = rep(NA, nrow(dataframe))
  
  for(i in 1:length(rwz)){
    t2<- dataframe[rwz[i],1] 
    t1 <- dataframe[rwz[i]-1,1]  
    
    # select TDR data within a dive 
    fdata <- filter(dataframe, date >= t1 & date <= t2)
    d <- fdata$depth_adj.x#depths
    t <- fdata$date #time (in seconds)
    
    xout <- seq(t1,t2,by="10 sec")
    t <- as.numeric(t-t1)
    xout <- as.numeric(xout-t1)
    
    pred <- aspline(t, d, xout=xout, n=2)  # prediction
    
    pt <- which(pred$y >= customDive_threshold)# get depths below dive threshold
    
    td1= xout[pt[1]]  # first time where depth is > than threshold
    if (td1 == 300){
      ptmin=pt-1
      td1= xout[ptmin[1]] 
    }else{
      td1=td1
    }
    
    desc.thresh.tim[rwz[i]] = td1
  }
  
  for(j in 1:length(rwz)){
    t2<- dataframe[rwz[j],1] 
    t1 <- dataframe[rwz[j]-1,1]  
    
    # select TDR data within a dive 
    fdata <- filter(dataframe, date >= t1 & date <= t2)
    d <- fdata$depth_adj.x#depths
    t <- fdata$date #time (in seconds)
    
    xout <- seq(t1,t2,by="10 sec")
    t <- as.numeric(t-t1)
    xout <- as.numeric(xout-t1)
    
    pred <- aspline(t, d, xout=xout, n=2)  # prediction
    
    ##
    pred<-as.data.frame(pred)
    
    pred<-filter(pred,y >= customDive_threshold)
    
    td12<-pred[1,2]
    td12<-round(td12, digits=1)
    
    ####
    desc.thresh.dep[rwz[j]] = td12
  }
  
  
  dataframenew = cbind(dataframe,desc.thresh.tim, desc.thresh.dep)
  
  dataframenew$descenttim <- with(dataframenew, lag(date)+ desc.thresh.tim)
  return(dataframenew)
  
}
#-------------------------------------------------------------------------------
#ascent.func # calculate the point turtles pass ascent thresholds for each dive. Add in artificial times. If predicted time same as next time then omit
#-------------------------------------------------------------------------------
ascent.func<-function(dataframe){
  
  rwz = which(dataframe[,9] == "ascent")
  asc.thresh.tim = rep(NA, nrow(dataframe))
  asc.thresh.dep = rep(NA, nrow(dataframe))
  
  for(i in 1:length(rwz)){
    t2<- dataframe[rwz[i],1] 
    t1 <- dataframe[rwz[i]-1,1]  
    
    # select TDR data within a dive 
    fdata <- filter(dataframe, date >= t1 & date <= t2)
    d <- fdata$depth_adj.x#depths
    t <- fdata$date #time (in seconds)
    
    xout <- seq(t1,t2,by="10 sec")
    t <- as.numeric(t-t1)
    xout <- as.numeric(xout-t1)
    
    pred <- aspline(t, d, xout=xout, n=2)  # prediction
    
    pt <- which(pred$y <= customDive_threshold)# get depths below dive threshold
    
    td1= xout[pt[1]]  # first time where depth is > than threshold
    
    if (td1 == 300){
      ptmin=pt-1
      td1= xout[ptmin[1]] 
    }else{
      td1=td1
    }
    
    asc.thresh.tim[rwz[i]] = td1
    
  }
  
 
    for(j in 1:length(rwz)){
      t2<- dataframe[rwz[j],1] 
      t1 <- dataframe[rwz[j]-1,1]  
      
      # select TDR data within a dive 
      fdata <- filter(dataframe, date >= t1 & date <= t2)
      d <- fdata$depth_adj.x#depths
      t <- fdata$date #time (in seconds)
      
      xout <- seq(t1,t2,by="10 sec")
      t <- as.numeric(t-t1)
      xout <- as.numeric(xout-t1)
      
      pred <- aspline(t, d, xout=xout, n=2)  # prediction
      
      ##
      pred<-as.data.frame(pred)
      
      pred<-filter(pred,y <= customDive_threshold)
      
      td12<-pred[1,2]
      td12<-round(td12, digits=1)
      
      ####
      asc.thresh.dep[rwz[j]] = td12
    }
  dataframenew = cbind(dataframe,asc.thresh.tim, asc.thresh.dep)
  
  dataframenew$ascenttim <- with(dataframenew, lag(date)+ asc.thresh.tim)
  return(dataframenew)
  
}

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
write.csv(ascDescFin_short, "C:/Users/user/Desktop/tort/test/preprocessed_dives_181762.csv")


################################################################################
#------------------------------------------------------------------------------#
# Step 6. Prepare data for manipulation #
#------------------------------------------------------------------------------#
################################################################################
ascDescFin_short<- read.csv("C:/Users/user/Desktop/tort/test/preprocessed_dives_181762.csv")

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
#rm otherwise this doesn't work down the line
rm(interval_list)
rm(interval_listbind)

################################################################################
#------------------------------------------------------------------------------#
# Step 7. Remove dives if NA's fall within dive intervals
#------------------------------------------------------------------------------#
################################################################################
#-------------------------------------------------------------------------------
# newint # function to get start and fin time of each dive
#-------------------------------------------------------------------------------
newint<-function(dataframe){
  ie<-nrow(dataframe)
  a<-dataframe[ie,2]
  b<-dataframe[1,2]
  # d<-dataframe[1,5]
  newdf<- as.data.frame(dataframe[1,5])
  newdf$start<- b
  newdf$fin<- a
  names(newdf)[1] = "diveno"
  return(newdf)
}

#apply funciton to find start and fin times of all dives
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
#-------------------------------------------------------------------------------
#split.dive # function to split dive if peaks detected and possible for turtle to surface based on trajectory
#-------------------------------------------------------------------------------
split.dive<-function(dataframe){

  
  p<-findpeaks(-dataframe$depth, threshold =3)
  j<-as.data.frame(p[,2])
  t<-nrow(j)
  fill=rep(NA, t)
  
  if (t>0){
    for (a in 1:nrow(j)){
      j2<-p[,2][a]  
      #depths& times of peak, n-1, n+1
      depthnm1<--dataframe[j2-1,3]
      depthpoint<-p[,1][a]
      depthnp1<--dataframe[j2+1,3]
      
      timenm1<-dataframe[j2-1,2]
      timepoint<-dataframe[j2,2]
      timenp1<-dataframe[j2+1,2]
      
      
      #timediff
      t1<-as.numeric(timepoint-timenm1)*60
      t2<-as.numeric(timenp1-timepoint)*60
      
      #depthdiff
      d1<-depthpoint-depthnm1
      d2<-depthpoint-depthnp1
      
      #calculate speed
      #first speed
      s1<-d1/t1
      #second speed
      s2<-d2/t2
      
      #could turtle cross threshold based on speed before/after peak
      distance_poss1<-s1*150
      distance_poss2<-s2*150
      
      distneeded<- -customDive_threshold-depthpoint
      
      if(distance_poss1>=distneeded || distance_poss2>=distneeded ){
        fill[a] = j[a,]
      }
      
    }
  }
  
  
  if (t>0 & sum(!is.na(fill))>0){
    
    for (k in 1:length(dataframe)){  
      rwz = fill[!is.na(fill)]
      de1=rep(NA, nrow(dataframe))
      newdv1<- .POSIXct(character())
      de2=rep(NA, nrow(dataframe))
      newdv2<- .POSIXct(character())
      de3=rep(NA, nrow(dataframe))
      newdv3<- .POSIXct(character())
      
      for (i in 1:length(rwz)){
        t1<- dataframe[rwz[i],2]
        t2 <- t1+150 
        t<-as.numeric(t1-t1)
        t<-c(t,150)
        d<-c(dataframe[rwz[i],3], 0)
        xout <- seq(t1,t2,by="2 sec") #was 5 secs but values 5.9 so changed for + precision. Was originally was 10 secs but changed becasue duplicate confusion
        xout <- as.numeric(xout-t1)
        pred <- aspline(t, d, xout=xout, n=2)  # prediction
        pt <- which(pred$y >= 1.00)# get depths below dive threshold
        pt<-rev(pt)
        td1= xout[pt[1]] 
        pt1d<-pt[1]
        
        date<-c(t1+td1)
        depth<- pred$y[pt1d]
        
        newdv1[i] = as.POSIXct(date)
        de1[i] = depth 
        
        #ab.thresh#
        pt2<- which(pred$y >= customDive_threshold)
        pt2<-rev(pt2)
        td1.1= xout[pt2[1]] 
        pt2d<-pt2[1]
        date2<-c(t1+td1.1)
        depth2<- pred$y[pt2d]
        newdv2[i] = as.POSIXct(date2)
        de2[i] = depth2
        
        ###########
        ##descent again, make sure it uses rwz[i], NOT J+1
        tda<- date
        tdb <- dataframe[rwz[i]+1,2] 
        ta<-0 #as.numeric(tdb-tda) #   difftime(tdb,tda, unit="secs"))
        ta<-c(ta,150)
        da<-c(dataframe[rwz[i]+1,3], depth)
        xout2 <- seq(tda, tdb,by="2 sec") #was 10 secs but changed becasue duplicate confusion
        xout2 <- as.numeric(xout2-tda)
        pred2 <- aspline(ta, da, xout=xout2, n=2)  # prediction
        pta <- which(pred2$y >= customDive_threshold)# get depths below dive threshold
        pta<-rev(pta)
        td1a= xout[pta[1]] 
        pt1da<-pta[1]
        
        date3<-c(tdb-td1a)
        depth3<- pred2$y[pt1da]
        newdv3[i] = as.POSIXct(date3)
        de3[i] = depth3
      }
      
      de1<-de1[!is.na(de1)]
      de2<-de2[!is.na(de2)]
      de3<-de3[!is.na(de3)]
      
      new_df <- data.frame(date1 = as.POSIXct(newdv1, format="%Y-%m-%d %H:%M:%S",tz="CET"), 
                           date2=as.POSIXct(newdv2, format="%Y-%m-%d %H:%M:%S",tz="CET"),
                           date3=as.POSIXct(newdv3, format="%Y-%m-%d %H:%M:%S",tz="CET"))
      
      new_df2 <- data.frame(dep1=de1,
                            dep2=de2,
                            dep3=de3)
      
      new1<-data.frame(date = as.POSIXct(newdv1, format="%Y-%m-%d %H:%M:%S",tz="CET"), 
                       dep=de1)
      new2<-data.frame(date=as.POSIXct(newdv2, format="%Y-%m-%d %H:%M:%S",tz="CET"),
                       dep=de2 )
      new3<-data.frame( date=as.POSIXct(newdv3, format="%Y-%m-%d %H:%M:%S",tz="CET"),
                        dep=de3)
      
      df<-rbind(new1,new2,new3)
      names(df)[names(df) == "dep"] <- "depth"
      
      
      #final DF!
      multdiv2 <- merge(dataframe,df, by=c("date", "depth"), all = TRUE)
      # multdiv2<- multdiv2[complete.cases(multdiv2[ , 1:2]),]
      multdiv2.2<-na.locf(multdiv2, na.rm = FALSE)
      return(multdiv2.2)
    }
  }
  
}

#screen all dives, apply when necessary
for(i in  1:length(names)){
  data =  eval(parse(text = paste0(names[i])))
  print(data)
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

write.csv(diveSurfBind_ord, "C:/Users/user/Desktop/tort/test/split_dive_181762.csv")

################################################################################
#------------------------------------------------------------------------------#
# Step 10. Relabel new processed data#
#------------------------------------------------------------------------------#
################################################################################

diveSurfBind_ord<- read.csv("C:/Users/user/Desktop/tort/test/split_dive_181762.csv")

##grep to add 00's to midnight
diveSurfBind_ord$date[grep("[0-9]{4}-[0-9]{2}-[0-9]{2}$",diveSurfBind_ord$date)] <- paste(
  diveSurfBind_ord$date[grep("[0-9]{4}-[0-9]{2}-[0-9]{2}$",diveSurfBind_ord$date)],"00:00:00")

diveSurfBind_ord$date <- as.POSIXct(diveSurfBind_ord$date,format="%Y-%m-%d %H:%M:%S",tz="CET") ## Your dates need to be in 

diveSurfBind_ord<- diveSurfBind_ord %>% drop_na(date)
colnames(diveSurfBind_ord)

diveSurfBind_ord<-subset(diveSurfBind_ord, select = c(date,depth))

#leave datetime and depth
diveSurfBind_ord

#again rounded to 5 for threshold
#-------------------------------------------------------------------------------
# phases # function to label surface/ dive phases
#-------------------------------------------------------------------------------
phases<-function(dataframe){
  
  #assign phase (dive/surface) for threshold value
  dataframe1<-dataframe %>%
    mutate(phase = case_when(round(depth) >= customDive_threshold ~ 'dive',
                             TRUE ~ 'surface'))
  
  #add phase of dives (surface etc.)
  dataframe2<-dataframe1 %>%
    mutate(thresh = case_when(lag(phase) == "surface" & phase == "dive" ~ 'descent',
                              lag(phase) == "dive" & phase == "surface" ~ 'ascent',
                              lag(phase) == "dive" & phase == "dive" ~ 'dive',
                              lag(phase) == "surface" & phase == "surface" ~ 'surface',
                              TRUE ~ ''))
  ##add dive numbers
  desc<-dataframe2 %>% filter(thresh == "descent")
  desc$diveno <- 1:nrow(desc) 
  
  #merge dataframes to label dives
  dataframe3<-merge(dataframe2, desc, by = 'date', all.x= TRUE)
  
  #remove bogus y columns but keep original dataframe and new dive numbers
  dataframe3.1 <- dataframe3[ -c(5: 7) ]
  
  dataframe3.2<-na.locf(dataframe3.1,na.rm = FALSE)
  dataframe4 <- within(dataframe3.2, diveno[thresh.x == 'surface'] <- NA)
  
  #change names#
  names(dataframe4)[2] <- 'depth'
  names(dataframe4)[3] <- 'phase'
  names(dataframe4)[4] <- 'thresh'
  
  return(dataframe4)
  
}

data_labsFinal<-phases(diveSurfBind_ord)

#-------------------------------------------------------------------------------
# depth.diff # function to calculate depth difference (m) from previous depth
#-------------------------------------------------------------------------------
depth.diff<-function(dataframe){
  #depth difference#
  diff = rep(NA,nrow(dataframe))
  for(i in 2:nrow(dataframe)){
    a = dataframe$depth[i]
    b = dataframe$depth[i - 1]
    c = a-b
    diff[i] = c
  }
  dataframe2<-cbind(dataframe,diff)
  return(dataframe2)
}
#-------------------------------------------------------------------------------
# time.diff # function to calculate time difference (s) from previous time
#-------------------------------------------------------------------------------
time.diff<-function(dataframe){
  tim = rep(NA,nrow(dataframe))
  for(i in 2:nrow(dataframe)){
    a = dataframe$date[i]
    b = dataframe$date[i - 1]
    c = difftime(a,b, unit="secs")
    tim[i] = c
  }
  dataframe2<-cbind(dataframe,tim)
  return(dataframe2)
}
#-------------------------------------------------------------------------------
# speed.diff # function to calculate speed difference (m/s) from previous speed
#-------------------------------------------------------------------------------
speed<-function(dataframe){
  spd = rep(NA,nrow(dataframe))
  for(i in 2:nrow(dataframe)){
    a = abs(dataframe$diff[i])
    b = dataframe$tim[i]
    if (b==0){
      c=0
    }
    else {
      c = a/b
    }
    if (c== "NaN") {
      c=0
    }
    spd[i] = c
    
  }
  dataframe2<-cbind(dataframe,spd)
  return(dataframe2)
}

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
#-------------------------------------------------------------------------------
# dive.phaselab # function to label phases of dive
#-------------------------------------------------------------------------------
dive.phaselab<-function(dataframe){
  #descent
  p<-dataframe$depth[which.max(dataframe$depth)]
  p2<-p*0.8
  t<-which(dataframe$depth >=p2)
  
  if (length(t)== nrow(dataframe)| p2 <= dataframe[1,2]){
    #for 3 point dives or very shallow dives (0.8Xmax depth < first reading)
    p<-which.max(dataframe$depth)
    p=p-1
    dataframe$descent<- NA
    dataframe[1:p,10]<- "descent"
    
    #ascent
    #find max depth and first occurance of max depth
    x<-dataframe$depth[which.max(dataframe$depth)]
    asc<-dataframe %>% filter(depth == x)
    asc$a <- +(!duplicated(asc$depth, fromLast=T))
    df <- asc[ -c(2:10) ]
    dataframe<-merge(dataframe, df, by = 'date', all.x= TRUE)
    
    dataframe<-dataframe %>%
      mutate(ph = case_when(descent == 'descent' ~ 'descent',
                            depth == x ~ 'bottom', 
                            lag(a) == 1 ~ 'ascent',
                            TRUE ~ as.character(NA)))
    
    dataframe$ph<-na.locf(dataframe$ph,na.rm = FALSE)
    
    dataframe<-dataframe %>%
      mutate(ph2 = case_when(ph == 'descent' ~ 'descent',
                             ph == 'ascent' ~ 'ascent', 
                             depth == x ~ 'bottom',
                             TRUE ~ 'n'))
    
    dataframe<-dataframe %>%
      mutate(ph3 = case_when(ph2 == 'n' & direction == 'positive' ~ 'minidesc',
                             ph2 == 'n' & direction == 'negative' ~ 'miniasc',
                             ph2 == 'n' & direction == 'no change' ~ 'hangout',
                             ph2 == 'descent'~ 'descent',
                             ph2 == 'ascent'~ 'ascent',
                             ph2 == 'bottom'~ 'bottom',
                             TRUE ~ as.character(NA)))
    
    
    dataframe <- dataframe[ -c(10:13) ]
    dataframe <- dataframe[ -c(6:8) ]
    
    dataframe<-depth.diff(dataframe)
    dataframe<-time.diff(dataframe)  
    dataframe<-speed(dataframe)
    return(dataframe)
  }
else {
  p3<-which(dataframe$depth >=p2)[1]
  p4=p3-1
  dataframe$descent<- NA
  dataframe[1:p4,10]<- "descent"
  
  #ascent
  #find max depth and first occurance of max depth
  x<-dataframe$depth[which.max(dataframe$depth)]
  asc<-dataframe %>% filter(depth <= (x*0.8) & is.na(descent))
  asc$a <- "ascent"
  df <- asc[ -c(2:10) ]
  dataframe<-merge(dataframe, df, by = 'date', all.x= TRUE)
  
  dataframe<-dataframe %>%
    mutate(ph = case_when(descent == 'descent' ~ 'descent',
                          depth >= p2 ~ 'bottom', 
                          a == 'ascent' ~ 'ascent',
                          TRUE ~ as.character(NA)))
  
  dataframe$ph<-na.locf(dataframe$ph,na.rm = FALSE)
  
  dataframe<-dataframe %>%
    mutate(ph2 = case_when(ph == 'descent' ~ 'descent',
                           ph == 'ascent' ~ 'ascent', 
                           depth >= p2 ~ 'bottom',
                           TRUE ~ 'n'))
  
  
  ##mini peaks##
  bs<-which(dataframe$depth >=p2)[1]
  
  be<-rev(which(dataframe$depth >=p2))[1]
  
  b<-dataframe[bs:be,]
  
  
  dataframe2<-b %>%
    mutate(ph3 = case_when(depth >= p2 ~ 'bottom',
                           ph2=='ascent' & direction == 'positive' ~ 'minidesc',
                           ph2 == 'ascent' & direction == 'negative' ~ 'miniasc',
                           ph2 == 'ascent' & direction == 'no change' ~ 'hangout',
                           TRUE ~ as.character(NA)))
  
  
  dataframe2.1<-dataframe2[,c(1,14)]
  
  
  dataframe3<-merge(dataframe, dataframe2.1, by = 'date', all.x= TRUE)
  
  dataframe4<-dataframe3 %>%
    mutate(ph2.1 = case_when(ph2=='ascent' & is.na(ph3) ~ 'ascent',
                             ph2 == 'descent' &  is.na(ph3) ~ 'descent',
                             ph3 == 'bottom' ~ 'bottom',
                             ph3 == 'miniasc' ~ 'miniasc',
                             ph3 == 'minidesc' ~ 'minidesc',
                             ph3 == 'hangout' ~ 'hangout',
                             TRUE ~ as.character(NA)))
  
  
  
  dataframenew<-dataframe4[c("date","depth","phase","thresh","diveno","direction","ph2.1", "diff", "tim", "spd")]
  names(dataframenew)[names(dataframenew) == "ph2.1"] = "ph3"
  
  dataframenew<-depth.diff(dataframenew)
  dataframenew<-time.diff(dataframenew)  
  dataframenew<-speed(dataframenew)
  
  dataframenew<- dataframenew[-c(8:10)]
  return(dataframenew)
}
}

#if transmission stops during the last dive then remove 'incomplete' dive
last.dive.name<-tail(names, n=1)
last.dive<-eval(parse(text = paste0(last.dive.name)))
if(tail(last.dive$depth,n=1)>=customDive_threshold){
 names<-names[1:length(names)-1]
}else{
  names<-names
}

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
#-------------------------------------------------------------------------------
# split.hangouts # function to split 'hangouts' if possible for turtle to reach surface based on trajectory. Note- Hangouts indicate ascent, followed by time spent stationary at depth and then descent
#-------------------------------------------------------------------------------
split.hangouts<- function(dataframe){
  
  if (length(which(dataframe$ph3 == "hangout" & dataframe$depth <10))>0){
    
    dataframe$Group <- NA
    # Find the indices where 'hangout' occurs
    hangout_indices <- which(dataframe$ph3 == "hangout"& dataframe$depth <10)
    
    # Assign group numbers to consecutive occurrences of 'hangout'
    group_num <- 1
    for (i in 1:length(hangout_indices)) {
      dataframe$Group[hangout_indices[i]:nrow(dataframe)] <- group_num
      if (i < length(hangout_indices) && hangout_indices[i + 1] != hangout_indices[i] + 1) {
        group_num <- group_num + 1
      }
    }
    dataframe$Group[!dataframe$ph3 %in% "hangout"] <- NA
    
    # Change 'Group' to NA for groups with more than 1 consecutive occurrences  
    df_filt <- as.data.frame(dataframe %>%
                               group_by(Group) %>%
                               mutate(group_count = sum(ph3 == "hangout")) %>%
                               mutate(Group = ifelse(ph3 == "hangout" & group_count > 1, NA, Group)) %>%
                               select(-group_count))
    
    
    hangout_indices <- which(df_filt$ph3 == "hangout")
    
    # Identify the indices where depth remains the same before 'hangout'
    same_depth_before_hangout <- hangout_indices[df_filt$depth[hangout_indices - 1] == df_filt$depth[hangout_indices]]
    
    # Initialize a new column 'Flag' with all FALSE values
    df_filt$Flag <- FALSE
    # Set TRUE for rows in 'Flag' column where depth remains the same before 'hangout'
    df_filt$Flag[same_depth_before_hangout] <- TRUE
    
    df_filt<-as.data.frame(  df_filt %>%
                               group_by(Group) %>%
                               mutate(Group = ifelse(any(Flag), NA, Group)) %>%
                               ungroup())
    
    t<-which(!is.na(df_filt$Group))
    
    
    if (length(t)>1){
      fill=rep(NA, length(t))
      for (a in 1:length(t)){
        #j2<-dataframe[t,2] 
        #depths& times of peak, n-1, n+1
        depthnm1<--dataframe[t[a]-1,2]
        depthpoint<- -dataframe[t[a],2]
        depthnp1<--dataframe[t[a]+1,2]
        
        timenm1<-dataframe[t[a]-1,1]
        timepoint<-dataframe[t[a],1]
        timenp1<-dataframe[t[a]+1,1]
        
        
        #timediff
        t1<-as.numeric(timepoint-timenm1)*60
        t2<-as.numeric(timenp1-timepoint)*60
        
        #depthdiff
        d1<-depthpoint-depthnm1
        d2<-depthpoint-depthnp1
        
        #calculate speed
        #first speed
        s1<-d1/t1
        #second speed
        s2<-d2/t2
        
        #could turtle cross threshold based on speed before/after peak
        distance_poss1<-s1*150
        distance_poss2<-s2*150
        
        distneeded<- -customDive_threshold-depthpoint
        
        if(distance_poss1>=distneeded || distance_poss2>=distneeded ){
          fill[a] = t[a]
        }
        
      }
      
      fill<- fill[!is.na(fill)]
      difcol2<-c(diff(fill),NA) 
      df2<-data.frame(difcol2, fill)
      if(length(which(df2$difcol2==1))>0){
        dfdrop<-df2[-(which(df2$difcol2 %in% "1")),]
      }
      
      fill2<-dfdrop$fill
    }else{
      fill2=t
    }
    
    
    mx<-max(dataframe$depth)
    #need to loop through output from t below
    if (length(fill2)>0){
      fill2<-fill2-1
      
      for (k in 1:length(dataframe)){  
        
        # newdv = rep(NA, nrow(dataframe))
        rwz = fill2
        de=rep(NA, nrow(dataframe))
        newdv<- .POSIXct(character(nrow(dataframe)))
        
        for (i in 1:length(rwz)){
          t1<- dataframe[fill2[i],1] 
          t2 <- t1+150 
          # t<-as.numeric(t1-t1)
          t<-c(0,150)
          d<-c(dataframe[rwz[i],2], 0)
          xout <- seq(t1,t2,by="2 sec") #was 10 secs but changed becasue duplicate confusion
          xout <- as.numeric(xout-t1)
          pred <- aspline(t, d, xout=xout, n=2)  # prediction
          pt <- which(pred$y >= 1.00)# get depths below dive threshold
          pt<-rev(pt)
          td1= xout[pt[1]] 
          pt1d<-pt[1]
          
          date<-c(t1+td1)
          depth<- pred$y[pt1d]
          
          newdv[rwz[i]] = as.POSIXct(date)
          de[rwz[i]] = depth 
          
          #ab.thresh#
          pt2<- which(pred$y >= customDive_threshold)
          pt2<-rev(pt2)
          td1.1= xout[pt2[1]] 
          pt2d<-pt2[1]
          date2<-c(t1+td1.1)
          depth2<- pred$y[pt2d]
          
          
          newdv[rwz[i]+1] = as.POSIXct(date2)
          de[rwz[i]+1] = depth2
          
          ###########
          ##descent again, make sure it uses rwz[i], NOT J+1
          tda<- date
          tdb <- dataframe[rwz[i]+1,1] 
          ta<-0 #as.numeric(tdb-tda) #   difftime(tdb,tda, unit="secs"))
          ta<-c(ta,150)
          da<-c(dataframe[rwz[i]+1,2], depth)
          xout2 <- seq(tda, tdb,by="2 sec") #was 10 secs but changed becasue duplicate confusion
          xout2 <- as.numeric(xout2-tda)
          pred2 <- aspline(ta, da, xout=xout2, n=2)  # prediction
          pta <- which(pred2$y >= customDive_threshold)# get depths below dive threshold
          pta<-rev(pta)
          td1a= xout[pta[1]] 
          pt1da<-pta[1]
          
          date3<-c(tdb-td1a)
          depth3<- pred2$y[pt1da]
          newdv[rwz[i]+2] = as.POSIXct(date3)
          de[rwz[i]+2] = depth3
        }
        
        dataframe.1<- cbind(dataframe,de,newdv)
        
        df <- dataframe.1[ -c(1:10) ]
        names(df)[names(df) == "newdv"] <- "date"
        names(df)[names(df) == "de"] <- "depth"
        
        
        #final DF!
        multdiv2 <- merge(dataframe,df, by=c("date", "depth"), all = TRUE)
        multdiv2<- multdiv2[complete.cases(multdiv2[ , 1:2]),]
        multdiv2.2<-na.locf(multdiv2, na.rm = FALSE)
        multdiv2.2<-select(multdiv2.2, -c(Group.x, Group.y))
        
        return(multdiv2.2)
      }
    }
  }else{
    
  }
} 

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
  
write.csv(final_dives, "C:/Users/user/Desktop/tort/test/final_forprocess_181762.csv")
################################################################################
#------------------------------------------------------------------------------#
# Step 13. Relabel all data#
#------------------------------------------------------------------------------#
################################################################################

proc_data<- read.csv("C:/Users/user/Desktop/tort/test/final_forprocess_181762.csv")

##grep to add 00's to midnight
proc_data$date[grep("[0-9]{4}-[0-9]{2}-[0-9]{2}$",proc_data$date)] <- paste(
  proc_data$date[grep("[0-9]{4}-[0-9]{2}-[0-9]{2}$",proc_data$date)],"00:00:00")

proc_data$date <- as.POSIXct(proc_data$date,format="%Y-%m-%d %H:%M:%S",tz="CET") ## Your dates need to be in 
proc_data<-proc_data[order(proc_data$date),]

#leave datetime and depth
#proc_data <- proc_data[ -c(1) ]

proc_data<-subset(proc_data, select = c(date,depth))

#reprocess
#rounded 5 in phases function
proc_dataPhase<-phases(proc_data)

#add in depth, time and speed diff from previous value

proc_dataPhase<-depth.diff(proc_dataPhase)
proc_dataPhase<-time.diff(proc_dataPhase)  
proc_dataPhase<-speed(proc_dataPhase)

#identify and label changes in direction, negative = decrease in depth
direc<-proc_dataPhase %>%
  mutate(direction = case_when(diff > 0 ~ 'positive',
                               diff < 0 ~ 'negative',
                               diff == 0 ~ 'no change',
                               TRUE ~ ''))
#remove duplicates
direc<-direc[!duplicated(direc), ]

################################################################################
#==============================================================================#
# Step 14. Compile final dives###
#==============================================================================#
################################################################################
#split dives by number
data_sep<-split(direc, direc$diveno)


#find total number of dives
lastDive<-tail((names(data_sep)), n=1)
divNum<-as.numeric(lastDive[1])

#print names of all dives
names<-sprintf("dive%d", 1:divNum) 

split_names <- as.character(names)

#make each dive a seperate dataframe
for (i in 1:length(data_sep)) {        # Run for-loop
  assign(split_names[i], data_sep[[i]])
}

last.dive.name<-tail(names, n=1)
last.dive<-eval(parse(text = paste0(last.dive.name)))
if(tail(last.dive$depth,n=1)>=customDive_threshold){
  names<-names[1:length(names)-1]
}

#for loop for each dive and get summary statistic
#did have phase lab function here, but already run previously
#apply to all dives
for(i in  1:length(names)){
  data =  eval(parse(text = paste0( names[i])))
  assign(paste(names[i],"phase",sep = "_"), t<-dive.phaselab(data))
  print(paste(names[i],"phase",sep = "_"))
}

##new data frame with all dives
dive_list = mget(ls(pattern = "_phase"))
dives_complete = do.call(what = rbind, args = dive_list)
dives_complete<-dives_complete[order(dives_complete$date),]

write.csv(dives_complete, "C:/Users/user/Desktop/tort/test/processed_dives_181762.csv")
################################################################################
###----------------------------------------------------------------------------#
### Step 15. Extract dive summaries #
###----------------------------------------------------------------------------#
################################################################################
#-------------------------------------------------------------------------------
# dive_summary # function to extract dive summaries from individual dives
#-------------------------------------------------------------------------------
dive_summary<-function(dataframe){
  #descent time
  bstart<-which(dataframe$ph3=='bottom')
  bstart<-bstart[1]
  descensionTime<-as.numeric(difftime(dataframe[bstart,1], dataframe[1,1], units="secs" ))
  #descent distance
  descd<- rev(which(dataframe$ph3=='descent'))[1]
  descd2<-descd+1
  desc_d<-abs(dataframe[c(1:descd2),8])
  descensionDistance<-sum(desc_d,na.rm=TRUE)
  
  
  maximumDepth<-dataframe$depth[which.max(dataframe$depth)]
  #regre<-0.5172* maximumDepth - 0.4106 from sandras 1999 green turtle paper
  #mode
  v<-dataframe$depth
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  mode<-getmode(dataframe$depth)
  ###
  
  fluc<- filter(dataframe,ph3=='minidesc'| ph3=='miniasc'| ph3=='hangout')
  
  if (nrow(fluc)>0) {
    
    peakfin<- rev(which(dataframe$ph3=='minidesc' | dataframe$ph3=='miniasc'| dataframe$ph3=='hangout'))[1]
    peakfin2<-peakfin+1
    peakstart<- which(dataframe$ph3=='minidesc' | dataframe$ph3=='miniasc'| dataframe$ph3=='hangout')[1]
    
    peak_d<-abs(dataframe[c(peakstart:peakfin2),8])
    peak_dist<-sum(peak_d)
    
    peakTime<-sum(fluc$tim)
    peakHeight<-fluc$depth[which.min(fluc$depth)] #height of peak
    #peak_bot<-fluc$depth[which.max(fluc$depth)] #lowest part of peak
    
    #bottom phase
    bfin<- rev(which(dataframe$ph3 == 'bottom'))[1]
    bottomTime<-as.numeric(difftime(dataframe[bfin,1], dataframe[bstart,1], units= "secs"))
    bottomTime<- bottomTime-peakTime
    bot2<- filter(dataframe,ph3=='bottom')
    botd<-abs(diff(bot2$depth)) #added abs here
    bottomDistance<-sum(botd)
    
    #count number of peaks
    
    dataframe<-dataframe %>%
      mutate(peak_num = case_when(lag(ph3) == "miniasc" & ph3 =="minidesc" ~ '1',
                                  lag(ph3) == "hangout" & ph3 =="minidesc" ~ '1',
                                  lag(ph3) == "minidesc" & ph3 =="miniasc" ~ '1',
                                  TRUE ~ '0'))
    
    dataframe$peak_num<-as.numeric(dataframe$peak_num)
    
    numberOfPeaks<-sum(dataframe$peak_num)
    
    last<-tail(fluc, n=1)
    last[1,1]
    last1<-which(dataframe[,1]== last[1,1])
    last2<-last1+1
    depthafter<-dataframe[last2,2]
    
    #prop of max depth depth after peak is at
    proportionMaxDepthAfterPeak<- depthafter/maximumDepth
    peakSD<- sd(fluc$depth)
    
    
    
  }else {
    peakTime<-NA
    peakHeight<-NA
    # peak_bot<-NA
    numberOfPeaks<-0
    #bottom phase
    bfin<- rev(which(dataframe$ph3 == 'bottom'))[1]
    bottomTime<-as.numeric(difftime(dataframe[bfin,1], dataframe[bstart,1], units= "secs"))
    bot2<- filter(dataframe,ph3=='bottom')
    botd<-abs(diff(bot2$depth))
    bottomDistance<-sum(botd)
    proportionMaxDepthAfterPeak<- NA
    peakSD<- NA
  }
  
  #ascent time
  ascfin<-rev(which(dataframe$ph3 == 'ascent'))[1]
  ascensionTime<-as.numeric(difftime(dataframe[ascfin,1], dataframe[bfin,1], units= "secs"))
  ascent<-dataframe %>% filter(ph3== "ascent")
  ascensionDistance<-sum(abs(ascent$diff),na.rm=TRUE)
  #change point
  a<- filter(dataframe,ph3=='ascent')
  if (nrow(a)>2) {
    asc_cp = cpt.mean((a$spd*100), method="PELT") #mean changepoints using PELT
    numberAscentionChangepoints<-ncpts(asc_cp)
  }else {
    numberAscentionChangepoints<-0
  }
  #add in depth of changepoints
  
  if (numberAscentionChangepoints > 0){
    points<-cpts(asc_cp)
    depth_acp<-a[points,2]
    
    depthAscentionChangepoint1<-depth_acp[1]
    depthAscentionChangepoint2<-depth_acp[2]
    secondAscentionChangepointeighty<-as.integer(depthAscentionChangepoint2> depthAscentionChangepoint1*0.8)
    
  }else {
    depthAscentionChangepoint1<-NA
    depthAscentionChangepoint2<-NA
    secondAscentionChangepointeighty<-NA
  }
  #descent change points
  w<-which(dataframe$ph3 == "descent") 
  if (length(w)>1) {
    first_bot<-w[1]
    desc_filt<-dataframe[2:first_bot,]
    desc_cp = cpt.mean(desc_filt$depth, method="PELT") #mean changepoints using PELT
    numberDescensionChangepoints<-ncpts(desc_cp)
  }else {
    numberDescensionChangepoints<-0
  }
  
  ##depth of desc change points##
  if (numberDescensionChangepoints > 0){
    points2<-cpts(desc_cp)
    depth_dcp<-a[points2,2]
  }else {
    depth_dcp<-NA
  }
  
  diveStartTime<-dataframe[1,1]
  
  #total dive time
  divefin<- rev(dataframe[,1])
  diveEndTime<-divefin[1]
  diveTime<-as.numeric(diveEndTime-diveStartTime)*60
  
 
  nr<- nrow(dataframe)
  maxround<-round(maximumDepth, digits=1) #can get rid of for entire dive cp
  
  if (nr>2 & maxround>customDive_threshold){
    nodesc<- dataframe[dataframe$ph3 != "descent", ]  # #can get rid of for entire dive cp
    
    tot_cp = cpt.mean((nodesc$spd*100), method="PELT") #mean changepoints using PELT
    totalChangepoints<-ncpts(tot_cp)
  }else{
    totalChangepoints=0
  }
  #flag dives where rounded max depth =5
  if(maxround<=customDive_threshold){
    exclude="Y"
  } else{
    exclude="N"
  }
  
  #flag platue in ascent
  plat<-which(a$direction == "no change")
  if (length(plat) > 0 ){
    flatAscentionChangepoint<-"1"
  } else {
    flatAscentionChangepoint<-"0"
  }
  
  #if down again in ascent phase
  down1<-which(a$direction == "positive")
  #if decrease is larger than X meters then dont
  down2<-which(a[down1,8]>2) #change the 2 here for threshold values
  
  if (length(down2) > 0 ){
    descentAfterStatAscent<-length(down2)
  } else {
    descentAfterStatAscent<-"0"
  }
  
  #is difference first part of ascent> asc1:2??
  differenceBetweenFlatAscentPoints<-as.integer(abs(a[1,8])> abs(a[2,8]))
  
  #is first ascent point in 60% of dive depth
  firstAscentionPointSixtyMaxDepth<-as.integer(a[1,2]<0.6*maximumDepth)
  
  #timing of maximumDepth
  max_time<-as.numeric(diveEndTime-dataframe[which.max(dataframe$depth),1])*60
  maxDepthPropOfDive<-max_time/diveTime *100
  
  #bin dives by number of points in each quater of depth
  quat<-maximumDepth/5 #was 4
  quaters<-seq(0,maximumDepth, by= quat)
  bin<-dataframe %>% mutate(depth_bin = cut(depth, breaks=c(quaters)))
  fre<-as.data.frame(table(bin$depth_bin))
  numberOfReadingsFirstFifth<-fre[1,2]
  numberOfReadingsSecondFifth<-fre[2,2]
  numberOfReadingsThirdFifth<-fre[3,2]
  numberOfReadingsFourthFifth<-fre[4,2]
  numberOfReadingsFifthFifth<-fre[5,2]
  
  nm<-c("numberOfReadingsFirstFifth", "numberOfReadingsSecondFifth", "numberOfReadingsThirdFifth", "numberOfReadingsFourthFifth", "numberOfReadingsFifthFifth")
  totpt<-c(numberOfReadingsFirstFifth,numberOfReadingsSecondFifth,numberOfReadingsThirdFifth,numberOfReadingsFourthFifth,numberOfReadingsFifthFifth)
  prop2<-data.frame(nm,totpt)
  whichFifthTime<-prop2[which.max(prop2$totpt),1]
  #whichFifthTime= what prop spent most of dive
  #under 10
  maxDepthUnderTenM<-as.integer(maximumDepth> 10)
  
  descent<-dataframe %>% filter(ph3== "descent")
  
  botsd<-dataframe%>% filter(!ph3== "descent" & !ph3== "ascent")#anything not asc or desc
  
  bottom<-dataframe %>% filter(ph3== "bottom")
  
  descentSD<-sd(descent$depth)
  restSD<-sd(botsd$depth) #anything not asc or desc
  ascentSD<-sd(ascent$depth)
  bottomSD<- sd(bottom$depth)
  
  ##indication of if ascent traj curved

  
  #extrapolate straight line for ascent
  if(nrow(a)>2 & nrow(bot2)>0){
    t1<- dataframe[bfin,1] 
    t2 <- dataframe[ascfin,1] 
    
    dep1<-dataframe[bfin,2]
    dep2<-dataframe[ascfin,2]
    
    xout <- seq(t1,t2,by="60 sec") #maybe change to 300 secs
    
    xout2 <- as.numeric(xout-t1)
    
    t<-c(xout2[1], rev(xout2)[1])
    
    d<-c(-dep1, dep2)
    
    pred <- aspline(t, d, xout=xout2)  # prediction
    round(pred$y, digits=2)
    
    newdf<-cbind(xout,round(pred$y,digits=2))
    
    df<-as.data.frame(xout)
    dep<-round(pred$y,digits=2)
    newdf<-cbind(df, abs(dep))
    names(newdf)[1] <- "date"
    #remove first and last values as dont care
    newdf2<-newdf[-(1),]
    newdf3<-slice(newdf2, 1:(n() - 1))  
    
    #merge new predictions with existing times
    newdf4<-merge(dataframe, newdf3, by="date")
    
    crv<-newdf4$`abs(dep)`-newdf4$depth
    
    ascentOffsetStraightTrajectory<-mean(crv)
    ascentTrajectoryDirection<-sign(ascentOffsetStraightTrajectory)
    
  }else{
    ascentOffsetStraightTrajectory<- NA
    ascentTrajectoryDirection<-sign(ascentOffsetStraightTrajectory)
  }
  proportionDiveBottom<-bottomTime/diveTime*100
  proportionDiveAscent<-ascensionTime/diveTime*100
  
  if (nrow(a)>1){
    stationaryAscent<-a[1,2]/maximumDepth*100} else{
      stationaryAscent<-NA
    }
  
  #how much of max depth is first ascent?
  #a1/maximumDepth*100
  
  #put all in dataframe
  sum_df<- data.frame(diveStartTime, descensionTime,descensionDistance, bottomTime, ascensionTime, ascensionDistance, maximumDepth, peakTime, peakHeight, numberOfPeaks, mode, numberAscentionChangepoints, numberDescensionChangepoints, depthAscentionChangepoint1, depthAscentionChangepoint2, bottomDistance, diveTime, totalChangepoints, exclude, flatAscentionChangepoint, descentAfterStatAscent, differenceBetweenFlatAscentPoints, maxDepthPropOfDive, numberOfReadingsFirstFifth, numberOfReadingsSecondFifth, numberOfReadingsThirdFifth, numberOfReadingsFourthFifth, numberOfReadingsFifthFifth, maxDepthUnderTenM, whichFifthTime, proportionMaxDepthAfterPeak, descentSD, restSD, ascentSD, peakSD, bottomSD, firstAscentionPointSixtyMaxDepth, ascentOffsetStraightTrajectory, ascentTrajectoryDirection, proportionDiveAscent, proportionDiveBottom, secondAscentionChangepointeighty, stationaryAscent,diveEndTime)
  
  return(sum_df)
  
}

final<- paste(names, "phase", sep = "_")

#apply to all dives

for(i in  1:length(final)){
  data =  eval(parse(text = paste0( final[i])))
  assign(paste(final[i],"final",sep = "_"), t<-dive_summary(data))
  print(paste(final[i],"final",sep = "_"))
}
final2<- paste(final, "final", sep = "_")

fin_list = mget(ls(pattern = "_phase_final"))
dive_summaries = do.call(what = rbind, args = fin_list)
dive_summaries<-dive_summaries[order(dive_summaries$diveStartTime),]


dive_summaries$diveID <- c(1:nrow(dive_summaries))

#replace id with turtleID
dive_summaries$organismID<-181762
dive_summaries$concatID<-paste0("181762","_",dive_summaries$diveID)

#-------------------------------------------------------------------------------
# surf.time # function to calculate surface period following a dive
#-------------------------------------------------------------------------------
surf.time<-function(dataframe){
  tim = rep(NA,nrow(dataframe))
  for(i in 2:nrow(dataframe)){
    a = dataframe$diveStartTime[i]
    b = dataframe$diveEndTime[i - 1]
    c = difftime(a,b, unit="secs")
    tim[i] = c
  }
  tim2<-tim[-1]
  tim3<-append(tim2,NA)
  
  dataframe2<-cbind(dataframe,tim3)
  return(dataframe2)
}

dive_summaries<-surf.time(dive_summaries)
colnames(dive_summaries)[colnames(dive_summaries) == "tim3"] ="surface_int"

write.csv(dive_summaries, "C:/Users/user/Desktop/tort/test/dive_summarys_181762.csv")
head(dive_summaries)
beep()

################################## Finish!######################################
