#calibration benchmarks

#load packages
library(rgdal)
library(rgeos)
library(raster)
library(tidyr)
library(dplyr)
#install.packages('OneR')
library(OneR)
library(ggplot2)
#install.packages('patchwork')
library(patchwork)
#install.packages('hrbrthemes')
library(hrbrthemes)
library(purrr)
library(plyr)

#read in files

#historic files
clim0 <- raster("C:/Users/13146/Documents/ReburnsAK/FULL_Landscape_Inputs_AK_062521_COPYME/AK_ClimateMap_10_Regions.tif")

#dalton_perims <- readOGR("C:/Users/13146/Dropbox/FireData_SCRPPLEParameterization/Fire Perimeters/perim_dalton_buff.shp") ## this is from Geomac
#plot(dalton_perims)
 
fire_hist <- readOGR("C:/Users/13146/Documents/ReburnsAK/SCRPPLE_AK/AKstudyarea_Fires_Lightning.shp")
#dalton_perim_latlong <- spTransform(dalton_perims, crs(fire_hist))
clim0_outline <- as(extent(clim0), 'SpatialPolygons')
crs(clim0_outline) <- crs(clim0)
clim0_outlinelatlong <- spTransform(clim0_outline, crs(fire_hist))

#fix self intersection with gBuff with width=0
#fire_hist <- gBuffer(fire_hist, byid=T, width=0)
#now crop to study landscape
dalton_hist <- crop(fire_hist, (clim0_outlinelatlong))
dalton_hist <- spTransform(dalton_hist, crs(clim0))
plot(dalton_hist);plot(clim0_outline, add=T)


#bring in mtbs files and bring together... as a stack? 
mtbs_dir <- "C:/Users/13146/Dropbox/FireData_SCRPPLEParameterization/MTBS/composite_data/MTBS_BSmosaics/"
clim0_proj <- projectRaster(clim0, crs="+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
#i=1989
possibly_function <- possibly(function(i){
  
          mtbs <- raster(paste0(mtbs_dir, i, "/mtbs_AK_", i, "/mtbs_AK_", i, ".tif")) %>% crop(clim0_proj) %>% resample(clim0_proj, method='ngb')
          #mtbs_all <- stack(mtbs_all, mtbs)
          mtbs_df <- as.data.frame(mtbs) 
          mtbs_df <- data.frame(mtbs_df[complete.cases(mtbs_df),])
          names(mtbs_df) <- "severity"
          return(mtbs_df)
          
        }, otherwise="something wrong here")

mtbs_all <- map(1984:2017, possibly_function)
mtbs_all_df <- ldply(mtbs_all, data.frame)[,1] %>% na.omit()

hist(mtbs_all_df)

#### Fire Severity ####

#combine severity into three classes
# is <- c(1,2,3,4,5,6)
# becomes <- c(1,1,2,2,3,3)
# rcl <- as.matrix(data.frame(is, becomes))
# mtbs_all_df <- as.numeric(mtbs_all_df)
# rcl_sev <- reclassify(mtbs_all_df, rcl)

mtbs_all_df <- data.frame(mtbs_all_df)
severity_df <- mtbs_all_df %>% mutate(mtbs_all_df=recode(mtbs_all_df, '1'='1', '2'='1', '3'='2', '4'='2', '5'='3', '6'='3'))
hist(as.numeric(severity_df$mtbs_all_df))
severity_df$Name <- "Historical"
names(severity_df)[1] <- "Severity"
write.csv(severity_df, "C:/Users/13146/Dropbox/NSF ReburnsAK LANDIS Results/landscape/R/Fire_Stats/OG_full_landscape_severity_df.csv")

####fire rotation (ie total time / proportion of area burned during that time)####


#### Limit to fires over 4ha####
dalton_hist_over4 <- subset(dalton_hist, dalton_hist$FIRE_SIZE>=9.88422)
#dalton_hist_over4 <- dalton_hist

###Historic fire rotation
cells_per_ha <-4
#rasterize dalton_hist
timeperiod_hist <- max(as.numeric(paste(dalton_hist_over4$FIRE_YEAR))) - min(as.numeric(paste(dalton_hist_over4$FIRE_YEAR)))
totalarea <- ncell(clim0)*cells_per_ha #gives total area in hectares

burned_acres <- sum(dalton_hist_over4$FIRE_SIZE)
burned_hectares <- burned_acres * 0.404686

rotation_hist <- timeperiod_hist/(burned_hectares/totalarea)
rotation_hist 

####Number of fires ####

###Historic number of fires
#total fires
total_fires_hist <- nrow(unique(cbind(dalton_hist_over4$FIRE_NAME, dalton_hist_over4$LATITUDE, dalton_hist_over4$LONGITUDE, dalton_hist_over4$FIRE_YEAR)))

##fires per year
firesperyr_hist <- total_fires_hist/timeperiod_hist
firesperyr_hist 

#fires per year variability
head(dalton_hist_over4)
dalton_hist_over4$FIRE_YEAR <- as.numeric(dalton_hist_over4$FIRE_YEAR)
#dalton_hist_over4 <- subset(dalton_hist_over4, dalton_hist_over4$FIREYEAR >1969 & dalton_hist_over4$FIREYEAR < 1991)
annual_fires <- NULL
for (i in 1:length(unique(dalton_hist_over4$FIRE_YEAR))) {
  #i=48
  fires <- subset(dalton_hist_over4, dalton_hist_over4$FIRE_YEAR == unique(dalton_hist_over4$FIRE_YEAR)[i])
  num_fires <- nrow(unique(cbind(fires$FIRE_NAME, fires$LATITUDE, fires$LONGITUDE, fires$FIRE_YEAR)))
  ha_burned <- sum(fires$FIRE_SIZE)*0.404686 #convert to ha
  year <- as.numeric(paste(fires$FIRE_YEAR[1]))
  num_fires <- c(year, num_fires, ha_burned)
  annual_fires <- data.frame(rbind(annual_fires, num_fires))
  
}
#add in 0 fires years
names(annual_fires) <- c("Year", "Num_Fires", "Ha_burned")
annual_fires$Year <- as.numeric(paste(annual_fires$Year))
annual_fires <- annual_fires[order(annual_fires$Year),]
#find missing years
for (i in 1:nrow(annual_fires)){
  #i=4
  if (annual_fires$Year[i + 1] != annual_fires$Year[i] + 1){
    
    missing <- data.frame((annual_fires$Year[i] + 1), 0, 0)
    names(missing) <- c("Year", "Num_Fires", "Ha_burned")
    annual_fires <- rbind(annual_fires, missing)
    annual_fires <- annual_fires[order(annual_fires$Year),]
  } 
  
}

#Num fires per year
firesperyr_hist <- mean(annual_fires$Num_Fires)
firesperyr_hist_sd <- sd(annual_fires$Num_Fires)
#Hectares burned per year
haperyr_hist <- mean(annual_fires$Ha_burned)
haperyr_hist_sd <- sd(annual_fires$Ha_burned)


###Historic fire size stats
dalton_hist_df <- as.data.frame(dalton_hist_over4)
dalton_hist_df$hectares <- dalton_hist_df$FIRE_SIZE * 0.404686
meansize_hist <- mean(dalton_hist_df$hectares)
meansize_hist

sdsize_hist <- sd(dalton_hist_df$hectares)

maxsize_hist <- max(dalton_hist_df$hectares)
maxsize_hist


#histogram of fire size in hectares
hist(dalton_hist_df$hectares)
hist_histogram_df <- dalton_hist_df[,c("FOD_ID","hectares")]
hist_histogram_df$Name <- "Historical"
names(hist_histogram_df) <- c("ID", "Hectares", "Name")

#create hectare bins
#options(scipen=999)
#dalton_hist_df$hectare_bins <- cut(dalton_hist_df$hectares, breaks=50)
#head(dalton_hist_df)

# #getting uncertainty metrics for fire rotation. Calculate on a decadal basis
# #rotation_hist <- timeperiod/(burnedarea_hist/totalarea)
# #rotation_hist #146 years
# 
# decadal_rotation <- NULL
# #years <- c(1946, 1956, 1966, 1976, 1986, 1996, 2006)
# for (i in 1992:2015) {
#   #i=1956
#   decade <- subset(annual_fires, annual_fires$Year >= i & annual_fires$Year < i + 10)
#   onerotation <- 10/((sum(decade$Ha_burned))/totalarea)
#   decadal_rotation <- rbind(onerotation, decadal_rotation)
#   
# }
# decade_rotation <- mean(decadal_rotation)
# decade_rotation_sd <- sd(decadal_rotation)

#summary table with proportions
dalton_hist_size_summ <- hist_histogram_df %>% 
  na.omit() %>%
  group_by(Hectares) %>%
  dplyr::summarise(count=n()) 

max(dalton_hist_size_summ$Hectares)
breaks <- 100 * c(0:5570.91)
dalton_hist_size_summ$hectbins <- cut(dalton_hist_size_summ$Hectares, breaks=breaks, lables=TRUE)
dalton_hist_size_grouped <- stats::aggregate(dalton_hist_size_summ$count, by=list((dalton_hist_size_summ$hectbins)), FUN=sum)

dalton_hist_size_grouped$percent <- (dalton_hist_size_grouped$x/(sum(dalton_hist_size_grouped$x)))*100
dalton_hist_size_grouped #summary table showing distribution
names(dalton_hist_size_grouped) <- c("Bins", "Count", "Percent")
dalton_hist_size_grouped$Bins <- as.character(dalton_hist_size_grouped$Bins)
tags <- strsplit(dalton_hist_size_grouped$Bins, c(",")) %>% unlist() %>% matrix(ncol=2, byrow=TRUE) %>% as.data.frame()
tags <- tags[,1] 
tags <- sub('.', '', tags) %>% as.numeric()
dalton_hist_size_grouped$Bins <- tags
plot(dalton_hist_size_grouped$Percent ~ dalton_hist_size_grouped$Bins)
dalton_hist_size_grouped$Name <- "Historical"

#mean(dalton_hist_df$hectares)
#min(dalton_hist_df$hectares)

fire.stat.hist <- data.frame(firesperyr_hist, firesperyr_hist_sd, haperyr_hist, haperyr_hist_sd, maxsize_hist, meansize_hist, 
                             sdsize_hist, decade_rotation, decade_rotation_sd)

colnames(fire.stat.hist) <- c("AVG num of fires/yr", "SD num of fires/yr", "AVG ha burned/yr", "SD ha burned/yr", 
                         "Max ha burned", "Mean ha/fire", "SD ha/fire", "Fire Rotation", "Fire Rotation SD")

fire.stat.hist$Name <- c("Historic")
write.csv(fire.stat.hist, "C:/Users/13146/Documents/Fire_History_Stats_FULL_062221.csv")

####SIMULATION RESLUTS####

##This calculates basic fire stats from the DFFS summary log file
##These data re to be compared to historical stats for calibration purposes
#Alec wrote this script and it was modified by ML

fire_dat_original<- read.csv("C:/Users/13146/Documents/ReburnsAK/SCRPPLE_AK/Run_061821/scrapple-summary-log.csv")
fire_dat <- fire_dat_original[,c("SimulationYear","TotalBurnedSitesLightning", "NumberFiresLightning", "TotalBiomassMortalityLightning")]

AvgFireSizes <- (fire_dat[,"TotalBurnedSitesLightning"]* cells_per_ha)/fire_dat[,"NumberFiresLightning"]
TotalFireSizes <- (fire_dat[,"TotalBurnedSitesLightning"]* cells_per_ha)
all.firedat <- cbind(fire_dat, AvgFireSizes, TotalFireSizes)
all.firedat <- all.firedat[complete.cases(all.firedat),]

fire_events <- read.csv("C:/Users/13146/Documents/ReburnsAK/SCRPPLE_AK/Run_061821/scrapple-events-log.csv")
firedat_histogram_df <- fire_events[,c(1,10)]

firedat_histogram_df$TotalSitesBurned <- firedat_histogram_df$TotalSitesBurned * cells_per_ha
firedat_histogram_df$Name <- "Simulated"
names(firedat_histogram_df) <- c("ID", "Hectares", "Name")
firedat_histogram_df$ID <- as.character(firedat_histogram_df$ID)

#firedat_histogram_summ <- firedat_histogram_df %>% 
#  group_by(Hectares) %>% 
#  summarise(count=n()) 

firedat_histogram_summ <- firedat_histogram_df

#max(firedat_histogram_summ$Hectares)
#breaks <- 100 * c(0:43416)
breaks <- c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130)

firedat_histogram_summ$hectbins <- base::cut(firedat_histogram_summ$Hectares, breaks=breaks)
firedat_histogram_grouped <- stats::aggregate(firedat_histogram_summ$Hectares, by=list(firedat_histogram_summ$hectbins), FUN=sum)

firedat_histogram_grouped$percent <- (firedat_histogram_grouped$x/(sum(firedat_histogram_grouped$x)))*100
#firedat_histogram_grouped #summary table showing distribution
names(firedat_histogram_grouped) <- c("Bins", "Count", "Percent")
firedat_histogram_grouped$Bins <- as.character(firedat_histogram_grouped$Bins)
tags <- strsplit(firedat_histogram_grouped$Bins, c(",")) %>% unlist() %>% matrix(ncol=2, byrow=TRUE) %>% as.data.frame()
tags <- tags[,1] 
tags <- sub('.', '', tags) %>% as.numeric()
firedat_histogram_grouped$Bins <- tags
#plot(firedat_histogram_grouped$Percent ~ firedat_histogram_grouped$Bins)
firedat_histogram_grouped$Name <- "Simulated"

#Manually adjust to years the sim ran and the resolution of your cell size.
sim_years<- max(fire_dat_original$SimulationYear)

#summary of all fire data
area.sum <- sum(all.firedat$TotalFireSizes)
avg.area <-area.sum/sim_years ##avg sites burned per year
tot.fires<- sum(all.firedat$NumberFiresLightning)
avg.fires <- tot.fires/sim_years  ##average # of fires per year
fire.max <- max(all.firedat$TotalBurnedSites* cells_per_ha)
fire.mean <- mean(all.firedat$AvgFireSizes)
rotation <- sim_years/(area.sum/totalarea) #fire rotation
fire.stat<-data.frame(avg.fires, 0, avg.area, 0, fire.max, fire.mean, 0, rotation, 0)

colnames(fire.stat) <- c("AVG num of fires/yr", "SD num of fires/yr", "AVG ha burned/yr", "SD ha burned/yr",
                         "Max ha burned", "Mean ha/fire", "SD ha/fire", "Fire Rotation", "Fire Rotation SD")

fire.stat$Name <- "Simulated"

#severities
sim_severity <- fire_dat_original[,11:13]
#head(sim_severity)
names(sim_severity) <- c("1","2","3")
#wide to long
sim_severity <- gather(sim_severity, severity, cells)
sim_severity_new <- NULL
for (i in 1:nrow(sim_severity)) {
  #i=88
  sim_severity_new <- rbind(sim_severity_new, sim_severity[rep(i, sim_severity[i,2]),])
  
}
sim_severity_new$Name <- "Simulated"
sim_severity_new <- sim_severity_new[,c(1,3)]
sim_severity_new$severity <- as.numeric(sim_severity_new$severity)
names(sim_severity_new)[1] <- "Severity"
#head(sim_severity_new)

#comparison charts 
compare <- rbind(fire.stat, fire.stat.hist)

compare_histogram <- rbind(dalton_hist_size_grouped, firedat_histogram_grouped)

compare_severity <- rbind(sim_severity_new, severity_df)

names(compare) <- c("FiresPerYr", "FiresPerYr_sd", "HAburnedPerYr", "HAburnedPerYr_sd", "MaxHABurned", "AvgHAPerFire", "SDHAPerFire", "FireRotation", "FireRotationSD", "Name")


if ((fire.stat.hist$`AVG num of fires/yr` - fire.stat.hist$`SD num of fires/yr`)>=0) {
  ymin=(fire.stat.hist$`AVG num of fires/yr` - fire.stat.hist$`SD num of fires/yr`)} else {ymin=0}
compare_FiresPerYr <- ggplot(data=compare, aes(fill=Name, y=FiresPerYr, x=Name)) + 
  geom_bar(stat='identity', color="#e9ecef") +
  geom_errorbar(aes(x=Name, y=FiresPerYr, ymin=ymin, ymax=FiresPerYr + FiresPerYr_sd), width=0.25, alpha=0.5, size=.8) +
  ggtitle("Average # Fires / Year") +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() +
  theme(legend.position = "none")

if ((fire.stat.hist$`AVG ha burned/yr` - fire.stat.hist$`SD ha burned/yr`)>=0) {
  ymin=(fire.stat.hist$`AVG ha burned/yr` - fire.stat.hist$`SD ha burned/yr`)} else {ymin=0}
compare_HAburnedPerYr <- ggplot(data=compare, aes(fill=Name, y=HAburnedPerYr, x=Name)) + 
  geom_bar(stat='identity', color="#e9ecef") + 
  geom_errorbar(aes(x=Name, y=HAburnedPerYr, ymin=ymin, ymax=HAburnedPerYr + HAburnedPerYr_sd), width=0.25, alpha=0.5, size=.8) +
  ggtitle("Average Hectares Burned / Year") +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() +
  theme(legend.position = "none")

compare_MaxHABurned <- ggplot(data=compare, aes(fill=Name, y=MaxHABurned, x=Name)) + 
  geom_bar(stat='identity', color="#e9ecef") + 
  ggtitle("Maximum Fire Size (Ha)") +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() +
  theme(legend.position = "none")

if ((fire.stat.hist$`Mean ha/fire` - fire.stat.hist$`SD ha/fire`)>=0) {
  ymin=(fire.stat.hist$`Mean ha/fire` - fire.stat.hist$`SD ha/fire`)} else {ymin=0}
compare_AvgHAPerFire <- ggplot(data=compare, aes(fill=Name, y=AvgHAPerFire, x=Name)) + 
  geom_bar(stat='identity', color="#e9ecef") + 
  geom_errorbar(aes(x=Name, y=AvgHAPerFire, ymin=ymin, ymax=AvgHAPerFire + SDHAPerFire), width=0.25, alpha=0.5, size=.8) +
  ggtitle("Average Fire Size (Ha)") +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() +
  theme(legend.position = "none")

if ((fire.stat.hist$`Fire Rotation` - fire.stat.hist$`Fire Rotation SD`)>=0) {
  ymin=(fire.stat.hist$`Fire Rotation` - fire.stat.hist$`Fire Rotation SD`)} else {ymin=0}
compare_FireRotation <- ggplot(data=compare, aes(fill=Name, y=FireRotation, x=Name)) +
  geom_bar(stat='identity', color="#e9ecef") + 
  geom_errorbar(aes(x=Name, y=FireRotation, ymin=ymin, ymax=FireRotation + FireRotationSD), width=0.25, alpha=0.5, size=.8) +
  ggtitle("Fire Rotation") +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() +
  theme(legend.position = "none")

compare_histogram$Perc_Rnd <-ceiling(compare_histogram$Percent)
#rebuild frequency df based on percentages
compare_histogram_new <- NULL
for (i in 1:nrow(compare_histogram)) {
  #i=88
  compare_histogram_new <- rbind(compare_histogram_new, compare_histogram[rep(i, compare_histogram[i,5]),])
  
}

compare_histogram_new <- compare_histogram_new[,c(1,4)]

compare_firesizedistr <- ggplot(compare_histogram_new, aes(x=Bins, fill=Name)) + 
  geom_histogram(color="#e9ecef", alpha=0.6, position = 'dodge') +
  ggtitle("Fire Size Distribution") +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  ylab("Percentage of Fires") +
  theme_ipsum()

#fire severity distribution 
severity_summ <- data.frame(with(compare_severity, table(Severity, Name)))
sev_summ_hist <- subset(severity_summ, severity_summ$Name=="Historical")
sev_summ_hist$Percent <- sev_summ_hist$Freq/sum(sev_summ_hist$Freq)*100
sev_summ_sim <- subset(severity_summ, severity_summ$Name=="Simulated")
sev_summ_sim$Percent <- sev_summ_sim$Freq/sum(sev_summ_sim$Freq)*100
severity_summ <- rbind(sev_summ_hist, sev_summ_sim)

compare_severitydistr <- ggplot(severity_summ, aes(x=Name, y=Percent, fill=Severity)) +
  geom_bar(color="#e9ecef", position='fill', stat='identity') +
  ggtitle("Fire Severity Distribution") +
  scale_fill_manual(values=c("cornsilk", "#69b3a2", "#404080")) +
  ylab("Severity Proportions") +
  theme_ipsum()


#plot altogether                                
((compare_FiresPerYr|compare_HAburnedPerYr)/
  (compare_MaxHABurned|compare_AvgHAPerFire)/
  (compare_FireRotation|compare_severitydistr)) /
  compare_firesizedistr

