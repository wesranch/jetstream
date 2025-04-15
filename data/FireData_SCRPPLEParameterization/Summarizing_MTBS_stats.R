### Fire stats from MTBS then combined with AKLFDB ###
dir <- "/Users/shelbyweiss/Documents/FireData_SCRPPLEParameterization/FullLandscape_11milha_072221/"
clim0 <- raster(paste0(dir, "AK_ClimateMap_10_Regions.tif"))
outline <- readOGR(paste0(dir, "StudyArea_11milha_072221.shp")) 

# bring in mtbs files and bring together
mtbs_dir <- "/Users/shelbyweiss/Documents/FireData_SCRPPLEParameterization/MTBS_BSmosaics/"
clim0_proj <- projectRaster(clim0, crs="+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
+ellps=GRS80 +towgs84=0,0,0")

#i=1990
 possibly_function <- possibly(function(i){
   
   mtbs <- raster(paste0(mtbs_dir, i, "/mtbs_AK_", i, "/mtbs_AK_", i, ".tif")) %>% crop(clim0_proj) %>% resample(clim0_proj, method='ngb')
   mtbs_df <- as.data.frame(mtbs) 
   mtbs_df <- data.frame(mtbs_df[complete.cases(mtbs_df),])
   names(mtbs_df) <- "severity"
   mtbs_df$year <- paste(i)
   return(mtbs_df)
   
 }, otherwise="something wrong here")
 
 mtbs_all <- map(1984:2017, possibly_function)
 mtbs_all_df <- ldply(mtbs_all, data.frame)[,1:2] %>% na.omit()
 
 hist(mtbs_all_df)
 head(mtbs_all_df)
 
saveRDS(mtbs_all_df, paste0(mtbs_dir, "mtbs_all_df.rds"))
mtbs_all_df <- readRDS(paste0(mtbs_dir, "mtbs_all_df.rds"))

# summarize extent stats by year
mtbs_extents <- NULL
#i=1985
for (i in 1984:2017){
  
  oneyr <- subset(mtbs_all_df, mtbs_all_df$year==i)
  nburncells <- nrow(oneyr)
  mtbs_extent <- cbind(i, nburncells)
  mtbs_extents <- rbind(mtbs_extents, mtbs_extent)
  
}
mtbs_extents <- as.data.frame(mtbs_extents)
names(mtbs_extents) <- c("year", "nburncells")
cells_per_ha <-4
mtbs_extents$extent <- mtbs_extents$nburncells*cells_per_ha
meanhaburnedperyr <- mean(mtbs_extents$extent)
meanhaburnedperyr
sdhaburnedperyr <- sd(mtbs_extents$extent)
sdhaburnedperyr

# fire rotation
totalarea <- 2903305*ha_per_cell #total landscape area in hectares
timeperiod <- max(as.numeric(paste(mtbs_extents$year))) - min(as.numeric(paste(mtbs_extents$year))) + 1
burnedarea <- nrow(mtbs_all_df)*cells_per_ha
rotation <- timeperiod/(burnedarea/totalarea)
rotation 

# Restrict years to 1984-1990
mtbs_extents_84to90 <- subset(mtbs_extents, mtbs_extents$year>=1984 & mtbs_extents$year<=1990) 
mean(mtbs_extents_84to90$extent)
sd(mtbs_extents_84to90$extent)
burnedarea <- sum(mtbs_extents_84to90$extent)
timeperiod <- max(as.numeric(paste(mtbs_extents_84to90$year))) - min(as.numeric(paste(mtbs_extents_84to90$year)))
rotation <- timeperiod/(burnedarea/totalarea)
rotation 

# Restrict years to 1992-2015
mtbs_extents_92to2015 <- subset(mtbs_extents, mtbs_extents$year>=1992 & mtbs_extents$year<=2015) 
mean(mtbs_extents_92to2015$extent)
sd(mtbs_extents_92to2015$extent)
burnedarea <- sum(mtbs_extents_92to2015$extent)
timeperiod <- max(as.numeric(paste(mtbs_extents_92to2015$year))) - min(as.numeric(paste(mtbs_extents_92to2015$year)))
rotation <- timeperiod/(burnedarea/totalarea)
rotation 

 #### Fire Severity ####
#mtbs_all_df_84to90 <- subset(mtbs_all_df, mtbs_all_df$year>=1984 & mtbs_all_df$year<=1990) 
mtbs_all_df_92to2015 <- subset(mtbs_all_df, mtbs_all_df$year>=1992 & mtbs_all_df$year<=2015) 

#mtbs_all_df <- data.frame(mtbs_all_df_84to90)
mtbs_all_df <- data.frame(mtbs_all_df_92to2015)

mtbs_all_df$severity[mtbs_all_df$severity==2] <- "1"
mtbs_all_df$severity[mtbs_all_df$severity==3] <- "2"
mtbs_all_df$severity[mtbs_all_df$severity==4] <- "2"
mtbs_all_df$severity[mtbs_all_df$severity==5] <- "3"
mtbs_all_df$severity[mtbs_all_df$severity==6] <- "3"

summary(as.factor(mtbs_all_df$severity))
hist(as.numeric(mtbs_all_df$severity))
mtbs_all_df$Name <- "Historical"

write.csv(mtbs_all_df, paste0("/Users/shelbyweiss/Documents/FireData_SCRPPLEParameterization/FULL11milha_severity_df_1992_2015.csv"))
