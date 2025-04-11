###############################################################################
### Compare fire stats from model run to Karen Short and AK Large Fire Dbase###
###############################################################################

#load packages
library(rgdal)
#install.packages('rgeos')
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

cells_per_ha <- 4
totalarea <- 5242453*cells_per_ha

#folder with model run inputs/outputs
model_dir <- "C:/Users/sweiss2/Documents/FULL_Landscape_Inputs_AK_063021_results/"

# read in summarized historic data
fire.stat.hist <- read.csv(paste0(model_dir,"Fire_History_Stats_FULL_useforplotting.csv"))
colnames(fire.stat.hist) <- c("AVG num of fires/yr", "SD num of fires/yr", "AVG ha burned/yr", "SD ha burned/yr",
                         "Max ha burned", "Mean ha/fire", "SD ha/fire", "Fire Rotation", "Name")

hist_size_grouped_AKLFDB_70to90 <- read.csv(paste0(model_dir,"FULL_hist_size_grouped_AKLFDB_1970_to_1990.csv"))
hist_size_grouped_AKLFDB_70to90 <- hist_size_grouped_AKLFDB_70to90[,-1]
hist_size_grouped_AKLFDB_70to90$Name <- "Historic AKLFDB 1970-1990"

hist_size_grouped_AKLFDB_92to15 <- read.csv(paste0(model_dir,"FULL_hist_size_grouped_AKLFDB_1992_to_2015.csv"))
hist_size_grouped_AKLFDB_92to15 <- hist_size_grouped_AKLFDB_92to15[,-1]
hist_size_grouped_AKLFDB_92to15$Name <- "Historic AKLFDB 1992-2015"

hist_size_grouped_KarenShort_92to15 <- read.csv(paste0(model_dir,"FULL_hist_size_grouped_KarenShort_1992_to_2015.csv"))
hist_size_grouped_KarenShort_92to15 <- hist_size_grouped_KarenShort_92to15[,-1]
hist_size_grouped_KarenShort_92to15$Name <- "Historic Karent Short 1992-2015"

severity_84to90 <- read.csv(paste0(model_dir, "FULL_severity_df_1984_1990.csv"))
severity_84to90 <- severity_84to90[,-1]

severity_92to15 <- read.csv(paste0(model_dir, "FULL_severity_df_1992_2015.csv"))
severity_92to15 <- severity_92to15[,-1]


# read in simulation results
fire_dat_original<- read.csv(paste0(model_dir,"Outputs/DGS/scrapple-summary-log.csv"))
fire_dat <- fire_dat_original[,c("SimulationYear","TotalBurnedSitesLightning", "NumberFiresLightning", "TotalBiomassMortalityLightning")]

AvgFireSizes <- (fire_dat[,"TotalBurnedSitesLightning"]* cells_per_ha)/fire_dat[,"NumberFiresLightning"]
TotalFireSizes <- (fire_dat[,"TotalBurnedSitesLightning"]* cells_per_ha)
all.firedat <- cbind(fire_dat, AvgFireSizes, TotalFireSizes)
all.firedat <- all.firedat[complete.cases(all.firedat),]

fire_events <- read.csv(paste0(model_dir,"Outputs/DGS/scrapple-events-log.csv"))
firedat_histogram_df <- fire_events[,c(1,10)]

firedat_histogram_df$TotalSitesBurned <- firedat_histogram_df$TotalSitesBurned * cells_per_ha
firedat_histogram_df$Name <- "Simulated"
names(firedat_histogram_df) <- c("ID", "Hectares", "Name")
firedat_histogram_df$ID <- as.character(firedat_histogram_df$ID)

firedat_histogram_summ <- firedat_histogram_df

breaks <- 100 * c(0:max(firedat_histogram_summ$Hectares))
firedat_histogram_summ$hectbins <- cut(firedat_histogram_summ$Hectares, breaks=breaks, lables=TRUE)
firedat_histogram_grouped <- stats::aggregate(firedat_histogram_summ$Hectares, by=list(firedat_histogram_summ$hectbins), FUN=sum)

firedat_histogram_grouped$percent <- (firedat_histogram_grouped$x/(sum(firedat_histogram_grouped$x)))*100

names(firedat_histogram_grouped) <- c("Bins", "Count", "Percent")
firedat_histogram_grouped$Bins <- as.character(firedat_histogram_grouped$Bins)
tags <- strsplit(firedat_histogram_grouped$Bins, c(",")) %>% unlist() %>% matrix(ncol=2, byrow=TRUE) %>% as.data.frame()
tags <- tags[,1] 
tags <- sub('.', '', tags) %>% as.numeric()
firedat_histogram_grouped$Bins <- tags

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
fire.stat<-data.frame(avg.fires, 0, avg.area, 0, fire.max, fire.mean, 0, rotation, 0, "Simulated")

colnames(fire.stat) <- c("AVG num of fires/yr", "SD num of fires/yr", "AVG ha burned/yr", "SD ha burned/yr",
                         "Max ha burned", "Mean ha/fire", "SD ha/fire", "Fire Rotation", "Years.covered", "Name")

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


#comparison charts 
fire.stat.hist <- fire.stat.hist[,-11]
fire.stat.hist[,10] <- c("Historic AKLFDB 1970-1990", "Historic Karen Short 1992-2015", "Historic AKLFDB 1992-2015", "MTBS 1992-2015")
names(fire.stat.hist) <- c("AVG num of fires/yr", "SD num of fires/yr", "AVG ha burned/yr", "SD ha burned/yr",
                           "Max ha burned", "Mean ha/fire", "SD ha/fire", "Fire Rotation", "Years.covered", "Name")
compare <- rbind(fire.stat, fire.stat.hist)

compare_histogram <- rbind(hist_size_grouped_AKLFDB_70to90, hist_size_grouped_AKLFDB_92to15, hist_size_grouped_KarenShort_92to15, firedat_histogram_grouped)

severity_84to90 <- severity_84to90[,-2]
names(severity_84to90)[1] <- "Severity"
severity_92to15 <- severity_92to15[,-2]
names(severity_92to15)[1] <- "Severity"

names(compare) <- c("FiresPerYr", "FiresPerYr_sd", "HAburnedPerYr", "HAburnedPerYr_sd", "MaxHABurned", "AvgHAPerFire", "SDHAPerFire", "FireRotation", "Years", "Name")


### Plot 1970-1990 comparisons ###
compare_70to90 <- subset(compare, compare$Name=="Simulated"|compare$Name=="Historic AKLFDB 1970-1990")

compare_FiresPerYr_70to90 <- ggplot(data=compare_70to90, aes(fill=Name, y=FiresPerYr, x=Name)) + 
  geom_bar(stat='identity', color="#e9ecef") +
  geom_errorbar(aes(x=Name, y=FiresPerYr, ymin=FiresPerYr, ymax=FiresPerYr + FiresPerYr_sd), width=0, alpha=0.5, size=.8) +
  ggtitle("Average # Fires / Year") +
  scale_fill_manual(values=c("#568956", "#F7E925")) +
  theme_classic() +
  theme(legend.position = "none")
compare_FiresPerYr_70to90

compare_HAburnedPerYr_70to90 <- ggplot(data=compare_70to90, aes(fill=Name, y=HAburnedPerYr, x=Name)) + 
  geom_bar(stat='identity', color="#e9ecef") + 
  geom_errorbar(aes(x=Name, y=HAburnedPerYr, ymin=HAburnedPerYr, ymax=HAburnedPerYr + HAburnedPerYr_sd), width=0, alpha=0.5, size=.8) +
  ggtitle("Average Hectares Burned / Year") +
  scale_fill_manual(values=c("#568956", "#F7E925")) +
  theme_classic() +
  theme(legend.position = "none")
compare_HAburnedPerYr_70to90

compare_MaxHABurned_70to90 <- ggplot(data=compare_70to90, aes(fill=Name, y=MaxHABurned, x=Name)) + 
  geom_bar(stat='identity', color="#e9ecef") + 
  ggtitle("Maximum Fire Size (Ha)") +
  scale_fill_manual(values=c("#568956", "#F7E925")) +
  theme_classic() +
  theme(legend.position = "none")
compare_MaxHABurned_70to90

compare_AvgHAPerFire_70to90 <- ggplot(data=compare_70to90, aes(fill=Name, y=AvgHAPerFire, x=Name)) + 
  geom_bar(stat='identity', color="#e9ecef") + 
  geom_errorbar(aes(x=Name, y=AvgHAPerFire, ymin=AvgHAPerFire, ymax=AvgHAPerFire + SDHAPerFire), width=0.0, alpha=0.5, size=.8) +
  ggtitle("Average Fire Size (Ha)") +
  scale_fill_manual(values=c("#568956", "#F7E925")) +
  theme_classic() +
  theme(legend.position = "none")
compare_AvgHAPerFire_70to90

compare_FireRotation_70to90 <- ggplot(data=compare_70to90, aes(fill=Name, y=FireRotation, x=Name)) +
  geom_bar(stat='identity', color="#e9ecef") + 
  ggtitle("Fire Rotation") +
  scale_fill_manual(values=c("#568956", "#F7E925")) +
  theme_classic() +
  theme(legend.position = "none")
compare_FireRotation_70to90

compare_histogram$Perc_Rnd <-ceiling(compare_histogram$Percent)
#rebuild frequency df based on percentages
compare_histogram_new <- NULL
for (i in 1:nrow(compare_histogram)) {
  #i=88
  compare_histogram_new <- rbind(compare_histogram_new, compare_histogram[rep(i, compare_histogram[i,5]),])
  
}

compare_histogram_new <- compare_histogram_new[,c(1,4)]

compare_histogram_70to90 <- subset(compare_histogram_new, compare_histogram_new$Name=="Simulated"|compare_histogram_new$Name=="Historic AKLFDB 1970-1990")

compare_firesizedistr_70to90 <- ggplot(compare_histogram_70to90, aes(x=Bins, fill=Name)) + 
  geom_histogram(color="#e9ecef", alpha=0.8, position = 'dodge') +
  ggtitle("Fire Size Distribution") +
  scale_fill_manual(values=c("#568956", "#F7E925")) +
  ylab("Percentage of Fires") +
  theme_classic()
compare_firesizedistr_70to90

#fire severity distribution 
compare_severity_84to90 <- rbind(sim_severity_new, severity_84to90)

severity_summ <- data.frame(with(compare_severity_84to90, table(Severity, Name)))
sev_summ_hist <- subset(severity_summ, severity_summ$Name=="Historical")
sev_summ_hist$Percent <- sev_summ_hist$Freq/sum(sev_summ_hist$Freq)*100
sev_summ_sim <- subset(severity_summ, severity_summ$Name=="Simulated")
sev_summ_sim$Percent <- sev_summ_sim$Freq/sum(sev_summ_sim$Freq)*100
severity_summ <- rbind(sev_summ_hist, sev_summ_sim)

compare_severitydistr_84to90 <- ggplot(severity_summ, aes(x=Name, y=Percent, fill=Severity)) +
  geom_bar(color="#e9ecef", position='fill', stat='identity') +
  ggtitle("Fire Severity Distribution") +
  scale_fill_manual(values=c("cornsilk", "orange", "maroon")) +
  ylab("Severity Proportions") +
  theme_classic()
compare_severitydistr_84to90

#plot altogether                                
((compare_FiresPerYr_70to90|compare_HAburnedPerYr_70to90)/
    (compare_MaxHABurned_70to90|compare_AvgHAPerFire_70to90)/
    (compare_FireRotation_70to90|compare_severitydistr_84to90)) /
  compare_firesizedistr_70to90 +
  plot_annotation(title = "Simulation results compared to historic stats 1970-1990")



### Plot 1992-2015 comparisons ###
compare_92to15 <- subset(compare, compare$Name=="Simulated"|compare$Name=="Historic AKLFDB 1992-2015"|
                           compare$Name=="Historic Karen Short 1992-2015")

compare_92to15_wMTBS <- subset(compare, compare$Name=="Simulated"|compare$Name=="Historic AKLFDB 1992-2015"|
                           compare$Name=="Historic Karen Short 1992-2015"|compare$Name=="MTBS 1992-2015")


compare_FiresPerYr_92to15 <- ggplot(data=compare_92to15, aes(fill=Name, y=FiresPerYr, x=Name)) + 
  geom_bar(stat='identity', color="#e9ecef") +
  geom_errorbar(aes(x=Name, y=FiresPerYr, ymin=FiresPerYr, ymax=FiresPerYr + FiresPerYr_sd), width=0, alpha=0.5, size=.8) +
  ggtitle("Average # Fires / Year") +
  scale_fill_manual(values=c("#568956", "#67D567", "#F7E925")) +
  theme_classic() +
  theme(legend.position = "none")
compare_FiresPerYr_92to15

compare_HAburnedPerYr_92to15 <- ggplot(data=compare_92to15_wMTBS, aes(fill=Name, y=HAburnedPerYr, x=Name)) + 
  geom_bar(stat='identity', color="#e9ecef") + 
  geom_errorbar(aes(x=Name, y=HAburnedPerYr, ymin=HAburnedPerYr, ymax=HAburnedPerYr + HAburnedPerYr_sd), width=0, alpha=0.5, size=.8) +
  ggtitle("Average Hectares Burned / Year") +
  scale_fill_manual(values=c("#568956", "#67D567", "#A7F5A7", "#F7E925")) +
  theme_classic() +
  theme(legend.position = "none")
compare_HAburnedPerYr_92to15

compare_MaxHABurned_92to15 <- ggplot(data=compare_92to15, aes(fill=Name, y=MaxHABurned, x=Name)) + 
  geom_bar(stat='identity', color="#e9ecef") + 
  ggtitle("Maximum Fire Size (Ha)") +
  scale_fill_manual(values=c("#568956", "#67D567", "#F7E925")) +
  theme_classic() +
  theme(legend.position = "none")
compare_MaxHABurned_92to15

compare_AvgHAPerFire_92to15 <- ggplot(data=compare_92to15, aes(fill=Name, y=AvgHAPerFire, x=Name)) + 
  geom_bar(stat='identity', color="#e9ecef") + 
  geom_errorbar(aes(x=Name, y=AvgHAPerFire, ymin=AvgHAPerFire, ymax=AvgHAPerFire + SDHAPerFire), width=0.0, alpha=0.5, size=.8) +
  ggtitle("Average Fire Size (Ha)") +
  scale_fill_manual(values=c("#568956", "#67D567", "#F7E925")) +
  theme_classic() +
  theme(legend.position = "none")
compare_AvgHAPerFire_92to15

compare_FireRotation_92to15 <- ggplot(data=compare_92to15_wMTBS, aes(fill=Name, y=FireRotation, x=Name)) +
  geom_bar(stat='identity', color="#e9ecef") + 
  ggtitle("Fire Rotation") +
  scale_fill_manual(values=c("#568956", "#67D567", "#A7F5A7", "#F7E925")) +
  theme_classic() +
  theme(legend.position = "none")
compare_FireRotation_92to15

compare_histogram$Perc_Rnd <-ceiling(compare_histogram$Percent)
#rebuild frequency df based on percentages
compare_histogram_new <- NULL
for (i in 1:nrow(compare_histogram)) {
  #i=88
  compare_histogram_new <- rbind(compare_histogram_new, compare_histogram[rep(i, compare_histogram[i,5]),])
  
}

compare_histogram_new <- compare_histogram_new[,c(1,4)]

compare_histogram_92to15 <- subset(compare_histogram_new, compare_histogram_new$Name=="Simulated"|compare_histogram_new$Name=="Historic AKLFDB 1992-2015"|
                                     compare_histogram_new$Name=="Historic Karent Short 1992-2015")

compare_firesizedistr_92to15 <- ggplot(compare_histogram_92to15, aes(x=Bins, fill=Name)) + 
  geom_histogram(color="#e9ecef", alpha=0.8, position = 'dodge') +
  ggtitle("Fire Size Distribution") +
  scale_fill_manual(values=c("#568956", "#67D567", "#F7E925")) +
  ylab("Percentage of Fires") +
  theme_classic()
compare_firesizedistr_92to15

#fire severity distribution 
compare_severity_92to15 <- rbind(sim_severity_new, severity_92to15)

severity_summ <- data.frame(with(compare_severity_92to15, table(Severity, Name)))
sev_summ_hist <- subset(severity_summ, severity_summ$Name=="Historical")
sev_summ_hist$Percent <- sev_summ_hist$Freq/sum(sev_summ_hist$Freq)*100
sev_summ_sim <- subset(severity_summ, severity_summ$Name=="Simulated")
sev_summ_sim$Percent <- sev_summ_sim$Freq/sum(sev_summ_sim$Freq)*100
severity_summ <- rbind(sev_summ_hist, sev_summ_sim)

compare_severitydistr_92to15 <- ggplot(severity_summ, aes(x=Name, y=Percent, fill=Severity)) +
  geom_bar(color="#e9ecef", position='fill', stat='identity') +
  ggtitle("Fire Severity Distribution") +
  scale_fill_manual(values=c("cornsilk", "orange", "maroon")) +
  ylab("Severity Proportions") +
  theme_classic()
compare_severitydistr_92to15

#plot altogether                                
((compare_FiresPerYr_92to15|compare_HAburnedPerYr_92to15)/
    (compare_MaxHABurned_92to15|compare_AvgHAPerFire_92to15)/
    (compare_FireRotation_92to15|compare_severitydistr_92to15)) /
  compare_firesizedistr_92to15 +
  plot_annotation(title = "Simulation results compared to historic stats 1992-2015")

