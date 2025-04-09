# mosaic gee outputs
# Wesley Rancher


library(terra)
library(dplyr)
library(foreach)
library(doParallel)
library(future)

# directories
file_dir <- "F:/remoteSensing/data/output/gee/tiles/"
output_dir <- "F:/remoteSensing/data/output/gee/mosaics/"
#setwd(file_dir)

# read in the rasters 
list_of_files <- list.files(file_dir, pattern = "landis-tile-.*\\.tif$", full.names = TRUE)


# iterate over a sequence of years and pull out files specific to the year in the sequence
years <- c(2008,2019,2023)
years_string <- as.character(years)

#parallel
#num_cores <- detectCores() - 2
#cl <- makeCluster(num_cores)
#registerDoParallel(cl)

#foreach(year = years_string, .packages = c("terra", "dplyr")) %dopar% {
plan(multisession, workers = availableCores() - 1)
for (i in seq_along(years_string)) {
  
  #pull out unique year
  year <- years_string[[3]]
  files_one_year <- list_of_files[grepl(year, list_of_files)]
  #files_one_year <- list_of_files
  
  list_of_rasters <- lapply(files_one_year, function(file) {
    r <- rast(file)
    r
  })
  for (i in seq_along(list_of_rasters)){
    plot(list_of_rasters[[i]]$slope)
  }
  # get band names for retention
  band_names <- names(list_of_rasters[[5]])
  #flattened_band_names <- unlist(band_names)
  
  # turn list in sprc and mosaic
  rsrc <- terra::sprc(list_of_rasters)
  m <- mosaic(rsrc)#2019 1:29pm 4/2
  #coll_of_rasters <- sprc(list_of_rasters)
  print(paste0(year, " start: ", Sys.time()))
  #mosaiced_raster <- merge(coll_of_rasters) 
  names(m) <- band_names
  
  #fire to categorical
  m$tsf_0_10 <- (m >= 0) & (m < 10)
  m$tsf_10_20 <- (m >= 10) & (m < 20)
  m$tsf_20_25 <- (m >= 20) & (m < 25)
  
  m2 <- subset(m, 'tsf', negate=TRUE)
  
  #save it
  output_filename <- paste0(output_dir, "landisReg-mosaic-", year, ".tif")
  print(output_filename)
  writeRaster(m2, filename = output_filename, filetype = "GTiff", overwrite = TRUE)
  print(paste0(year, " finish: ", Sys.time()))
  rm(files_one_year, list_of_rasters, coll_of_rasters, m, m2)
  gc()
}
stopCluster(cl)

