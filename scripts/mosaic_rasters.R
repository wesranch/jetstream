# mosaic gee outputs
# Wesley Rancher

install.packages("terra")
install.packages("dplyr")
library(terra)
library(dplyr)

# dir and outdir
file_dir <- "/home/wrancher/rs-data/tiles/"
output_dir <- "/home/wrancher/rs-data/mosaics/"
setwd(file_dir)

# read in the rasters 
list_of_files <- list.files(file_dir, pattern = "FullComp_.*\\.tif$", full.names = TRUE)


# iterate over a sequence of years and pull out files specific to the year in the sequence
years <- seq(2000, 2023)
years_string <- as.character(years)
for (i in seq_along(years_string)) {
  
  #pull out unique year
  year <- years_string[[i]]
  files_one_year <- list_of_files[grepl(year, list_of_files)]

  
  list_of_rasters <- lapply(files_one_year, function(file) {
    r <- rast(file)
    r
  })
  
  # get band names for retention
  band_names <- names(list_of_rasters[[1]])
  #flattened_band_names <- unlist(band_names)
  
  # turn list in sprc and mosaic
  coll_of_rasters <- sprc(list_of_rasters)
  print(paste0(year, " start: ", Sys.time()))
  mosaiced_raster <- merge(coll_of_rasters) 
  print(paste0(year, " finish: ", Sys.time()))
  names(mosaiced_raster) <- band_names
  
  #save it
  output_filename <- paste0(output_dir, "comp-mosaic-", year, ".tif")
  print(output_filename)
  writeRaster(mosaiced_raster, filename = output_filename, filetype = "GTiff", overwrite = TRUE)
  rm(list_of_rasters, coll_of_rasters, mosaiced_raster)
  gc()
}


