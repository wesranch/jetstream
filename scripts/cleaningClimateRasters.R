# cleaning climate rasters
# Wesley Rancher

library(terra)

#files
dir <- "F:/remoteSensing/data/output/raster/climate/"
years <- seq(2025, 2050)
files <- list.files(dir, 
                    pattern = paste0("annual_.*", years, collapse = "|"), 
                    full.names = TRUE)#2025-2050
bounds <- vect("F:/randomForest_Landis/input/shp/FullLandscapeV3_082722.shp")
ref <- rast("F:/remoteSensing/data/output/gee/mosaics/landisReg-mosaic-2015.tif")
perm <- rast("F:/remoteSensing/data/output/raster/climate/permafrost/mu_permafrost_0_100_2.img")
names(perm) <- "dep2_permafrost"

#function to clip future climate rasters and append randomly sampled landsat composites 2000-2024
clean_raster <- function(file){
  r <- rast(file)
  bname <- basename(file)
  var <- strsplit(bname, "_")[[1]][2]
  scenario <- strsplit(bname, "_")[[1]][3]
  yr <- gsub(".*_(\\d{4}).*", "\\1", bname)
  
  #clip
  #crs(r) <- crs(ref)
  r_resampled <- project(r, ref)
  clipped <- mask(r_resampled, bounds)
  names(clipped) <- var
  
  #randomly sample landsat image and append
  random_file <- sample(landsat_files, 1, replace = TRUE)
  random_r <- rast(random_file)
  stacked <- c(random_r, clipped, perm)
  
  return(stacked)
}
clipped_climate_rasters <- lapply(files, clean_raster)


#function to add historic climate rasters to landsat composites
pr_files_historic <- list.files(dir, pattern = "pr_.*1971_2000.tif", full.names = TRUE, recursive = TRUE)
tas_files_historic <- list.files(dir, pattern = "tas_.*1971_2000.tif", full.names = TRUE, recursive = TRUE)

#append historic to landsat comps 2000-2024
clean_hist_climate <- function(files){
  stack <- rast(files)
  bname <- basename(files[[1]])
  var <- strsplit(bname, "_")[[1]][1]
  annual_mean_normal <- mean(stack)
  r_resampled <- resample(annual_mean_normal, ref, method = "near")
  #r_resampled <- resample(perm, ref, method = "near")
  masked <- mask(r_resampled, bounds)
  #perm <- masked
  names(masked) <- var
  return(masked)
}

pr_historic_raster <- clean_hist_climate(pr_files_historic)
tas_historic_raster <- clean_hist_climate(tas_files_historic)

mosaics_dir <- "F:/remoteSensing/data/output/gee/mosaics/"
landsat_files <- list.files(mosaics_dir, pattern = ".tif", full.names = TRUE)

actually_append_climate <- function(file){
  r <- rast(file)
  bname <- basename(file)
  yr <- gsub(".*-(\\d{4}).*", "\\1", bname)
  stacked <- c(r, pr_historic_raster, tas_historic_raster, perm)
  writeRaster(stacked, paste0(mosaics_dir, "landisReg-FULLmosaic-", yr, ".tif"), filetype = "GTiff", overwrite = TRUE)
}
lapply(landsat_files, actually_append_climate)
