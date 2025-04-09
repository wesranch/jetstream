# cleaning climate rasters
# Wesley Rancher

library(terra)

#files
years <- seq(2010, 2024)
dir <- "F:/remoteSensing/data/output/raster/climate/"
files <- list.files(dir, 
                    pattern = paste0("annual_.*", years, collapse = "|"), 
                    full.names = TRUE)#2006-2024
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
  
  return(clipped)
}

#grab specific variables and scenarios
pr_files <- files[grep("_pr_", basename(files))] 
tas_files <- files[grep("_tas_", basename(files))]
gfdl_pr_files <- pr_files[grep("_GFDL-CM3_", basename(pr_files))]
gfdl_tas_files <- tas_files[grep("_GFDL-CM3_", basename(tas_files))]
ccsm_pr_files <- pr_files[grep("_CCSM4_", basename(pr_files))]
ccsm_tas_files <- tas_files[grep("_CCSM4_", basename(tas_files))]

#loop over years an grab each raster,clean it, and append to random landsat raster
for (year in years){
  #grab file for each year
  gfdl_pr <- gfdl_pr_files[grep(paste0(year), basename(gfdl_pr_files))]
  gfdl_tas <- gfdl_tas_files[grep(paste0(year), basename(gfdl_tas_files))]
  ccsm_pr <- ccsm_pr_files[grep(paste0(year), basename(ccsm_pr_files))]
  ccsm_tas <- ccsm_tas_files[grep(paste0(year), basename(ccsm_tas_files))]
  
  #clean
  gfdl_pr_clean <- clean_raster(gfdl_pr)
  gfdl_tas_clean <- clean_raster(gfdl_tas)
  ccsm_pr_clean <- clean_raster(ccsm_pr)
  ccsm_tas_clean <- clean_raster(ccsm_tas)
  
  #sample (2000-2024 landsat composites)
  yearly_landsat <- landsat_files[grep(paste0(year), basename(landsat_files))]
  landsat_r <- rast(yearly_landsat)
  stacked_gfdl <- c(landsat_r, gfdl_pr_clean, gfdl_tas_clean, perm)
  stacked_ccsm <- c(landsat_r, ccsm_pr_clean, ccsm_tas_clean, perm)
  
  #save scenarios separately
  writeRaster(stacked_gfdl, paste0(mosaics_dir, "landisReg-FULLmosaic-gfdl-", paste0(year), ".tif"), filetype = "GTiff", overwrite = TRUE)
  writeRaster(stacked_ccsm, paste0(mosaics_dir, "landisReg-FULLmosaic-ncar-", paste0(year), ".tif"), filetype = "GTiff", overwrite = TRUE)
  
  rm(yearly_landsat, landsat_r, stacked_gfdl, stacked_ccsm,
     gfdl_pr_clean, gfdl_tas_clean, ccsm_pr_clean, ccsm_tas_clean)
  gc()
}


################################################################################
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
