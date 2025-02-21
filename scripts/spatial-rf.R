#spatial rf


library(dplyr)
library(spatialRF)


################################################################################
file_dir <- "/Users/wancher/Downloads/"
out_dir <- "/Users/wancher/Documents/thesis/data/output/"
pattern <- "pixel-vals-"

#files
files <- list.files(file_dir, pattern = pattern , full.names = TRUE, recursive = TRUE)
files_to_predict <- list.files(file_dir, pattern = "all-pixels-", full.names = TRUE, recursive = TRUE)

#prediction data stuff
df_for_predicting <- function(file){
  #grab year from filename
  file_name <- basename(file)
  year <- sub('.*-(\\d{4})\\.csv', '\\1', file_name)
  
  #read file and rm index
  df <- read.csv(file)
  df <- df[-c(1)]
  
  #add xy
  coords <- sub('.*\\[([-0-9.]+),([-0-9.]+)\\].*', '\\1,\\2', df$.geo)
  splitxy <- strsplit(coords, ",")
  x <- as.numeric(sapply(splitxy, `[`, 1))
  y <- as.numeric(sapply(splitxy, `[`, 2))
  df$x <- x
  df$y <- y
  df <- df[-c(31)]
  df$year <- year
  return(df)
}
#add year as name for element in list
dfs_predict <- list()
for (file in files_to_predict) {
  df <- df_for_predicting(file)
  dfs_predict[[df$year[1]]] <- df
}

#lapply to df for training data
dfs <- lapply(files, read.csv)
for (i in seq_along(dfs_predict)){
  print(ncol(dfs_predict[[i]]))
}

#remove dfs with less than max columns and rbind
dfs <- dfs[sapply(dfs, function(df) ncol(df) == 37)]
main_df <- do.call(rbind, dfs)

names(main_df)
clean_df <- main_df[-c(1)]#rm sys index

#convert geo column to xy
coords <- sub('.*\\[([-0-9.]+),([-0-9.]+)\\].*', '\\1,\\2', clean_df$.geo)
splitxy <- strsplit(coords, ",")
x <- as.numeric(sapply(splitxy, `[`, 1))
y <- as.numeric(sapply(splitxy, `[`, 2))

#distance matrix
xy <- c(x,y)
clean_df$x <- x
clean_df$y <- y
no_duplicate_xy <- clean_df[!duplicated(clean_df[c("x","y")]),]
locations <- clean_df[c("x", "y")]#need psp?
plot(locations$x, locations$y)
#write.csv(locations, "/Users/wancher/Documents/thesis/data/biomass-input/locations.csv", row.names = FALSE)


################################################################################
response_variables <- c("black cottonwood", "resin birch", "black spruce",
                        "white spruce", "quaking aspen")
training_df <- clean_df %>%
  #mutate(across(-Species, ~ na_if(.x, -9999))) %>%
  filter(Species %in% response_variables) %>%
  rename(year=Year)%>%
  #group_by(site, year, common, across(c(1,2,4:18,20:34))) %>%
  #summarise(total_plot_mass = sum(mass, na.rm = TRUE), .groups = 'drop') %>%
  #left_join(clim_perma, by = c("site","year"))%>%
  
  #convert to g C / m2
  #mutate(ag_carbon = as.numeric(((total_plot_mass / 10000) / 10000) * 0.47)) %>%
  select(-year, -PSP, -SPP_count, -.geo) %>%
  drop_na()

library(terra)

prediction_maps <- list()
summary_dfs <- list()
#response <- response_variables[[2]]
for (response in response_variables) {
    
  training_df_one_spp <- training_df %>% filter(Species == response)
  distance.matrix <- as.matrix(dist(subset(training_df_one_spp, select = -c(x, y))))
  #distance thresholds (same units as distance_matrix)
  thresholds <- quantile(distance.matrix, probs = c(0.25, 0.5, 0.75))
  predictors <- setdiff(names(training_df_one_spp), c("Species", "Biogm2"))
  
  # non spatial 
  model.non.spatial <- spatialRF::rf(
    data = training_df_one_spp,
    dependent.variable.name = "Biogm2",
    predictor.variable.names = predictors,
    distance.matrix = distance.matrix,
    distance.thresholds = thresholds,
    xy = training_df_one_spp[c("x","y")],
    verbose = FALSE
  )
  
  #eval
  model.non.spatial <- spatialRF::rf_evaluate(
    model = model.non.spatial,
    xy = training_df_one_spp[c("x","y")],
    repetitions = 50,         #number of spatial folds
    training.fraction = 0.90, #training data fraction on each fold
    metrics = "r.squared",
    verbose = FALSE
  )
  
  model.spatial <- spatialRF::rf_spatial(
    model = model.non.spatial,
    method = "mem.moran.sequential", #default method
    verbose = FALSE,
  )
  
  print(response)
  spatialRF::print_performance(model.non.spatial)
  spatialRF::print_performance(model.spatial)
  
  ################################################################################
  # predictions
  years <- seq(2000,2011)
  for (year in years){
    
    year_chr <- as.character(year)
    predicted <- stats::predict(
      object = model.non.spatial,
      data = dfs_predict[[year_chr]][, -33],
      type = "response"
    )$predictions
    
    dfs_predict[[year_chr]]$predicted <- predicted #add as column
    dfs_predict[[year_chr]]$Species <- response
    
    summary_df <- dfs_predict[[year_chr]] %>%
      group_by(year, Species)%>%
      summarise(PredAvg = mean(predicted, na.rm = T),
                PredSD = sd(predicted, na.rm = T))
      
    #store in the list
    summary_dfs[[year_chr]] <- summary_df
    
    #for maps
    #coordinates <- vect(dfs_predict[[1]], geom = c("x", "y"))
    #pred_raster <- rasterize(coordinates, terra::rast(), field = "predicted", background = NA)
    #crs(pred_raster) <- "EPSG:3338"
  }
}

################################################################################
#training df plot // do this on cleaned vars
spatialRF::plot_training_df(
  data = training_df_one_spp,
  dependent.variable.name = "Biogm2",
  predictor.variable.names = predictors,
  ncol = 6,
  point.color = viridis::viridis(100, option = "F"),
  line.color = "gray30"
)

#partial dependency
#https://cran.r-project.org/web/packages/pdp/vignettes/pdp-intro.pdf
library(pdp)
pdp::partial(
  model.non.spatial, 
  train = training_df_one_spp, 
  pred.var = "DSM", 
  plot = TRUE, 
  grid.resolution = 200
)

################################################################################
#importance plots
p1 <- spatialRF::plot_importance(
  model.non.spatial, 
  verbose = FALSE) + 
  ggplot2::ggtitle("Non-spatial model") 

p2 <- spatialRF::plot_importance(
  model.spatial,
  verbose = FALSE) + 
  ggplot2::ggtitle("Spatial model")

p1 | p2 


################################################################################
# investigating spatial predictors
library(tigris)
library(sf)
us_states <- tigris::states(cb = TRUE, resolution = "500k")
alaska <- us_states %>% filter(STUSPS == "AK")
alaska_3338 <- st_transform(alaska, crs = 3338)
plot(alaska_3338["STUSPS"]) 

#spatial predictors explained here:
# https://blasbenito.github.io/spatialRF/#fitting-a-spatial-model-with-rf_spatial
spatial.predictors <- spatialRF::get_spatial_predictors(model.spatial)
pr <- data.frame(spatial.predictors, training_df_one_spp[, c("x", "y")])
pr_sf <- st_as_sf(pr, coords = c("x", "y"), crs = 4326)
pr_3338 <- st_transform(pr_sf, crs = 3338)
pr2 <- as.data.frame(st_coordinates(pr_3338))
pr2 <- cbind(pr2, pr_sf[, !(names(pr_sf) %in% c("geometry"))]) %>%
  select(-geometry)

p1 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = alaska_3338, fill = "white") +
  ggplot2::geom_point(
    data = pr2,
    ggplot2::aes(
      x = X,
      y = Y,
      color = spatial_predictor_357.9323_2
    ),
    size = 2.5
  ) +
  ggplot2::scale_color_viridis_c(option = "F") +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "Eigenvalue") +
  #ggplot2::scale_x_continuous(limits = c(-170, -30)) +
  #ggplot2::scale_y_continuous(limits = c(-58, 80))  +
  ggplot2::ggtitle("Variable: spatial_predictor_0_2") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")

p2 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = alaska_3338, fill = "white") +
  ggplot2::geom_point(
    data = pr2,
    ggplot2::aes(
      x = X,
      y = Y,
      color = spatial_predictor_357.9323_5,
    ),
    size = 2.5
  ) +
  ggplot2::scale_color_viridis_c(option = "F") +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "Eigenvalue") +
  #ggplot2::scale_x_continuous(limits = c(-155, 154)) +
  #ggplot2::scale_y_continuous(limits = c(50, 70))  +
  ggplot2::ggtitle("Variable: spatial_predictor_0_5") + 
  ggplot2::theme(legend.position = "bottom") + 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("")

p1 | p2
