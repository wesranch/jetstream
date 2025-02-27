#spatial rf


library(dplyr)
library(spatialRF)
library(tidyr)
library(terra)
library(ggplot2)
library(tidymodels)
library(patchwork)

################################################################################
file_dir <- "/Users/wancher/Downloads/"
out_dir <- "/Users/wancher/Documents/thesis/data/output/"
pattern <- "pixel-vals-3seas"

#files
files <- list.files(file_dir, pattern = pattern , full.names = TRUE, recursive = T)
files_to_predict <- list.files(file_dir, pattern = "all-pixels-", full.names = TRUE, recursive = TRUE)
raster_files <- list.files(out_dir, "dalton-clean-", full.names = T, recursive = T)

#rasters
shp <- vect("/Users/wancher/Documents/thesis/data/landis/dalton_input/Dalton_Landscape.shp")
slope <- rast("/Users/wancher/Documents/thesis/data/output/dalton-slope.tif")
elev <- rast("/Users/wancher/Documents/thesis/data/output/dalton-elev.tif")
aspect <- rast("/Users/wancher/Documents/thesis/data/output/dalton-aspect.tif")

raster_to_df <- function(file){
  file_name <- basename(file)
  year <- sub('.*-(\\d{4})\\.tif', '\\1', file_name)
  
  r <- rast(file)
  r_full <- c(r, elev, slope, aspect)
  df <- as.data.frame(r_full, xy = TRUE)
  df$year <- year
  return(df)
}
raster_dataframes <- lapply(raster_files, raster_to_df)
rasters_clean <- list()
for (df in raster_dataframes) {
  rasters_clean[[df$year[1]]] <- df
}

#comment out if not mapping
dfs_predict <- rasters_clean

#prediction data wrangling
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
for (i in seq_along(dfs)){
  print(ncol(dfs[[i]]))
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
#xy <- c(x,y)
clean_df$x <- x
clean_df$y <- y
#no_duplicate_xy <- clean_df[!duplicated(clean_df[c("x","y")]),]
#locations <- clean_df[c("x", "y")]#need psp?
#plot(locations$x, locations$y)
#write.csv(locations, "/Users/wancher/Documents/thesis/data/biomass-input/locations.csv", row.names = FALSE)


################################################################################
# response_variables <- c("AKbirch", "BCotton",
#                         "WhiteSpruce", "QAspen", "BlackSpruce")
response_variables <- c("resin birch", "black spruce",
                       "white spruce", "quaking aspen")
training_df <- clean_df %>%
  #mutate(across(-Species, ~ na_if(.x, -9999))) %>%
  filter(Species %in% response_variables) %>%
  rename(year=Year)%>%
  #group_by(site, Biogm2, Species, across(c(1:33))) %>%
  #summarise(total_plot_mass = sum(Biogm2, na.rm = TRUE), .groups = 'drop') %>%
  #left_join(clim_perma, by = c("site","year"))%>%
  
  #convert to g C / m2
  #mutate(ag_carbon = as.numeric(((total_plot_mass / 10000) / 10000) * 0.47)) %>%
  select(-year, -PSP, -SPP_count, -.geo) %>%
  drop_na()

################################################################################
# Iterate
models <- list()
spatial_models <- list()
prediction_maps <- list()
summary_dfs <- list()
training_dataframes <- list()
#response <- response_variables[[3]]
for (response in response_variables) {
  
  #filter df to species  
  training_df_one_spp <- training_df %>% 
    filter(Species == response) %>%
    select(-Species)
  
  # tidy models to rm correlated vars
  tm_recipe <- recipe(Biogm2 ~ ., data = training_df_one_spp) %>%
    step_corr(all_numeric_predictors(), threshold = 0.80, 
              use = "pairwise.complete.obs", method = "pearson")
  prepped_recipe <- prep(tm_recipe, training = training_df_one_spp)
  data_corrRM <- bake(prepped_recipe, new_data = training_df_one_spp)
  
  split <- initial_split(data_corrRM, prop = .7)
  train <- training(split)
  test  <- testing(split)
  
  #create clean training df
  train_df <- train %>% 
    #append biomass column
    left_join(training_df_one_spp %>% select(Biogm2,x,y), 
              by = c("x", "y", "Biogm2"), relationship = "many-to-many")
  predictors <- setdiff(names(train_df), "Biogm2")
  
  #create distance matrix
  distance.matrix <- as.matrix(dist(subset(train_df, select = -c(x, y))))
  thresholds <- quantile(distance.matrix, probs = c(0.25, 0.5, 0.75))
  #thresholds <- c(0, 2500, 5000, 7500, 12000)
  # non spatial random forest
  model.non.spatial <- spatialRF::rf(
    data = train_df,
    dependent.variable.name = "Biogm2",
    predictor.variable.names = predictors,
    distance.matrix = distance.matrix,
    distance.thresholds = thresholds,
    xy = train_df[c("x","y")],
    verbose = FALSE) 
  
  #eval
  model.non.spatial <- spatialRF::rf_evaluate(
    model = model.non.spatial,
    xy = train_df[c("x","y")],
    repetitions = 30,         #number of spatial folds
    training.fraction = 0.80, #training data fraction on each fold
    metrics = "r.squared",
    verbose = FALSE)
  
  model.spatial <- spatialRF::rf_spatial(
    model = model.non.spatial,
    method = "mem.moran.sequential", #default method
    verbose = FALSE)
  
  print(response)
  spatialRF::print_performance(model.non.spatial)
  spatialRF::print_performance(model.spatial)
  
  models[[response]] <- model.non.spatial
  spatial_models[[response]] <- model.spatial
  training_dataframes[[response]] <- train_df
  
  ##############################################################################
  #library(ranger)
  #library(SpatialML)
  #geographic random forest package here:
  formula <- as.formula(paste("Biogm2", "~", 
                              paste(predictors, collapse = "+")))
  grf_model <- grf(formula,
                 dframe = train_df, kernel = "adaptive", bw = 50, coords = train_df[c("x","y")], 
                 ntree = 500, forests = TRUE, geo.weighted = TRUE, print.results = TRUE)
  
  ##############################################################################
  # predictions
  prediction_maps_years <- list()
  years <- c(2000,2024)
  for (year in years){
    
    year_chr <- as.character(year)
    #keep columns after baking from above but join to plot locations df or raster
    df_predict <- left_join(dfs_predict[[year_chr]], train_df)%>%
      select(intersect(names(dfs_predict[[year_chr]]), names(train_df)))
  
    
    predicted <- stats::predict(
      object = grf_model$Global.Model,
      data = df_predict,
      type = "response"
    )$predictions
    
    
    df_predict$predicted <- predicted #add as column
    df_predict$Species <- response
    
    summary_df <- df_predict %>%
      #group_by(year, Species)%>%
      summarise(PredAvg = mean(predicted, na.rm = T),
                PredSD = sd(predicted, na.rm = T))
      
    #store in the list
    summary_dfs[[length(summary_dfs) + 1]] <- summary_df
    
    #df to matrix
    #m <- as.matrix(df_predict[, c("x","y","predicted")])
    
    #empty_r <- rast(xmin = min(m[,1]), xmax = max(m[,1]), 
    #          ymin = min(m[,2]), ymax = max(m[,2]), 
    #          resolution = 30) 
    
    #pred_raster <- rasterize(m[, 1:2], empty_r, values = m[,3], background = NA)
    #names(pred_raster) <- "Biogm2_Pred"
    #prediction_maps_years[[year_chr]] <- pred_raster
    #rm(predicted)
    #rm(pred_raster)
    #rm(m, empty_r)
    gc()
  }
  #prediction_maps[[response]] <- prediction_maps_years
  
  rm(model.spatial)
  rm(model.non.spatial)
  gc()
}
for (response in response_variables){
  predictions <- prediction_maps[[response]]
  r2000 <- predictions$`2000`
  r2024 <- predictions$`2024`
  #plot(r2024-r2000)
  response2 <- gsub(" ", "-", response)
  writeRaster(r2000, paste0("/Users/wancher/Documents/thesis/data/output/", response2, "-biomass-2000.tif"), overwrite = TRUE)
  writeRaster(r2024, paste0("/Users/wancher/Documents/thesis/data/output/", response2, "-biomass-2024.tif"), overwrite = TRUE)
  

}
for (i in seq_along(prediction_maps)){
  writeRaster(i, )
}
plotting_df <- do.call(rbind, summary_dfs)
plotting_df$year <- as.numeric(plotting_df$year)
ggplot(plotting_df, aes(y= PredAvg, x = year, color = Species, group = Species))+
  geom_line(linewidth = 2)+
  geom_point(size = 3)+
  scale_color_paletteer_d("lisa::EdwardHopper")+
  #geom_smooth()+
  #facet_wrap(~Species, ncol = 1)+
  theme_minimal()+
  scale_x_continuous(breaks = seq(2000, 2025, by = 5))+
  scale_y_continuous(breaks = seq(0, 1000, by = .1), limits = c(0, 1000))

################################################################################
#training df plot // do this on cleaned vars
spatialRF::plot_training_df(
  data = train_df,
  dependent.variable.name = "Biogm2",
  predictor.variable.names = predictors,
  ncol = 5,
  method = "lm",
  point.color = viridis::viridis(100, option = "plasma"),
  line.color = "black"
)

library(pdp)
#partial dependency
#https://cran.r-project.org/web/packages/pdp/vignettes/pdp-intro.pdf
for (response in response_variables){
  model.non.spatial <- models[[response]]
  train_df <- training_dataframes[[response]]
  depplot <- pdp::partial(
    model.non.spatial, 
    train = train_df, 
    pred.var = "slope", 
    grid.resolution = 200,
    plot = TRUE)
  plot(depplot)
}

################################################################################
#importance plots
for (response in response_variables){
  model.non.spatial <- models[[response]]
  model.spatial <- spatial_models[[response]]
  
  p1 <- spatialRF::plot_importance(
    model.non.spatial, 
    verbose = FALSE) + 
    ggplot2::ggtitle("Non-spatial model")
  
  p2 <- spatialRF::plot_importance(
    model.spatial,
    verbose = FALSE) + 
    ggplot2::ggtitle("Spatial model")
  
  double <- p1 | p2
  plot_this <- double + plot_annotation(title = response)
  print(plot_this)
}

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
