#
#
#
#
#
#
#
#
#
#
# title: Random Forest Modeling
# author: Wesley Rancher
# date: "2025-1-22"

setwd("F:/jetstream/")
dir <- "F:/remoteSensing/data/output/"
#
#
#
#
#
#
#
#
#
#
#
#
library(dplyr)
r_files <- list.files(dir, pattern="-mosaic-", full.names = TRUE, recursive = TRUE)

#Training Dataframe was prepared in (../scripts/grf_ag_carbon.py)
#these are in the same order as the raster bands
training_df <- read.csv("data/RandomForest_TrainingDATA_04082025.csv") %>% select(-PSP) %>% select(-starts_with("tsf_"))
#training_df <- read.csv("data/RandomForest_TrainingDATA_04102025_V2.csv") %>% select(-PSP)#one with Karen Short fire size added
#
#
#
#
#
library(terra)

#files
ref <- rast(paste0(dir, "gee/mosaics/landisReg-FULLmosaic-hist-2015.tif"), lyr=33)#elevation
bounds <- vect("F:/randomForest_Landis/input/shp/FullLandscapeV3_082722.shp")
fieldSites <- vect("F:/randomForest_Landis/input/shp/cafiPlotsCropped.shp")

plot(bounds, main = "Field Site Locations")
plot(fieldSites, add = TRUE, pch = 15, col = "purple")

#sample the landscape (very small proportion of the landscape but should be geographically representative)
regularPoints <- spatSample(ref, size=200000, method = "regular", as.points = TRUE, ext = bounds, na.rm=TRUE)
#randomPoints <- spatSample(ref, size=180000, method = "random", as.points = TRUE, ext = bounds, na.rm=TRUE)

plot(bounds, main = "Regular Points for predictions")
plot(sample(regularPoints, 1000), add = TRUE, pch = 15, col = "deeppink2")#just plotting 1000
#plot(bounds, main = "Random Points for predictions")
#plot(sample(randomPoints, 1000), add = TRUE, pch = 15, col = adjustcolor("blue", alpha.f = 0.85))
#
#
#
#
#
library(future)
#define function to add xy bands
read_raster <- function (file) {
  r <- rast(file)
  
  #remove bands that I managed to duplicate 
  #if (nlyr(r) > 41) r <- terra::subset(r, 1:41)
  
  #filter to our systematic points scheme
  r_vals <- terra::extract(r, regularPoints, xy = TRUE) %>% select(-ID)
  #r_vals$dep2_permafrost <- as.integer(r_vals$dep2_permafrost)#dbl to int
  
  #entire image with no climate
  r_clean <- r[-c(29:32,36:38)]
  blocks(x, n=4)
  #something going on with tsf variables (sub r_vals here for predicting through time)
  #r_vals <- r_vals %>% select(-starts_with("tsf"))
  # r_vals$tsf_0_10 <- as.character(as.logical(r_vals$tsf_0_10))
  # r_vals$tsf_10_20 <- as.character(as.logical(r_vals$tsf_10_20))
  # r_vals$tsf_20_25 <- as.character(as.logical(r_vals$tsf_20_25))
  # r_vals <- r_vals[complete.cases(r_vals), ]
  return(r_clean)
}
#
#
#
#
#
library(data.table)
pfun <- function(object, new_data){
  predicted <- predict(object, new_data)$predictions
}
  
predict_from_model <- function(file, fit, species){
  
  #naming conventions
  spp <- tools::toTitleCase(strsplit(species, " ")[[1]][1])%>%paste0(" ", strsplit(species, " ")[[1]][2])
  scenario <- strsplit(basename(file), "-")[[1]][3]
  yr <- gsub(".*-(\\d+)\\.[a-zA-Z]+$", "\\1", basename(file))
  raster <- read_raster(file)


  #prediction
  print(paste0("Processing ", scenario, " scenario ", "for ", spp, " ", "YEAR: ", yr))
  
  #2000000 pixels (I intentionally left this number lower for making predictions through time)
  predicted <- predict(fit, na.omit(raster)) #%>% 
    # summarise(AvgCarbon_gm2=mean(.pred, na.rm = TRUE),
    #           SdCarbon_gm2=sd(.pred, na.rm = TRUE))
  


  dt <- setDT(as.data.frame(regularPoints))[, AG_C_gm2 :=rep(predicted$.pred, length.out = nrow(regularPoints))]
  #dt$AG_C_gm2 <- predicted
  dt$ClimateScenario <- scenario
  dt$Species <- spp
  dt$Year <- yr
  dt$x <- raster$x
  dt$y <- raster$y
  rm(raster)
  return(dt)
}

make_validation_df <- function(fit, test_data, species, mean_rmse, mean_rsq){
  spp <- tools::toTitleCase(strsplit(species, " ")[[1]][1])%>%paste0(" ", strsplit(species, " ")[[1]][2])
  predicted <- predict(fit, test_data)
  observed <- test_data$Carbon_gm2
  vals <- predicted[[1]]
  df <- rbind(data.frame(Observed=observed, Predicted=vals, Species=spp, RMSE=mean_rmse, RSQ=mean_rsq))
  return(df)
}
#
#
#
#
#
#
#
library(purrr)
library(tidymodels)
library(randomForest)
library(kknn)
library(kernlab)
#library(partykit)
library(vip)

species <- c("black spruce", "white spruce", "resin birch", "quaking aspen")

train_model <- function(df, species, response="Carbon_gm2") {
  print(paste0("running model for: ", species))
  
  #filter to one species
  training_df_one_spp <- training_df %>% filter(Species == species)
  predictors <- setdiff(names(training_df_one_spp), c("Species", "Carbon_gm2"))

  #train and validation set (stratify the predictor dataset based on nature of the data measured in situ)
  split_data <- initial_split(training_df_one_spp %>% select(all_of(response), all_of(predictors)), prop = 0.95, strata = all_of(response))
  
  #training and test set
  train_data <- training(split_data)
  test_data <- testing(split_data)
  formula <- as.formula(paste(response, "~", paste(predictors, collapse = "+")))
  
  #tidy models recipe (swapping in entire dataset)
  recipe <- recipe(formula = formula, data = train_data) #%>% 
    #step_nzv(all_predictors()) %>%
    #step_corr(all_numeric_predictors(), threshold = 0.80, use = "pairwise.complete.obs", method = "pearson")
  
  #define models to run (nts: cite packages)
  #use additional script to set up cforest to be tidymodels friendly
  #source("scripts/define_cforest_model.R")
  m_rf <- rand_forest(mtry = tune(), trees = tune(), min_n = tune()) %>% set_mode("regression") %>% set_engine("randomForest")
  m_ranger <- rand_forest(mtry = tune(), trees = tune(), min_n = tune()) %>% set_mode("regression") %>% set_engine("ranger")
  m_xgb <- boost_tree(mtry = tune(), trees = tune(), min_n = tune()) %>% set_mode("regression") %>% set_engine("xgboost")
  m_kknn <- nearest_neighbor(neighbors=tune(), weight_func = tune(), dist_power = tune()) %>% set_mode("regression") %>% set_engine("kknn")
  m_svm <- svm_rbf(cost= tune(), rbf_sigma = tune()) %>% set_mode("regression") %>% set_engine("kernlab")
  #m_cforest <- cforest() %>% set_mode("regression") %>% set_engine("partykit")
  models <- list(rf=m_rf,ranger=m_ranger)
  
  #params and grids
  #top3_params <- parameters(mtry(), trees(), min_n())
  # I think these two need recipes where the data is normalized
  kknn_params <- parameters(neighbors(range = c(4, 25)), weight_func(values = c("rectangular", "triangular", "biweight", 
                                                                                "triweight", "epanechnikov", "gaussian")),dist_power(range = c(1, 2)))  
  svm_params <- parameters(cost(range = c(-8, 4), trans = scales::log10_trans()), rbf_sigma(range = c(-3, 0), trans = scales::log10_trans()))
    
  #rf, ranger, and xgb use the same params
  rf_grid <- expand.grid(mtry = seq(5, 20, 5), trees = seq(200,500,100), min_n = seq(1, 12, 3))

  kknn_grid <- grid_regular(kknn_params, levels = c(4,3,2))
  svm_grid <- grid_regular(svm_params, levels = 10)
  grids <- list(rf_grid=rf_grid, ranger_grid=rf_grid)
  
  #cross validation
  nfold <- 5
  print(paste0("performing cross validation across ", nfold, " folds"))
  pv_folds <- vfold_cv(train_data, v = nfold, strata = all_of(response))

  #store the results for all models (selecting the best model across the cross-val folds based on RMSE)
  #working outside of the map2 function to generate predictions
  
  #results <- map2(models, grids, function(model, grid) {
  #trying to leave this dynamic
  tuned_results <- workflow() %>% add_recipe(recipe) %>% add_model(m_rf) %>% 
    tune_grid(resamples = pv_folds, grid = rf_grid, metrics = metric_set(rmse, rsq, mae))
  autoplot(tuned_results)
  best_params <- select_best(tuned_results, metric="rmse")
  wflow <- finalize_workflow(workflow() %>% add_recipe(recipe) %>% add_model(m_rf), best_params)
  
  #fit across the resamples
  fit <- fit_resamples(wflow, resamples = pv_folds, metrics = metric_set(rmse, rsq, mae))
  mean_rmse <- collect_metrics(fit) %>% filter(.metric == "rmse") %>% pull(mean)
  mean_rsq <- collect_metrics(fit) %>% filter(.metric == "rsq") %>% pull(mean)
  
  #print(paste0("MODEL: ", model$engine))
  print(paste("Mean RMSE across folds: ", mean_rmse))
  print(paste("Mean R² across folds: ", mean_rsq))
  
  final_fit <- fit(wflow, training_df_one_spp)
  finalized_model <- extract_fit_parsnip(final_fit)
  
  # importance scores (permuation based)
  importance <- vi(finalized_model)
  vip(importance, num_features = 15, geom = "point")
  
  #  return(final_fit)
  #})
  
  #generate predictions (using the regular RF model fit)
  all_predictions <- lapply(unique(r_files), function(file) {predict_from_model(file, final_fit, species)})
  all_predictions_df <- bind_rows(all_predictions)
  
  #generate validation data for plotting
  validation_fit <- fit(wflow, test_data)
  validation_pred_df <- make_validation_df(validation_fit, test_data, species, mean_rmse, mean_rsq)
    
  return(list(train_dfs=train_data, test_dfs=test_data, final_models=finalized_model, 
              all_predictions=all_predictions_df, validation_data=validation_pred_df))
}

#run the models
model_output <- lapply(species, function(spp){train_model(training_df,species=spp)})
predictions <- do.call(rbind, lapply(model_output[1:4], function(x) x$all_predictions))
validations <- do.call(rbind, lapply(model_output[1:4], function(x) x$validation_data))
#write.csv(predictions, "F:/jetstream/data/PredictedCarbon_RasterValues.csv", row.names = FALSE)
#
#
#
#
#
library(ggplot2)
library(viridis)

#plotting validation data vs predictions for each species with a 60/40 split
validation_plot <- ggplot(validations, aes(x = Observed, y = Predicted, color = abs(Predicted-Observed))) +
  #geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE)+
  geom_point(size = 5, alpha = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid") +
  scale_color_viridis(option = "H", direction=-1)+
  facet_grid(~ Species, nrow(1)) +
  theme_bw() +
  scale_x_continuous(limits=c(0,500))+
  scale_y_continuous(limits=c(0,500))+
  labs(x = "Observed", y = "Predicted") +
  theme(legend.position = "none", strip.text = element_text(size = 14), axis.title = element_text(size = 16)) +
  geom_text(aes(x = 0.9 * max(Observed), y = 0.1 * max(Predicted), label = paste("r2:", round(RSQ, 2), "\nRMSE:", round(RMSE, 2))),
            size = 10, hjust = 1, vjust = 0, color = "black")

plot(validation_plot)
#
#
#
#
#

#
#
#
#
#
years <- seq(2006, 2024)
dir <- "F:/Rancher_Sims_Full_Landscape/FinalSimulationsThesis/ClimateData/"
climate_files <- list.files(dir, pattern = paste0("^(tas.*(", paste(years, collapse = "|"), "))\\.tif$"),full.names = TRUE,recursive = TRUE)#2006-2024

#grab predictions from one species
predictions_Bspruce <- predictions %>% filter(ClimateScenario == "hist", Species == "Black spruce")

# Read in the climate data
join_climate_to_predictions <- function(file){
  r <- rast(file)
  
  bname <- basename(file)
  var <- strsplit(bname, "_")[[1]][1]
  scenario <- strsplit(bname, "_")[[1]][5]
  pway <- strsplit(bname, "_")[[1]][6]
  
  yr <- gsub(".*_(\\d{4}).*", "\\1", bname)
  month <- gsub(".*_(\\d{2})_.*", "\\1", bname)
  #clip
  #crs(r) <- crs(ref)
  r_resampled <- project(r, ref)
  clipped <- mask(r_resampled, bounds)
  res(clipped) <- res(r)
  names(clipped) <- paste0(var, month)
  layer_name <-  paste0(var, month)
  
  #extract a subset of the raster values
  r_vals <- terra::extract(clipped, regularPoints, xy = TRUE) %>% select(-ID)
  fname <- paste0("Subset_", bname)
  
  regularPoints[[paste0(layer_name)]] <- r_vals
  smaller_r <- rasterize(regularPoints, clipped, field = layer_name)
  names(smaller_r) <- layer_name
  writeRaster(smaller_r, filename = fname, overwrite = TRUE)
  
  #now we treat the predicted as our "observed"
  joined_vals <- merge(r_vals, predictions_Bspruce, by = c("x","y"))%>%
    select(-paste0(layer_name,".x"))%>%
    rename(paste0(layer_name) == paste0(layer_name, ".y"))
  
  #join
  r_vals$prediction <- predictions_Bspruce$prediction

                       
  return(joined_vals)
}


joined_predictions_new_climate <- lapply(unique(climate_files), function(climate_files) {join_climate_to_predictions(file, final_fit, species)})


#
#
#
#
