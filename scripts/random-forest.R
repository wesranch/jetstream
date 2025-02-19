# title: Random Forest Modeling
# author: Wesley Rancher
# date: "2025-1-22"

################################################################################
#install packages on jetstream
install.packages("dplyr", lib = "~/R/Library")
install.packages("tidymodels", lib = "~/R/Library")
install.packages("randomForest", lib = "~/R/Library")
install.packages("vip", lib = "~/R/Library")
install.packages("ggplot2", lib = "~/R/Library")

#load libraries
library(dplyr)
library(tidymodels)
library(randomForest)
library(vip)
library(ggplot2)
#library(Metrics)
library(viridis)

#environment variables
# library(argparse)
# parser <- ArgumentParser()
# parser$add_argument("--input_file", required = TRUE)
# parser$add_argument("--output_dir", required = TRUE)
# args <- parser$parse_args()
# input_file <- args$input_file
# output_dir <- args$output_dir

################################################################################
# read in training data
dir <- "/Users/wancher/Documents/thesis/"
landsat_topo <- read.csv(paste0(dir, "data/output/pixel-vals-corrected-some.csv"))
#landsat_topo <- read.csv(paste0(dir, "data/output/pixel-vals-indices-topo.csv"))
clim_perma <- read.csv(paste0(dir, "data/output/pixel-vals-climate-perma.csv"))
clim_perma <- clim_perma[c(1,2,3,4,7)]#gfdl only
clim_perma <- clim_perma[c(2,3,5,6,7)]#ccsm only

# data cleaning
response_variables <- c("Alaskan birch", "Black spruce",
                        "White spruce", "Trembling aspen")
training_df <- landsat_topo %>%
  mutate(across(-common, ~ na_if(.x, -9999))) %>%
  filter(common %in% response_variables) %>%
  group_by(site, year, common, across(c(1,2,4:18,20:34))) %>%
  summarise(total_plot_mass = sum(mass, na.rm = TRUE), .groups = 'drop') %>%
  left_join(clim_perma, by = c("site","year"))%>%
  
  #convert to g C / m2
  #mutate(ag_carbon = as.numeric(((total_plot_mass / 10000) / 10000) * 0.47)) %>%
  select(-year, -site, -mass) %>%
  drop_na()

#write.csv(training_df, "data/output/rf-df-2006-2015.csv", row.names = FALSE)


#empty lists
metrics <- list()
autoplots <- list()
rf_workflows <- list()
training_dfs <- list()
test_dfs <- list()

################################################################################
# train models for each species
#response <- response_variables[[4]]
for (response in response_variables) {
  print(paste0("running model for: ", response))
  #filter to one species
  training_df_one_spp <- training_df %>% filter(common == response)
  
  #should the response be carbon or biomass
  predictors <- setdiff(names(training_df_one_spp), c("common", "total_plot_mass"))

  #90/20 split training/validation
  pv_split <- initial_split(training_df_one_spp %>% 
                              select(all_of("total_plot_mass"), all_of(predictors)),
                            prop = 0.90,
                            strata = all_of("total_plot_mass"))
  #train and test set
  pv_train <- training(pv_split)
  pv_test <- testing(pv_split)
  training_dfs[[response]] <- pv_train
  test_dfs[[response]] <- pv_test
  
  formula <- as.formula(paste("total_plot_mass", "~", 
                              paste(predictors, collapse = "+")))

  #recipe
  rf_recipe <- recipe(formula = formula, data = pv_train) %>%
    step_nzv(all_predictors()) %>%
    step_corr(all_numeric_predictors(), threshold = 0.80,
              use = "pairwise.complete.obs", method = "pearson")
  #model specs
  rf_spec_ranger <- rand_forest(mtry = tune(),
                                trees = tune(), min_n = tune()) %>%
    set_mode("regression") %>%
    set_engine("randomForest")
  
  #tuning grid
  ranger_grid <- expand.grid(
    mtry = seq(2, 20, 5),
    trees = 500,
    min_n = seq(3, 9, 3)
  )
  ranger_workflow <- workflow() %>%
    add_recipe(rf_recipe) %>%
    add_model(rf_spec_ranger)
  #cross fold validation and resampling error for Aspen
  pv_folds <- vfold_cv(pv_train, v = 10, strata = all_of("total_plot_mass"))

  #store the results
  ranger_tune_results <- tune_grid(
    ranger_workflow,
    resamples = pv_folds,
    grid = ranger_grid,
    metrics = metric_set(rmse, rsq, mae)
  )

  best_ranger_params <- select_best(ranger_tune_results, metric = "rmse")
  final_ranger_workflow <- finalize_workflow(ranger_workflow,
                                             best_ranger_params)
  rf_workflows[[response]] <- final_ranger_workflow
}

################################################################################
# save models to disk
scenarios <- c("gfdl", "ccsm")
for (scenario in scenarios){
    for (response in response_variables){
    model_fit <- rf_workflows[[response]]
    saveRDS(model_fit, file = paste0("data/output/models/", scenario, "-", response, "-model.rds"))
    }
}

################################################################################
#save variable importance plots
for (response in response_variables){
  pv_train <- training_dfs[[response]]
  final_ranger_workflow <- rf_workflows[[response]]
  final_fit <- fit(final_ranger_workflow, data = pv_train)
  fitted_model <- extract_fit_parsnip(final_fit)$fit #fitted model
  var_imp <- vip::vi(fitted_model) #%>% arrange(desc(Importance))
  
  imp_plot <- ggplot(var_imp, aes(x = reorder(Variable, Importance), y = Importance)) +
    geom_segment(aes(xend = Variable, y = 0, yend = Importance),
                 color = "deepskyblue4", linewidth = 2) +
    geom_point(size = 4,  color = "black") +
    coord_flip() +
    theme_minimal() +
    labs(
      title = paste0("Feature Importance ", response),
      x = "Variable",
      y = "Importance"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      axis.title = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 16, color = "black"),
      axis.title.y = element_text(size = 18, color = "black"),
      axis.title.x = element_text(size = 18, color = "black"),
    )

  ggsave(paste0("writing/figures/vip-", response, ".jpeg"), 
         plot = imp_plot, width = 9, height = 10)
}

################################################################################
# generate predicted vs observed dfs
plot_data <- data.frame()
fits <- list()
for (response in response_variables) {
  pv_train <- training_dfs[[response]]
  pv_test <- test_dfs[[response]]
  final_ranger_workflow <- rf_workflows[[response]]
  final_fit <- fit(final_ranger_workflow, data = pv_train)
  fits[[response]] <- final_fit
  predicted <- predict(final_fit, new_data = pv_test)
  
  #store metrics
  observed <- pv_test$total_plot_mass
  predicted_values <- predicted[[1]]
  rmse_value <- rmse(observed, predicted_values)
  r2_value <- cor(observed, predicted_values)^2
  
  plot_data <- bind_rows(plot_data, data.frame(
    Observed = observed,
    Predicted = predicted_values,
    Response = response,
    RMSE = rmse_value,
    R2 = r2_value
  ))
}

################################################################################
# create predicted vs observed plot
pvo_plot <- ggplot(plot_data, aes(x = Observed, y = Predicted, color = abs(Predicted-Observed))) + # nolint
  geom_smooth(method = "lm", color = "black", linetype = "solid", se = FALSE)+
  geom_point(size = 4, alpha = 0.75) +
  #geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_viridis(option = "plasma", direction = 1)+
  facet_wrap(~Response) +
  theme_minimal() +
  labs(x = "Observed", y = "Predicted") +
  theme(legend.position = "none",
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  geom_text(aes(x = 0.9 * max(Observed), y = 0.1 * max(Predicted),
                label = paste("r2:", round(R2, 2), "\nRMSE:", round(RMSE, 2))),
            size = 10, hjust = 1, vjust = 0, color = "black")
ggsave(paste0("writing/figures/rf-metrics-ccsm.jpeg"), plot = pvo_plot, width = 10, height = 10) #nolint
