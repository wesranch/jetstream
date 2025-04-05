#%% libraries
import os
import glob
import pandas as pd
import numpy as np
import dask.array as da
import rioxarray
import rasterio as rast
import xarray as xr
from affine import Affine
from rasterio import features
from pyproj import Proj, Transformer
import re
from rioxarray.merge import merge_arrays
from itertools import product
from scipy.stats import gaussian_kde
import seaborn as sns
from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import root_mean_squared_error, r2_score
from PyGRF import PyGRF #geographic rf

import pickle
from sklearn.inspection import PartialDependenceDisplay
import matplotlib.pyplot as plt
import seaborn as sns

"""

Author: @wesranch
Date: 2025-04-04

This script combines processed Landsat data (2000-2024), 
historic climate data from PRISM (normals 1971-2000), 
and topographic and permafrost data to train and fit random forest models.

Predictions are made using Landsat composites plus historic climate rasters to have trends under "historic climate"
+ predictions using Landsat composites plus future climate from CCSM4 (RCP8.5 CMIP5 "intermediate climate")
+ predictions using Landsat composites plus future climate from GFDL (RCP8.5 CMIP5 "extreme climate")

"""

#files (raster values at plot locations)
path = "F:/remoteSensing/data/output/csv/Sampled_Alaska/pixel-vals-fire*.csv"
files = glob.glob(path)
climate_df = pd.read_csv("F:/remoteSensing/data/output/csv/pixel-vals-hist-clim.csv")
climate_df['pr'] = climate_df['hist_pr_1971_2000']
climate_df['tas'] = climate_df['hist_tas_1971_2000']
climate_df.drop(['hist_tas_1971_2000', 'hist_pr_1971_2000'], axis=1, inplace=True)
permafrost_df = pd.read_csv("F:/remoteSensing/data/output/csv/pixel-vals-perm.csv")
y_variable = 'Carbon_gm2'
response_variables = ['black spruce', 'white spruce', 'quaking aspen', 'resin birch']

#function to save tiff
def save_to_tiff(data, output_path, raster_template):
    data = data.compute()
    transform = raster_template.transform
    crs = raster_template.crs
    height, width = data.shape
    data.rio.to_raster(output_path, transform=transform, crs=crs, dtype='float32', 
                       tiled=True, windowed=True)

#patterns
mosaics_historic = "F:/remoteSensing/data/output/gee/mosaics/landisReg-FULLmosaic-hist-*.tif"
mosaics_intermediate = "F:/remoteSensing/data/output/gee/mosaics/landisReg-FULLmosaic-ncar-*.tif"
mosaics_extreme = "F:/remoteSensing/data/output/gee/mosaics/landisReg-FULLmosaic-gfdl-*.tif"

#list files
raster_files_historic = glob.glob(mosaics_historic)
raster_files_intermediate = glob.glob(mosaics_intermediate)
raster_files_extreme = glob.glob(mosaics_extreme)

#%% read in rasters to predict from
def read_raster(file):
    """
    Reads in each band of landsat composites + climate rasters and adds xy as bands
    These are needed as x and y are used as predictors in the models
    """
    rds = rioxarray.open_rasterio(file, lock=False, cache=False, chunks={'x': 256, 'y': 256})
    #rds = rioxarray.open_rasterio(file, lock=False, cache=False)
    
    #metadata
    bands = len(rds.long_name)
    transform = rds.rio.transform()
    height, width = rds.shape[1], rds.shape[2]
    rds.name = "bands"
    
    #band names
    for i in range(bands):
        print(rds.long_name[i])
    
    #generate coords
    x_coords = xr.DataArray(rds.x, name="x")
    y_coords = xr.DataArray(rds.y, name="y")
  
    #convert to 2d array
    x_coords_2d = x_coords.expand_dims(dim={'y': height})
    y_coords_2d = y_coords.expand_dims(dim={'x': width})
    
    #3d array for stacking
    x_coords_3d = xr.DataArray(x_coords_2d.expand_dims(dim={'band': [0]}), name='x')
    y_coords_3d = xr.DataArray(y_coords_2d.expand_dims(dim={'band': [0]}), name='y')
    x_coords_3d = x_coords_3d.transpose('band', 'y', 'x')
    y_coords_3d = y_coords_3d.transpose('band', 'y', 'x')
    
    #stack bands
    stacked_array = rioxarray.merge.merge_arrays([rds, x_coords_3d, y_coords_3d])
    #stacked_array = xr.merge([rds, x_coords_expanded, y_coords_expanded])

    #export order of band names to make sure columns match
    band_names = list(rds.long_name)
    band_names.extend(["x", "y"])
    
    print(f"all band names: {band_names}")
    
    return stacked_array, band_names, transform, width, height




#%% process, clean, join
def clean_dataframe(file):
    df = pd.read_csv(file)
    df[f'{y_variable}'] = df['Biogm2']*0.47#units of carbon
    df[['x','y']] = df['.geo'].str.extract(r'\[([-0-9.]+),([-0-9.]+)\]').astype(float)

    #turn time since fire categorical
    df['tsf_0_10'] = (df['tsf'] >= 0) & (df['tsf'] < 10)
    df['tsf_10_20'] = (df['tsf'] >= 10) & (df['tsf'] < 20)
    df['tsf_20_25'] = (df['tsf'] >= 20)
    
    df.drop(['system:index','Biogm2','SPP_count','.geo', 'severity_background',
             'severity_non_mapping','severity_increased_greenness', 'tsf'], 
             axis=1, inplace = True)

    df.dropna(inplace=True)
    print(df.shape[1])
    return df

dfs = list(map(clean_dataframe,files))#lapply equivalent

#remove dfs that have less columns than they should have
dfs = [df for df in dfs if df.shape[1] == 44]
landsat_df = pd.concat(dfs,axis=0)#rbind
landsat_df = landsat_df.reset_index(drop=True)

#climate and perm
clim_perm_df = pd.concat([climate_df, permafrost_df],axis=1)
clim_perm_df = clim_perm_df.loc[:, ~clim_perm_df.columns.duplicated()].copy()
clim_perm_df = clim_perm_df.reset_index(drop=True)

#join
full_df = pd.concat([landsat_df, clim_perm_df], axis=1).dropna(inplace=False)
full_df = full_df.loc[:, ~full_df.columns.duplicated()].copy().drop(['year','Year'], axis=1)
print(full_df.columns)

#%% #standardize and remove correlated vars
def standardize(data, stats):
    std_replaced = stats['std'].replace(0, np.nan)  # Replace 0 std with NaN to avoid division error
    standardized_data = (data - stats['mean']) / std_replaced
    standardized_data = standardized_data.fillna(0)  # Replace NaNs (from zero std) with 0
    return standardized_data

# correlation
def remove_correlated(df, stats):
    corr_matrix = df.corr().abs()
    upperT = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
    drop = [col for col in upperT.columns if any(upperT[col]>0.80)]
    df_filtered = df.drop(columns=drop)
    df_standardized = standardize(df_filtered, stats)
    #plot
    sns.heatmap(df_filtered, annot=False, cmap='coolwarm', fmt=".2f")
    return df_standardized

#%% #hyperparam tuning grid
def hyperparam_grid(predictors, bw, lw, x_train_df, y_train_df, 
                    x_test_df, y_test_df, coords, coords_test):
    """
    Applies grid search combinations for number of trees and variables tried at each split.
    Used in tuning the model.
    """
    param_grid = {
        'n_estimators': np.arange(40, 100, 5),
        'max_features': np.arange(1, min(6, len(predictors)) + 1, 2)
    }
    param_combinations = list(product(param_grid['n_estimators'], param_grid['max_features']))
    best_rmse = float('inf')
    best_params = None

    for n_estimators, max_features in param_combinations:
        pygrf = PyGRF.PyGRFBuilder(band_width=bw, train_weighted=True, predict_weighted=True, 
                                   bootstrap=False, resampled=True,
                                   n_estimators=n_estimators, max_features=max_features)

        pygrf.fit(x_train_df, y_train_df, coords)
        predict_combined, _, _ = pygrf.predict(x_test_df, coords_test, local_weight=lw)

        rmse = root_mean_squared_error(y_test_df, predict_combined)  # RMSE
        r_squared = r2_score(y_test_df, predict_combined)

        print(f"Params: n_estimators={n_estimators}, max_features={max_features} --> RMSE: {rmse:.4f}, R²: {r_squared:.4f}")

        # Update best parameters if current RMSE is lower
        if rmse < best_rmse:
            best_rmse = rmse
            best_params = {'n_estimators': n_estimators, 'max_features': max_features}

    print(f"Best Params: {best_params}, Best RMSE: {best_rmse:.4f}")
    
    return best_params

#cross validation
def kfoldit(df, predictors, bw, best_params, x_train_df, y_train_df, 
            x_test_df, y_test_df, coords, coords_test, lw, response):
    """
    Cross validation to find optimal tuning params account for potential model variance
    Used to determine optimal params for tuning final model for each species
    """
    #cross validation
    folds = KFold(n_splits=10, shuffle=True, random_state=63)
    rmses = []
    for i, (train_index, test_index) in enumerate(folds.split(df)):
        print(response, ":")
        print(f'fold {i+1}')

        X_train_fold, X_test_fold = df.iloc[train_index], df.iloc[test_index]
        y_train_fold, y_test_fold = X_train_fold[y_variable], X_test_fold[y_variable]
        xy_coords_train_fold = X_train_fold[['x', 'y']]
        xy_coords_test_fold = X_test_fold[['x', 'y']]
        X_train_IV = X_train_fold[predictors]
        X_test_IV = X_test_fold[predictors]

        #standardize independent vars
        training_stat = X_train_IV.describe().transpose()
        X_train_scaled = standardize(X_train_IV, training_stat)
        X_test_scaled = standardize(X_test_IV, training_stat)
        
        #build and fit
        pygrf = PyGRF.PyGRFBuilder(band_width=bw, train_weighted=True, predict_weighted=True, 
                                   bootstrap=False, resampled=True, random_state=42,
                                   **best_params)
        #pygrf.fit(x_train_df, y_train_df, coords)
        predict_combined, _, _ = pygrf.predict(x_test_df, coords_test, local_weight=lw)

        #metrics
        rmse = root_mean_squared_error(y_test_df, predict_combined)  # RMSE
        rmses.append(rmse)#take an average of this
        
    return rmses

#%% #function to predict and store metrics from stored models
def predict_and_summarize(bw, best_params, lw, xtrain, ytrain, test_data, ytest, xy_coords, response):
    """
    Main function for generating results. Uses optimal parameters from grid search
    Generates prediction maps, predicted vs observed plots, and variable importance graphs for both local and global models from PyGRF
    """
    
    #empties
    y_predict = []
    y_true = []
    local_vi_df = pd.DataFrame()
    global_vi_df = pd.DataFrame()
    metrics_dfs = []
    models = []
    combined_predictions = []
    global_predictions = []
    local_predictions = []
    scaled_training_df = pd.DataFrame()
    
    #build and fit
    pygrf = PyGRF.PyGRFBuilder(band_width=bw, train_weighted=True, predict_weighted=True, 
                               bootstrap=False, resampled=True,
                               **best_params)
    #pygrf.fit(xtrain, ytrain, xy_coords)
    
    #predict_combined, predict_global, predict_local = pygrf.predict(test_data, xy_coords, local_weight=lw)
    
    #predict on raster
    def predict_and_save(file, response, build, xtrain, ytrain, xy_coords, lw):
        """
        Reads in raster files and applies a rasterize function for adding x and y bands
        Generates a fit from the model and predicts on the raster and stores the output
        """
        #read in the raster
        bname = os.path.basename(file)
        splitname_right = bname.rsplit('-', 1)
        splitname = bname.rsplit('-')
        yr = splitname_right[1].rsplit(".")[0]
        scenario = splitname[2]
        spp = response.replace(' ', '').lower()
        
        print(f"reading raster data for year {yr} for {scenario} climate")
        stacked_array, band_names, transform, width, height = read_raster(file)#3:01p Friday 4/4
        
        #fit and predict
        xtrain = xtrain[band_names]#ensure the order of the random forest data columns match the raster layer order
        build.fit(xtrain, ytrain, xy_coords)
        print(f"Making prediction map for {response} year {yr} using {scenario} climate")
        predict_combined, predict_global, predict_local = build.predict(stacked_array, xy_coords, local_weight=lw)
        
        output_path = f"{output_dir}/Carbon_{spp}_{scenario}_{yr}.tif" #gm2
        save_to_tiff(predict_combined.compute(), output_path, raster_template)
        print(f"Prediction map saved for {response} year {yr} with {scenario} climate")
        
        stacked_array = None
        predict_combined = None
    
    #save predictions
    output_dir = "F:/remoteSensing/data/output/raster/"
    raster_template = rast.open("F:/remoteSensing/data/output/gee/mosaics/landisReg-FULLmosaic-hist-2015.tif")
    
    predictions_historic = list(map(lambda file: predict_and_save(file, response, pygrf, xtrain, ytrain, xy_coords, 
                                                         lw, output_dir, raster_template), raster_files_historic))
    predictions_intermediate = list(map(lambda file: predict_and_save(file, response, pygrf, xtrain, ytrain, xy_coords, 
                                                         lw, output_dir, raster_template), raster_files_intermediate))
    predictions_extreme = list(map(lambda file: predict_and_save(file, response, pygrf, xtrain, ytrain, xy_coords, 
                                                         lw, output_dir, raster_template), raster_files_extreme))
    
    #local variable importance
    local_model_vip = pd.DataFrame(pygrf.get_local_feature_importance())
    #local_model_vip['Fold'] = i + 1
    local_model_vip['Species'] = response

    #global rf model var importance
    global_model_features = pygrf.global_model.feature_names_in_
    global_model_importance_scores = pygrf.global_model.feature_importances_
    global_vip_df = pd.DataFrame({
        'Feature': global_model_features,
        'Importance': global_model_importance_scores,
        #'Fold': i+1,
        'Species': response
    })  
    
    global_vip_df = global_vip_df.sort_values(by='Importance',ascending=False).reset_index(drop=True)
    local_vi_df = pd.concat([local_vi_df, local_model_vip], ignore_index=True)
    global_vi_df = pd.concat([global_vi_df, global_vip_df], ignore_index=True)

    #global feature importance
    numeric_df = global_vip_df.loc[global_vip_df['Species']==f'{response}'].drop(columns=['Species'], errors='ignore')
    indices= numeric_df.sort_values(by='Importance',ascending=False).head(25)
    #local feature importance
    numeric_df_local = local_model_vip.loc[local_model_vip['Species']==f'{response}'].drop(columns=['Species', 'model_index'], errors='ignore')
    indices_local= numeric_df_local.mean().sort_values(ascending=False).head(25)
    #figure
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 12))  # Tall and narrow
    colors = plt.cm.cividis_r(indices['Importance'] / max(indices['Importance']))
    colors_local = plt.cm.cividis_r(indices_local.values / max(indices_local.values))
    # global on left
    axes[0].barh(indices['Feature'], indices['Importance'], color=colors)
    axes[0].set_xlabel("Global Importance Score")
    axes[0].set_ylabel("Features")
    #axes[0].set_title(f"LANDIS-II Feature Importance - {response}")
    axes[0].invert_yaxis()
    #local on the right
    axes[1].barh(indices_local.index, indices_local.values, color=colors_local)
    axes[1].set_xlabel("Local Importance Score")
    axes[1].set_ylabel("Features")
    #axes[1].set_title(f"Local Importance - {response}")
    axes[1].invert_yaxis()
    fig.suptitle(f"RS Feature Importance - {response}", fontsize=16, ha='center')
    plt.tight_layout()
    plt.show()
    
    #store predictions
    y_predict.extend(predict_combined)
    y_true.extend(ytest.tolist())

    rmse = root_mean_squared_error(ytest, predict_combined)
    rSquared = r2_score(ytest, predict_combined)
    
    #store metrics
    metrics_dict = {'RMSE': float(rmse), 'R2': float(rSquared), 'Species': response}
    metrics_dfs.append(metrics_dict)
    
    #print("rmse: " + str(round(rmse, 3)), "r2: " + str(round(rSquared, 3)))

    #store models and predictions
    models.append(pygrf)
    combined_predictions.extend(predict_combined)
    global_predictions.extend(predict_global)
    local_predictions.extend(predict_local)
    
    #predicted vs observed plot
    fig, ax = plt.subplots(figsize=(8, 6))
    xy = np.vstack([ytest, predict_combined])
    kde = gaussian_kde(xy)
    density = kde(xy)
    norm = plt.Normalize(vmin=np.min(density), vmax=np.max(density))
    colors = plt.cm.viridis(norm(density))
    ax.scatter(ytest, predict_combined, alpha=0.7, color=colors, edgecolors='b', s=50)
    ax.plot([min(ytest), max(ytest)], [min(ytest), max(ytest)], color='red', linestyle='--')
    ax.set_xlabel("Observed Values")
    ax.set_ylabel("Predicted Values")
    ax.set_title(f"Predicted vs Observed for {response}")
    ax.text(0.05, 0.95, f'R² = {rSquared:.4f}', transform=ax.transAxes, fontsize=12, verticalalignment='top', color='black')
    plt.tight_layout()
    plt.show()
    
    #clean
    pygrf = None
    predict_combined = None
    predict_global = None
    predict_local = None


    final_metrics_df = pd.DataFrame(metrics_dfs)

    return models, combined_predictions, global_predictions, local_predictions, final_metrics_df, local_vi_df, global_vi_df, scaled_training_df, predictions_historic, predictions_intermediate, predictions_extreme

#%% train, split, and fit
def split_train_fit(df, response):
    """
    Function to prepare cleaned predictor dataset for model training and evaluation
    
    """
    print(f'splitting data and training model for {response}')
    
    #filter to response
    df_filtered = df.loc[df['Species']==f'{response}']
    response_column = df_filtered[f'{y_variable}']
    
    df_filtered.drop('Species', axis=1, inplace = True)
    df_filtered2 = df_filtered.drop(f'{y_variable}', axis=1, inplace = False)

    #split
    X_train, X_test, y_train, y_test = train_test_split(df_filtered2, response_column, test_size=0.4, random_state=42)
    xy_coords = X_train[['x', 'y']]
    coords_test = X_test[['x', 'y']]

    #find optimal bandwidth and local model weight (incremental spatial autocorrelation)
    bandwidth, local_weight, p_val = PyGRF.search_bw_lw_ISA(y_train, xy_coords)
    
    not_predictors = [y_variable, 'PSP', 'x' 'y']
    X_columns = [col for col in df_filtered.columns if col not in not_predictors]

    #predictors only
    X_train_IV = X_train[X_columns]
    X_test_IV = X_test[X_columns]
    #standardize independent vars
    training_stat = X_train_IV.describe().transpose()
    X_train_scaled = standardize(X_train_IV, training_stat)
    X_test_scaled = standardize(X_test_IV, training_stat)
    
    #use these best params over folds after finding out which tuning grid is best
    best_params = hyperparam_grid(X_columns, bandwidth, local_weight, 
                                  X_train_scaled, y_train, X_test_scaled, 
                                  y_test, xy_coords, coords_test)

    #call the folding function using df_filtered bc it has our response column
    #rmses = kfoldit(df_filtered, X_columns, bandwidth, best_params, X_train_scaled, y_train,
    #                X_test_scaled, y_test, xy_coords, coords_test, local_weight, response)
    #print(f'rmses for {response}: {rmses}')#should have narrow range
    #average_rmse = sum(rmses) / len(rmses)
    #print(f'average RMSE for {response}: {average_rmse}')

    #predictions
    summary_data  = predict_and_summarize(bandwidth, best_params, local_weight, 
                                          X_train_scaled, y_train, X_test_scaled,
                                          y_test, xy_coords, response)
    
    #plots
    
    return summary_data

    
#%% specify hyperparams to append to dfs
version = 1
#trees = 100
#max_features = 1
random_seed = 'n'


# %% store metrics
os.chdir("F:/remoteSensing/data/output/csv/")

model_results = [split_train_fit(full_df, response) for response in response_variables]
list_of_models = [result[0] for result in model_results]
list_of_combined_predictions = [result[1] for result in model_results]
list_of_global_predictions = [result[2] for result in model_results]
list_of_local_predictions = [result[3] for result in model_results]
list_of_final_metrics = [result[4] for result in model_results]
list_of_local_varImp = [result[5] for result in model_results]
list_of_global_varImp = [result[6] for result in model_results]
list_of_scaled_training_dfs = [result[7] for result in model_results]

#with open(f'RFmodels_v{version}.pkl', 'wb') as model_file:
#    pickle.dump(list_of_models, model_file)

all_metrics_df = pd.concat(list_of_final_metrics, ignore_index=True)
local_vi_df = pd.concat(list_of_local_varImp, ignore_index=True)
global_vi_df = pd.concat(list_of_global_varImp, ignore_index=True)
scaled_train_df = pd.concat(list_of_scaled_training_dfs, ignore_index=True)

#all_metrics_df['trees'] = trees
#all_metrics_df['mtry'] = max_features
all_metrics_df['setSeed'] = random_seed

#save
all_metrics_df.to_csv(f'model_metrics_v{version}.csv', index=False)
#local_vi_df.to_csv(f'local_feature_importances_v{version}.csv', index=False)
#global_vi_df.to_csv(f'global_feature_importances_v{version}.csv', index=False)
