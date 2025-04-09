#%% libraries
import os
import glob
import pandas as pd
import numpy as np
import dask.array as da
import rioxarray
import rasterio as rast
import xarray as xr
import geopandas as gpd
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
from pyproj import Proj, transform #ensure training data xy are in same units as raster xy
import gc

import pickle
from sklearn.inspection import PartialDependenceDisplay
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings("ignore", message="X does not have valid feature names, but RandomForestRegressor was fitted with feature names")
warnings.filterwarnings("ignore", message="The weights matrix is not fully connected:")

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
#the sampled values are straight from GEE reduce region (mean reducer) at plot locations
#export type was as CSV to google drive 
path = "F:/remoteSensing/data/output/csv/Sampled_Alaska/pixel-vals-fire*.csv"
files = glob.glob(path)
climate_df = pd.read_csv("F:/remoteSensing/data/output/csv/pixel-vals-hist-clim.csv")
climate_df['pr'] = climate_df['hist_pr_1971_2000']
climate_df['tas'] = climate_df['hist_tas_1971_2000']
climate_df.drop(['hist_tas_1971_2000', 'hist_pr_1971_2000'], axis=1, inplace=True)
permafrost_df = pd.read_csv("F:/remoteSensing/data/output/csv/pixel-vals-perm.csv")
y_variable = 'Carbon_gm2'
response_variables = ['black spruce', 'white spruce', 'quaking aspen', 'resin birch']

#patterns
mosaics_historic = "F:/remoteSensing/data/output/gee/mosaics/landisReg-FULLmosaic-hist-*.tif"
mosaics_intermediate = "F:/remoteSensing/data/output/gee/mosaics/landisReg-FULLmosaic-ncar-*.tif"
mosaics_extreme = "F:/remoteSensing/data/output/gee/mosaics/landisReg-FULLmosaic-gfdl-*.tif"
output_dir = "F:/remoteSensing/data/output/raster/April_GRF_Predictions/"

#reference raster
raster_template = rast.open("F:/remoteSensing/data/output/gee/mosaics/landisReg-FULLmosaic-hist-2015.tif")
bandnames = list(raster_template.descriptions)
bandnames.append('x')
bandnames.append('y')

#shpfile
bounds = gpd.read_file('F:/randomForest_Landis/input/shp/FullLandscapeV3_082722.shp')

#list files
raster_files_historic = glob.glob(mosaics_historic)
raster_files_intermediate = glob.glob(mosaics_intermediate)
raster_files_extreme = glob.glob(mosaics_extreme)

#%%  #raster manipulation functions
def save_to_tiff(data, output_path, raster_template):
    data = data.compute()
    transform = raster_template.transform
    crs = raster_template.crs
    height, width = data.shape
    data.rio.to_raster(output_path, transform=transform, crs=crs, dtype='float32', 
                       tiled=True, windowed=True, creation_options={'BIGTIFF': 'YES'})
   
def read_raster(file):
    """
    Reads in each band of landsat composites + climate rasters and adds xy as bands
    These are needed as x and y are used as predictors in the models
    """
    rds = rioxarray.open_rasterio(file, lock=False, cache=False, chunks={'x': 100, 'y': 100})
    #rds = rds.rio.clip(bounds.geometry)
    #rds = rioxarray.open_rasterio(file, lock=False, cache=False)
    
    #metadata
    transform = rds.rio.transform()
    height, width = rds.shape[1], rds.shape[2]
    rds.name = "bands"

    
    #generate coords
    x_coords = xr.DataArray(rds.x, name="x")
    y_coords = xr.DataArray(rds.y, name="y")
  
    #convert to 2d array
    x_coords_2d = x_coords.expand_dims(dim={'y': height})
    y_coords_2d = y_coords.expand_dims(dim={'x': width})
    #x_coords_2d = x_coords_2d.transpose('y', 'x')
    #y_coords_2d = y_coords_2d.transpose('y', 'x')
    
    #3d array for stacking
    x_coords_3d = xr.DataArray(x_coords_2d.expand_dims(dim={'band': [0]}), name='x')
    y_coords_3d = xr.DataArray(y_coords_2d.expand_dims(dim={'band': [0]}), name='y')
    x_coords_3d = x_coords_3d.transpose('band', 'y', 'x')
    y_coords_3d = y_coords_3d.transpose('band', 'y', 'x')


    
    return rds, x_coords_3d, y_coords_3d, transform, width, height

def plot_single_band(rds, band_index):
    """
    Plots the specified slope band from the raster data.
    
    Parameters:
    - rds: The raster data (xarray.Dataset).
    - band_index: The index of the band to plot (default is 0).
    """
    # Extract the band
    band = rds.isel(band=band_index)
    data = band.values
    # Plot the band
    plt.figure(figsize=(10, 10))
    plt.imshow(data, cmap='viridis')  
    plt.title("slope")
    plt.show()

#%% process, clean, join predictor data for the model
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

#reorder columns to same as raster and convert xy from DD to UTM
proj_in = Proj(init='epsg:4326')
proj_out = Proj(init='epsg:32612')#utm zone 12n
transformer = Transformer.from_proj(proj_in, proj_out)
full_df['easting'], full_df['northing'] = transformer.transform(full_df['x'], full_df['y'])
full_df['x'] = full_df['easting']
full_df['y'] = full_df['northing']
full_df.drop(['easting', 'northing'], axis=1, inplace = True)

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
    Used in tuning the models.
    """
    param_grid = {
        'n_estimators': np.arange(40, 100, 5),
        'max_features': np.arange(1, min(6, len(predictors)) + 1, 1)
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

#predict on raster
def predict_and_save(file, response, build, xtrain, ytrain, lw):
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
    
    #apply read raster function
    print(f"reading raster data for year {yr} for {scenario} climate")
    rds, x_coords_3d, y_coords_3d, transform, width, height = read_raster(file)

    #splitting vertically
    size = 2000 
    #chunk_x = 0
    #chunk_y = 0
    for chunk_x in range(0, width, size):  
    #    for chunk_y in range(0, height, size):
        #omit for chunking and use the full size (too long of a runtime)
        chunk_width = size
        chunk_height = height
        if chunk_x + chunk_width > width:
            chunk_width = width - chunk_x  # Adjust the last chunk width
        window = (slice(0, chunk_height), slice(chunk_x, chunk_x + chunk_width))
        
        #chunk_width = min(size, width - chunk_x)
        #chunk_height = min(size, height - chunk_y) 
        #assert chunk_width == size or chunk_x + chunk_width == width
        #assert chunk_height == size or chunk_y + chunk_height == height
        
        #window = (slice(chunk_y, chunk_y + chunk_height), slice(chunk_x, chunk_x + chunk_width))
        chunk = rds.isel(x=window[1], y=window[0])
        #plot_single_band(chunk, band_index=33)
        x_chunk = x_coords_3d.isel(x=window[1], y=window[0])
        y_chunk = y_coords_3d.isel(x=window[1], y=window[0])
        print(f"chunk shape: {chunk.shape}")
        print(f"x_chunk shape: {x_chunk.shape}")
        print(f"y_chunk shape: {y_chunk.shape}")
        
        #convert to 2d array
        x_chunk_2d = x_chunk.squeeze(dim='band', drop=True)
        y_chunk_2d = y_chunk.squeeze(dim='band', drop=True)
        #x_broadcast, _ = xr.broadcast(chunk, x_chunk_2d)
        #y_broadcast, _ = xr.broadcast(chunk, y_chunk_2d)
        
        #stack bands with x and y bands we created
        band_indices = chunk.band.values
        x_broadcast = x_chunk.assign_coords(band=("band", [len(band_indices)+1]))
        y_broadcast = y_chunk.assign_coords(band=("band", [len(band_indices)+2]))
        stacked_chunk = xr.concat([chunk, x_broadcast, y_broadcast], dim="band")
        stacked_chunk.attrs["long_name"] = bandnames
        #x and y need to be clipped
        #clipped_stack = stacked_chunk.rio.clip(bounds.geometry)
        
        #convert stack to 2d array (n pixels, n features) for randomForest predictions
        transposed = stacked_chunk.transpose("y", "x", "band")
        nrows, ncols, nbands = transposed.shape
        stacked_2d = transposed.values.reshape(nrows*ncols, nbands)
        
        #xy coords must match the raster so these are different than those used for fitting
        xy_coords = np.column_stack([x_chunk_2d.values.ravel(),y_chunk_2d.values.ravel()])
        assert stacked_2d.shape[0] == xy_coords.shape[0], f"Shape mismatch in chunk {chunk_x}, {chunk_y}"
        
        #predictions (not sure what it will do with the NAs in raster data)
        output_path_combined = f"{output_dir}COMBINED/Carbon_{spp}_{scenario}_{yr}.tif"
        output_path_global = f"{output_dir}GLOBAL/Carbon_{spp}_{scenario}_{yr}.tif"
        output_path_local = f"{output_dir}LOCAL/Carbon_{spp}_{scenario}_{yr}.tif"
        
        #skip tiles that already exist
        # if (
        #     os.path.exists(output_path_combined)
        #     and os.path.exists(output_path_global)
        #     and os.path.exists(output_path_local)
        # ):
        #     print(f"Chunk {chunk_x}, {chunk_y} already processed. Skipping...")
        #     continue
        
        predict_combined, _, _ = build.predict(stacked_2d, xy_coords, local_weight=lw)
        #assuming we want combined (reshape to 2d then convert to xarray)
        predict_combined_array = np.array(predict_combined)
        #predict_local_array = np.array(predict_local)
        
        prediction_combined_2d = predict_combined_array.reshape((nrows, ncols))
        #prediction_local_2d = predict_local_array.reshape((nrows,ncols))
        #prediction_global_2d = predict_global.reshape((nrows, ncols))
        
        #xr data arrays
        prediction_combined_xr_array = xr.DataArray(
            prediction_combined_2d,
            dims=("y", "x"),
            coords={"x": stacked_chunk.coords["x"], "y": stacked_chunk.coords["y"]},
            attrs={"long_name": "Predicted_C_gm2"}
        )
        # just doing combined since comparing to landis so commenting out other returns
        # prediction_local_xr_array = xr.DataArray(
        #     prediction_local_2d,
        #     dims=("y", "x"),
        #     coords={"x": stacked_chunk.coords["x"], "y": stacked_chunk.coords["y"]},
        #     attrs={"long_name": "Predicted_C_gm2"}
        # )
        # prediction_global_xr_array = xr.DataArray(
        #     prediction_global_2d,
        #     dims=("y", "x"),
        #     coords={"x": stacked_chunk.coords["x"], "y": stacked_chunk.coords["y"]},
        #     attrs={"long_name": "Predicted_C_gm2"}
        # )
        
        #export
        save_to_tiff(prediction_combined_xr_array, output_path_combined, raster_template)
        #save_to_tiff(prediction_global_xr_array, output_path_global, raster_template)
        #save_to_tiff(prediction_local_xr_array, output_path_local, raster_template)
        
        #print(f"Processed and saved chunk: {chunk_x}, {chunk_y}")
        
        #clean things up
        gc.collect()
        predict_combined=None
        #predict_global=None
        #predict_local=None
        predict_combined_array=None
        #predict_local_array=None
        prediction_combined_2d=None
        #prediction_global_2d=None
        #prediction_local_2d=None
        prediction_combined_xr_array=None
        #prediction_global_xr_array=None
        #prediction_local_xr_array=None
        
        chunk=None
        x_chunk=None
        y_chunk=None
        x_chunk_2d=None
        y_chunk_2d=None
        x_broadcast=None
        y_broadcast=None
        stacked_chunk=None
        #clipped_stack=None
        stacked_2d=None
        xy_coords=None
        transposed=None

            
    print(f"processed {file} for {spp}")
    rds=None
    x_coords_3d=None
    y_coords_3d=None
    transform=None

#%% running in parallel when predicting
from concurrent.futures import ProcessPoolExecutor

def run_parallel_predictions(raster_files, response, model, xtrain, ytrain, lw):
    with ProcessPoolExecutor(max_workers=2) as executor:
        future_to_file = {
            executor.submit(predict_and_save, file, response, model, xtrain, ytrain, lw): file
            for file in raster_files
        }
        results = []
        for future, file in future_to_file.items():
            try:
                
                results.append(future.result())
            except Exception as e:
                print(f"Error processing {file}: {e}")
        return results
        
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
    pygrf.fit(xtrain, ytrain, xy_coords)
    
    #predict and save rasters
    #non parallel method
    predictions_historic = list(map(lambda file: predict_and_save(file, response, pygrf, xtrain, ytrain, lw), raster_files_historic))
    predictions_intermediate = list(map(lambda file: predict_and_save(file, response, pygrf, xtrain, ytrain, lw), raster_files_intermediate))
    predictions_extreme = list(map(lambda file: predict_and_save(file, response, pygrf, xtrain, ytrain, lw), raster_files_extreme))
    
    #predictions_historic = run_parallel_predictions(raster_files_historic, response, pygrf, xtrain, ytrain, lw)
    #predictions_intermediate = run_parallel_predictions(raster_files_intermediate, response, pygrf, xtrain, ytrain, lw)
    #predictions_extreme = run_parallel_predictions(raster_files_extreme, response, pygrf, xtrain, ytrain, lw)
    
    #predictions on test data for r2 graphs 
    predict_combined, predict_global, predict_local = pygrf.predict(test_data, xy_coords, local_weight=lw)
    
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

    return models, combined_predictions, global_predictions, local_predictions, final_metrics_df, local_vi_df, global_vi_df, scaled_training_df

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
    df_filtered2 = df_filtered2[bandnames]
    
    #split
    X_train, X_test, y_train, y_test = train_test_split(df_filtered2, response_column, test_size=0.2, random_state=42)
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
