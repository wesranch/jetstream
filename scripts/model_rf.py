#%% libraries
import os
import glob
import pandas as pd
import numpy as np
import rasterio as rast
import seaborn as sns
from sklearn.model_selection import train_test_split
from PyGRF import PyGRF #geographic rf

"""

Author: @wesranch
Date: 2025-03-03

This script combines processed landsat data (2000-2024), 
historic climate data from PRISM (normals 1971-2000), 
and topographic and permafrost data to train and fit random forest models.

"""

#files (raster values at plot locations)
path = "/Users/wancher/Documents/thesis/data/output/pixel-vals-3seas*.csv"
files = glob.glob(path)
climate_df = pd.read_csv("/Users/wancher/Documents/thesis/data/output/pixel-vals-hist-clim.csv")
permafrost_df = pd.read_csv("/Users/wancher/Documents/thesis/data/output/pixel-vals-perm.csv")

#%% process
def clean_dataframe(file):
    df = pd.read_csv(file)
    df['Carbon_gm2'] = df['Biogm2']*0.47
    df[['x','y']] = df['.geo'].str.extract(r'\[([-0-9.]+),([-0-9.]+)\]').astype(float)
    df.drop(['system:index','Biogm2','SPP_count','.geo'], axis=1, inplace = True)
    df.dropna(inplace=True)
    print(df.shape[1])
    return df

dfs = list(map(clean_dataframe,files))#lapply equivalent

#remove dfs that have less columns than they should have
dfs = [df for df in dfs if df.shape[1] == 36]
landsat_df = pd.concat(dfs,axis=0)#rbind
landsat_df = landsat_df.reset_index(drop=True)

#climate and perm
clim_perm_df = pd.concat([climate_df, permafrost_df],axis=1)
clim_perm_df = clim_perm_df.loc[:, ~clim_perm_df.columns.duplicated()].copy()
clim_perm_df = clim_perm_df.reset_index(drop=True)

#join
full_df = pd.concat([landsat_df, clim_perm_df], axis=1).dropna(inplace=False)
full_df = full_df.loc[:, ~full_df.columns.duplicated()].copy().drop('year', axis=1)
print(full_df.columns)

#%% correlation
def remove_correlated(df):
    corr_matrix = df.corr()

    #plot
    #sns.heatmap(df_corrRM, annot=False, cmap='coolwarm', fmt=".2f")

    #nothing removed while testing with this
    df_corrRM = corr_matrix[(corr_matrix < 0.70) & (corr_matrix > -0.70)]
    return df_corrRM

#%% split, train, and fit
def train_model(df, response):
    #filter to resonse
    df_filtered = df.loc[df['Species']==f'{response}']
    response_column = df_filtered['Carbon_gm2']

    df_filtered.drop('Species', axis=1, inplace = True)
    df_filtered.drop('Carbon_gm2', axis=1, inplace = True)

    #split
    X_train, X_test, y_train, y_test = train_test_split(df_filtered, response_column, test_size=0.2, random_state=42)
    xy_coords = X_train[['x', 'y']]
    coords_test = X_test[['x', 'y']]

    #Create a PyGRF model by specifying hyperparameters
    pygrf = PyGRF.PyGRFBuilder(n_estimators=500, max_features=1, band_width=39, train_weighted=True, predict_weighted=True, bootstrap=False,
                            resampled=True, random_state=42)

    #Fit the created PyGRF model based on training data and their spatial coordinates
    #X_train_corrRM = remove_correlated(X_train) #needs checked			  
    pygrf.fit(X_train, y_train, xy_coords)

    #Make predictions for testing data using the fitted PyGRF model
    predict_combined, predict_global, predict_local = pygrf.predict(X_test, coords_test, local_weight=0.46)
    
    return pygrf, predict_combined, predict_global, predict_local

# %% store the models
# plan to update to store plots and model metrics
response_variables = ['black spruce', 'white spruce', 'quaking aspen', 'resin birch']
models = [train_model(full_df, response) for response in response_variables]
global_fi = models[1][0].global_model.feature_importances_#add in the function
# %% predicting and storing
# hyperparameterize and do cross folds
