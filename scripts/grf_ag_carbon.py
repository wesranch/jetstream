#%% libraries
import os
import glob
import pandas as pd
import numpy as np
import rasterio as rast
import seaborn as sns
from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import root_mean_squared_error, r2_score
from PyGRF import PyGRF #geographic rf
import pickle

"""

Author: @wesranch
Date: 2025-03-03

This script combines processed Landsat data (2000-2024), 
historic climate data from PRISM (normals 1971-2000), 
and topographic and permafrost data to train and fit random forest models.

"""

#files (raster values at plot locations)
path = "/Users/wancher/Documents/thesis/data/output/pixel-vals-3seas*.csv"
files = glob.glob(path)
climate_df = pd.read_csv("/Users/wancher/Documents/thesis/data/output/pixel-vals-hist-clim.csv")
permafrost_df = pd.read_csv("/Users/wancher/Documents/thesis/data/output/pixel-vals-perm.csv")
y_variable = 'Carbon_gm2'
response_variables = ['black spruce', 'white spruce', 'quaking aspen', 'resin birch']

#%% process, clean, join
def clean_dataframe(file):
    df = pd.read_csv(file)
    df[f'{y_variable}'] = df['Biogm2']*0.47
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
full_df = full_df.loc[:, ~full_df.columns.duplicated()].copy().drop(['year','Year'], axis=1)
print(full_df.columns)

#%% #standardize and remove correlated vars
def standardize(data, stats):
    return (data - stats['mean']) / stats['std']

# correlation
def remove_correlated(df, stats):
    corr_matrix = df.corr().abs()
    upperT = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
    drop = [col for col in upperT.columns if any(upperT[col]>0.80)]
    df_filtered = df.drop(columns=drop)
    df_standardized = standardize(df_filtered, stats)
    #plot
    #sns.heatmap(df_filtered, annot=False, cmap='coolwarm', fmt=".2f")
    return df_standardized

#%% split, train, and fit
def train_and_fit_model(df, response):
    #filter to resonse
    df_filtered = df.loc[df['Species']==f'{response}']
    response_column = df_filtered[f'{y_variable}']

    df_filtered.drop('Species', axis=1, inplace = True)
    df_filtered2 = df_filtered.drop(f'{y_variable}', axis=1, inplace = False)

    #split
    X_train, X_test, y_train, y_test = train_test_split(df_filtered2, response_column, test_size=0.3, random_state=1999)
    xy_coords = X_train[['x', 'y']]
    #coords_test = X_test[['x', 'y']]

    #find optimal bandwidth and local model weight (incremental spatial autocorrelation)
    bandwidth, local_weight, p_val = PyGRF.search_bw_lw_ISA(y_train, xy_coords)

    #cross validation
    y_predict = []
    y_true = []
    local_vi_df = pd.DataFrame()
    global_vi_df = pd.DataFrame()
    metrics_dfs = []
    models = []
    combined_predictions = []
    global_predictions = []
    local_predictions = []

    kf = KFold(n_splits=5, shuffle=True, random_state=1999)
    not_predictors = [y_variable, 'PSP']
    X_columns = [col for col in df_filtered.columns if col not in not_predictors]

    for i, (train_index, test_index) in enumerate(kf.split(df_filtered)):
        print(response, ":")
        print(f'fold {i+1}')

        X_train_fold, X_test_fold = df_filtered.iloc[train_index], df_filtered.iloc[test_index]
        y_train_fold, y_test_fold = X_train_fold[y_variable], X_test_fold[y_variable]
        X_train_IV = X_train_fold[X_columns]
        X_test_IV = X_test_fold[X_columns]
        xy_coords_train = X_train_IV[['x', 'y']]
        xy_coords_test = X_test_IV[['x', 'y']]

        #standardize independent vars
        training_stat = X_train_IV.describe().transpose()
        X_train_scaled = remove_correlated(X_train_IV, training_stat)
        X_test_scaled = remove_correlated(X_test_IV, training_stat)
        print(f'non-correlated predictor vars for {response}: ')
        print(X_train_scaled.columns)
        
        pygrf = PyGRF.PyGRFBuilder(n_estimators=48, max_features=2, 
                                   band_width=bandwidth, train_weighted=True, 
                                   predict_weighted=True, bootstrap=False,
                                   resampled=True, random_state=1999)

        #Fit the created PyGRF model based on training data and their spatial coordinates
        #X_train_corrRM = remove_correlated(X_train) #needs checked			  
        pygrf.fit(X_train_scaled, y_train_fold, xy_coords_train)
        predict_combined, predict_global, predict_local = pygrf.predict(X_test_scaled, xy_coords_test, local_weight=local_weight)

        #local variable importance
        local_model_vip = pd.DataFrame(pygrf.get_local_feature_importance())
        local_model_vip['Fold'] = i + 1
        local_model_vip['Species'] = response

        #global rf model var importance
        global_model_features = pygrf.global_model.feature_names_in_
        global_model_importance_scores = pygrf.global_model.feature_importances_
        global_vip_df = pd.DataFrame({
            'Feature': global_model_features,
            'Importance': global_model_importance_scores,
            'Fold': i+1,
            'Species': response
        })  
        
        global_vip_df = global_vip_df.sort_values(by='Importance',ascending=False).reset_index(drop=True)
        local_vi_df = pd.concat([local_vi_df, local_model_vip], ignore_index=True)
        global_vi_df = pd.concat([global_vi_df, global_vip_df], ignore_index=True)

        #store predictions
        y_predict.extend(predict_combined)
        y_true.extend(y_test_fold.tolist())

        rmse = root_mean_squared_error(y_test_fold, predict_combined)
        rSquared = r2_score(y_test_fold, predict_combined)
        
        #store metrics
        metrics_dict = {'RMSE': float(rmse), 'R2': float(rSquared), 'Fold': int(i+1), 'Species': response}
        metrics_dfs.append(metrics_dict)
      
        #print("rmse: " + str(round(rmse, 3)), "r2: " + str(round(rSquared, 3)))

        #store models and predictions
        models.append(pygrf)
        combined_predictions.extend(predict_combined)
        global_predictions.extend(predict_global)
        local_predictions.extend(predict_local)

        pygrf = None
        predict_combined = None
        predict_global = None
        predict_local = None


    final_metrics_df = pd.DataFrame(metrics_dfs)

    return models, combined_predictions, global_predictions, local_predictions, final_metrics_df, local_vi_df, global_vi_df

#%% specify hyperparams to append to dfs
version = 2
trees = 48
max_features = 2
random_seed = 'y'

# %% store metrics
os.chdir("/Users/wancher/Documents/thesis/data/output/")
model_results = [train_and_fit_model(full_df, response) for response in response_variables]
list_of_models = [result[0] for result in model_results]
list_of_combined_predictions = [result[1] for result in model_results]
list_of_global_predictions = [result[2] for result in model_results]
list_of_local_predictions = [result[3] for result in model_results]
list_of_final_metrics = [result[4] for result in model_results]
list_of_local_varImp = [result[5] for result in model_results]
list_of_global_varImp = [result[6] for result in model_results]

#with open(f'RFmodels_v{version}.pkl', 'wb') as model_file:
#    pickle.dump(list_of_models, model_file)

all_metrics_df = pd.concat(list_of_final_metrics, ignore_index=True)
local_vi_df = pd.concat(list_of_local_varImp, ignore_index=True)
global_vi_df = pd.concat(list_of_global_varImp, ignore_index=True)
all_metrics_df['trees'] = trees
all_metrics_df['mtry'] = max_features
all_metrics_df['setSeed'] = random_seed

#save
#all_metrics_df.to_csv(f'model_metrics_v{version}.csv', index=False)
#local_vi_df.to_csv('local_feature_importances.csv', index=False)
#global_vi_df.to_csv('global_feature_importances.csv', index=False)
# %%
