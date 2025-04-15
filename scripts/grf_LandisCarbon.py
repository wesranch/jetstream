#%% libraries
import os
import glob
import pandas as pd
import numpy as np
import rasterio as rast
from itertools import product
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
Date: 2025-04-01

This script uses simulated output from LANDIS-II
to train and fit geographically weighted random forest models.

Not currently set to remove correlated vars!

"""

#files (raster values at plot locations)
path = "F:/randomForest_Landis/output/main_predictor_data.csv"
main_df = pd.read_csv(path)

y_variable = 'ag_carbon'
response_variables = ['BlackSpruce', 'WhiteSpruce', 'QuakingAspen', 'PaperBirch']

#%% process, clean, join
def clean_dataframe(df):
    #should be some duplicates (because of grouping by climate region)
    # but not as many that came out of R
    df_unique = df.drop_duplicates(subset=['spp', 'PSP', 'ClimateMap'], keep='first')
    #could add additional climate (just historic here)
    df_unique.drop(['PSP', 'scenario'], axis=1, inplace = True)
    df_unique.dropna(inplace=True)
    print(df_unique.shape[1])
    
    return df_unique

clean_df = clean_dataframe(main_df)

#reset index column
landis_df = clean_df.reset_index(drop=True)
print(landis_df.columns)


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
    param_grid = {
        'n_estimators': np.arange(40, 100, 10),
        'max_features': np.arange(1, min(5, len(predictors)) + 1, 2)
    }
    param_combinations = list(product(param_grid['n_estimators'], param_grid['max_features']))
    best_rmse = float('inf')
    best_params = None

    for n_estimators, max_features in param_combinations:
        pygrf = PyGRF.PyGRFBuilder(band_width=bw, train_weighted=True, predict_weighted=True, 
                                   bootstrap=False, resampled=True, random_state=42,
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
            x_test_df, y_test_df, coords, coords_test, lw):
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
        pygrf.fit(x_train_df, y_train_df, coords)
        predict_combined, _, _ = pygrf.predict(x_test_df, coords_test, local_weight=lw)

        rmse = root_mean_squared_error(y_test_df, predict_combined)  # RMSE
        rmses.append(rmse)#take an average of this
        
    return rmses

#%% train, split, and fit
def train_and_split(df, response):
    
    print(f'splitting data and training model for {response}')
    
    #filter to resonse
    df_filtered = df.loc[df['spp']==f'{response}']
    response_column = df_filtered[f'{y_variable}']
    
    df_filtered.drop('spp', axis=1, inplace = True)
    df_filtered2 = df_filtered.drop(f'{y_variable}', axis=1, inplace = False)

    #split
    X_train, X_test, y_train, y_test = train_test_split(df_filtered2, response_column, test_size=0.2, random_state=42)
    xy_coords = X_train[['x', 'y']]
    coords_test = X_test[['x', 'y']]

    #find optimal bandwidth and local model weight (incremental spatial autocorrelation)
    bandwidth, local_weight, p_val = PyGRF.search_bw_lw_ISA(y_train, xy_coords)
    
    not_predictors = [y_variable, 'x' 'y']
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

    #df_filtered b/c it has the ag_carbon column
    rmses = kfoldit(df_filtered, X_columns, bandwidth, best_params, X_train_scaled, y_train,
                    X_test_scaled, y_test, xy_coords, coords_test, local_weight)
    average_rmse = sum(rmses) / len(rmses)
    print(f'average RMSE for {response}: {average_rmse}')
    return average_rmse, local_weight

#store models
model_results = [train_and_split(landis_df, response) for response in response_variables]
final_rmses = [result[0] for result in model_results]
local_weights = [result[1] for result in model_results]


    
#%% #function to predict and store metrics from stored models
def predict_and_summarize(folded_models, lw, test_data, xy_coords):
    
    for i in range(folded_models):
        model = folded_models(i)
        test_data = test_data(i)
        coords = xy_coords(i)
        local_weight = lw(i)
        
        predict_combined, predict_global, predict_local = model.predict(test_data, 
                                                                        coords, 
                                                                        local_weight=local_weight)

        #local variable importance
        local_model_vip = pd.DataFrame(model.get_local_feature_importance()).mean()
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

    #store predictions
    y_predict.extend(predict_combined)
    y_true.extend(y_test.tolist())

    rmse = root_mean_squared_error(y_test, predict_combined)
    rSquared = r2_score(y_test, predict_combined)
    
    #store metrics
    metrics_dict = {'RMSE': float(rmse), 'R2': float(rSquared), 'Species': response}
    metrics_dfs.append(metrics_dict)
    
    #print("rmse: " + str(round(rmse, 3)), "r2: " + str(round(rSquared, 3)))

    #store models and predictions
    models.append(pygrf)
    combined_predictions.extend(predict_combined)
    global_predictions.extend(predict_global)
    local_predictions.extend(predict_local)
    
    #clean
    pygrf = None
    predict_combined = None
    predict_global = None
    predict_local = None


    final_metrics_df = pd.DataFrame(metrics_dfs)

    return models, combined_predictions, global_predictions, local_predictions, final_metrics_df, local_vi_df, global_vi_df, scaled_training_df

#isolate each spp
predictions = []
for response in response_variables:
    kfold_models = [result[0] for result in model_results]
    local_weights = [result[1] for result in model_results]
    kfold_test_data = [result[2] for result in model_results]
    kfold_xycoords = [result[3] for result in model_results]
    
    predictions = predict_and_summarize(kfold_models, local_weights, 
                                        kfold_test_data, kfold_xycoords)
    
#%% specify hyperparams to append to dfs
version = 3
trees = 100
max_features = 1
random_seed = 'y'


# %% store metrics
os.chdir("F:/randomForest_Landis/output/")
model_results = [train_and_fit_model(landis_df, response) for response in response_variables]
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

all_metrics_df['trees'] = trees
all_metrics_df['mtry'] = max_features
all_metrics_df['setSeed'] = random_seed

#save
all_metrics_df.to_csv(f'model_metrics_v{version}.csv', index=False)
#local_vi_df.to_csv(f'local_feature_importances_v{version}.csv', index=False)
#global_vi_df.to_csv(f'global_feature_importances_v{version}.csv', index=False)
# %% variable importance
indices = []
for response in response_variables:
    X = scaled_train_df.loc[scaled_train_df['Species']==f'{response}']
    

    #global feature importance
    numeric_df = global_vi_df.loc[global_vi_df['Species']==f'{response}'].drop(columns=['Species'], errors='ignore')
    indices= numeric_df.sort_values(by='Importance',ascending=False).head(25)
    
   #local feature importance
    numeric_df_local = local_vi_df.loc[local_vi_df['Species']==f'{response}'].drop(columns=['Species', 'model_index'], errors='ignore')
    indices_local= numeric_df_local.mean().sort_values(ascending=False).head(25)
    
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
    fig.suptitle(f"LANDIS-II Feature Importance - {response}", fontsize=16, ha='center')
    plt.tight_layout()
    plt.show()

# %%
