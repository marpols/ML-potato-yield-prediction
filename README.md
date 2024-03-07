# ML
R code for running Random Forest and hybrid model with STICS for potato yield prediction

RandomForest.R
Calls the below functions to train multiple random forest models for each individual field seperated out of the dataset as the test set.
Performs Hyperparameter tuning for nodesize if not set as NULL. 

rfFunctions.R
Functions used in RandomForest.R

Model Selection Functions:

calc_metrics(actual,predicted)
calculates the MSE, RMSE, RRMSE, R2, BIAS and PBIAS and returns them as a vector.

min_mse(feature selection results)
Find best model based on MSE for add-one-in feature selection models. The model with the lowest test MSE and also the smallest difference between testing and training MSE is selected.
returns the best model.

best_mod(models, nodesize=BOOL)
Find the best model from the results of gridsearch based on RMSE. If performing gridsearch for nodesize, set nodesize=T.
The model with the lowest test RMSE and also the smallest difference between testing and training RMSE is selected.
returns the best model.

Feature Selection:
featselec(train_x, train_y, train_x, train_y, ntree value, folds=INT)
folds - number of folds for CV
Performs an initial kfold cross validation to determine an optimum number of features from all features using the default R random forest hyperparameters.
Then performs an "add-one-in" selection method in which each feature determined from the kfold CV + 2 is added in one at a time and the errors bewtween the testing and training sets are compared.
returns a list of all models created for feature selection (Each key is the number of features and contains the model metrics and list of features used).

Hyperparameter tuning (for mtry and ntree):
hype_tune(training set (x and y values together), test_x, test_y, cross validation method="string", folds=INT,reps=INT, mtry_seq=vector, ntree_seq=vector, tune_search="string")
folds - number of folds for CV
reps - number of repetitions if doing a repeated CV
mtry_seq - sequence of mtry values to try
ntree_seq - sequence of ntree valuees to try
tune_search - type of search to perform. Grid or Random search.
Runs random forests for hyperparameter tuning based on the type of search performed. Returns a list of all models created.

rfPlots.R
generate_plot(datafile,run,method,labels="string")
datafile - .csv created in RandomForest.R containing all the actual and predicted yields for each field-Year
run - name given to current run for file creation (entered in RandomForest.R)
method - ML alone or hybrid approach (entered in RandomForest.R)
labels - label the points on the plot as individual field-years ("CFIDYr") or group by year ("Year")
Gerenerates Actual Vs. Predicted metrics and yield plots for the overall predictions.






