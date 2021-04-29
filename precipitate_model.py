from sklearn.ensemble import RandomForestClassifier as RF
from sklearn.model_selection import GridSearchCV, cross_val_score, KFold, train_test_split
from sklearn.preprocessing import StandardScaler

import pandas as pd
import numpy as np
from skopt import BayesSearchCV
from skopt.space import Real, Categorical, Integer

## ------ Get transforms for data ---------------------------------------
def transform_training(df, transform = None):
#
# standardize
  cols = df.columns
  if transform == None:
    transform = StandardScaler()
    transform.fit(df)
  df = pd.DataFrame(transform.transform(df), columns = cols)
#
  return(df, transform)

## -------- Bayesian Optimization ---------------------------------
def CV_opt(params, mod, X, y):
  ''' optimize fit '''
# make RF model and format y
  y = np.ravel(y)
#
  bayes_search = BayesSearchCV(mod, params, n_iter=32,
                               scoring="f1", n_jobs=10, cv=5)
  bayes_search.fit(X, y)
  return(bayes_search)

## --------- Model Evaluation -----------------------------------------  
from sklearn import metrics

def eval(model, test_X, test_y, train_X, train_y):
#
  y_pred = model.predict(test_X)
  accuracy = metrics.accuracy_score(test_y, y_pred)
  f1 = metrics.f1_score(test_y, y_pred)
  print("VALIDATION: accuracy: ", accuracy, " F1: ", f1)
  print("VALIDATION: recall: ", metrics.recall_score(test_y, y_pred),
        " precision: ", metrics.precision_score(test_y, y_pred))
#
  y_pred = model.predict(train_X)
  accuracy = metrics.accuracy_score(train_y, y_pred)
  f1 = metrics.f1_score(train_y, y_pred)
  print("TRAINING: accuracy: ", accuracy, " F1: ", f1)
  

## -------- load training data ------------------------
labels = "PAWIIR_precipitated_cells.csv"
Y_labels = pd.read_csv(labels)
Y = Y_labels["precip"].astype("int")
features = "2021-04-20_PAWIIR_features_10060x10.txt"
X = pd.read_csv(features, sep = "\t", header = 0).iloc[:,1:11]

#CART_labels = "CART4_Precipitated_cells.csv"
#Y_CART_labels = pd.read_csv(CART_labels)
#CART_Y = Y_CART_labels["precip"].astype("int")
#
#CART_features = "2021-04-13_CART4_features_4807x10.txt"
#CART_X = pd.read_csv(CART_features, sep = "\t", header = 0).iloc[:, 1:11]
#
#X = pd.concat([X, CART_X], axis = 0, ignore_index = True)
#Y = pd.concat([Y, CART_Y], axis = 0, ignore_index = True)

X_train, X_test, y_train, y_test = train_test_split(X,Y, test_size=0.25, random_state=19143)

## transfrom train set, then use same scale for test set
## doesn't need to be transformed for RF, but it just makes it easier for also doing SVM etc.
X_train, x_scale = transform_training(X_train)
X_test = pd.DataFrame(x_scale.transform(X_test), columns = X_train.columns)


## ------------- RANDOM FOREST ---------------------------
PARAMS = {"bootstrap": Categorical([True, False]), 
        "max_depth": Integer(3, 20), 
        "max_features": Categorical(['auto', 'sqrt','log2']), 
        "min_samples_leaf": Integer(2, 10),
        "min_samples_split": Integer(2, 10),
        "n_estimators": Integer(100, 500) }

## run Bayesian CV optimizer
untrained_RF = RF()
RF_CV = CV_opt(PARAMS, untrained_RF, X_train, y_train)
RF_mod = RF_CV.best_estimator_

eval(RF_mod, X_test, y_test, X_train, y_train)

pd.DataFrame(X_test.columns, RF_mod.feature_importances_)

## ----------- SVM ------------------------------------
from sklearn.svm import SVC

PARAMS = {'C': (1e-3, 1e+3, 'log-uniform'),
        'gamma': (1e-3, 1e+1, 'log-uniform'),
        'degree': (1, 3),  
        'kernel': ['linear', 'poly', 'rbf']} 

## run Bayesian CV optimizer
untrained_SVM = SVC()
SVM_CV = CV_opt(PARAMS, untrained_SVM, X_train, y_train)
SVM_mod = SVM_CV.best_estimator_

eval(SVM_mod, X_test, y_test, X_train, y_train)

## ---------- Adaboost -----------------------------
from sklearn.ensemble import AdaBoostClassifier as ABC

param_grid = { "learning_rate": (1e-3, 1e0, 'log-uniform'),
              "n_estimators": Integer(5, 150) }

base_SVC = SVC(kernel = "linear", random_state = 19002)

## run Bayesian CV optimizer
untrained_ABC = ABC(base_estimator = base_SVC, random_state = 19143, algorithm = "SAMME")
ABC_CV = CV_opt(param_grid, untrained_ABC, X_train, y_train)
ABC_mod = ABC_CV.best_estimator_

eval(ABC_mod, X_test, y_test, X_train, y_train)


## ---------- Adaboost W/ decision stump -----------------------------
from sklearn.ensemble import AdaBoostClassifier as ABC
from sklearn.tree import DecisionTreeClassifier as DTC

param_grid = {"base_estimator__criterion" : ["gini", "entropy"],
              "base_estimator__splitter" :   ["best", "random"],
              "learning_rate": (1e-4, 1e0, 'log-uniform'),
              "n_estimators": Integer(5, 500)
             }

base_DTC = DTC(random_state = 19002, max_depth = 1)

## run Bayesian CV optimizer
untrained_ABC_stump = ABC(base_estimator = base_DTC, random_state = 19143)
ABC_stump_CV = CV_opt(param_grid, untrained_ABC_stump, X_train, y_train)
ABC_stump_mod = ABC_stump_CV.best_estimator_

eval(ABC_stump_mod, X_test, y_test, X_train, y_train)


## ============= ALL PREDICTIONS ========================
predictions = pd.DataFrame(np.asarray( (RF_mod.predict(X_test), SVM_mod.predict(X_test),
                           ABC_mod.predict(X_test), y_test ) ).T )
predictions.columns = ["RF", "SVM", "Adaboost", "true"]
predictions["cell_name"] = Y_labels.iloc[y_test.index]["cells"].tolist()

predictions[predictions["true"] == 1]

predictions.to_csv("each_model_predictions.csv")

## ----------------------------------------------
#labels = "CART4_Precipitated_cells.csv"
labels2 = "PAUMXB_precipitated_cells.csv"
Y2_labels = pd.read_csv(labels2)
Y2 = Y2_labels["precip"].astype("int")
#features = "2021-04-13_CART4_features_4807x10.txt"
features2 = "2021-04-20_PAUMXB_features_21982x10.txt"
X2 = pd.read_csv(features2, sep = "\t", header = 0).iloc[:,1:11]

## transfrom train set, then use same scale for test set
## doesn't need to be transformed for RF, but it just makes it easier for also doing SVM etc.
X2, x_scale = transform_training(X2)

y2_pred = RF_mod.predict(X2)
accuracy = metrics.accuracy_score(Y2, y2_pred)
f1 = metrics.f1_score(Y2, y2_pred)
print("TEST: accuracy: ", accuracy, " F1: ", f1)
print("TEST: recall: ", metrics.recall_score(Y2, y2_pred),
        " precision: ", metrics.precision_score(Y2, y2_pred))

y2_pred = ABC_mod.predict(X2)
accuracy = metrics.accuracy_score(Y2, y2_pred)
f1 = metrics.f1_score(Y2, y2_pred)
print("TEST: accuracy: ", accuracy, " F1: ", f1)
print("TEST: recall: ", metrics.recall_score(Y2, y2_pred),
        " precision: ", metrics.precision_score(Y2, y2_pred))

y2_pred = ABC_stump_mod.predict(X2)
accuracy = metrics.accuracy_score(Y2, y2_pred)
f1 = metrics.f1_score(Y2, y2_pred)
print("TEST: accuracy: ", accuracy, " F1: ", f1)
print("TEST: recall: ", metrics.recall_score(Y2, y2_pred),
        " precision: ", metrics.precision_score(Y2, y2_pred))



