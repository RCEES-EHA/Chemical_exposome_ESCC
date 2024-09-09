# -*- coding: utf-8 -*-
import numpy as np
from numpy.core.numeric import indices
from numpy.lib.index_tricks import s_
import pandas as pd
import sys
import os
import pickle

from sklearn.feature_selection import VarianceThreshold, SelectFromModel, chi2
from collections import Counter
from sklearn.model_selection import StratifiedKFold, GridSearchCV,cross_validate, cross_val_predict
from sklearn.metrics import classification_report
from sklearn.svm import LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process.kernels import RBF
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, ExtraTreesClassifier,GradientBoostingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import confusion_matrix, accuracy_score, make_scorer, cohen_kappa_score
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from imblearn.pipeline import Pipeline
from imblearn.over_sampling import SMOTENC, SMOTE
from pathlib import Path
#from sklearn.pipeline import Pipeline


# import warnings filter
from warnings import simplefilter
# ignore all future warnings
simplefilter(action='ignore', category= FutureWarning)
simplefilter(action='ignore', category= UserWarning)
simplefilter(action='ignore', category= DeprecationWarning)

def tn(y_true, y_pred): return confusion_matrix(y_true, y_pred)[0, 0]
def fp(y_true, y_pred): return confusion_matrix(y_true, y_pred)[0, 1]
def fn(y_true, y_pred): return confusion_matrix(y_true, y_pred)[1, 0]
def tp(y_true, y_pred): return confusion_matrix(y_true, y_pred)[1, 1]
scoring = {'tp': make_scorer(tp), 'tn': make_scorer(tn),
           'fp': make_scorer(fp), 'fn': make_scorer(fn),
           'auc': 'roc_auc',
           'acc': make_scorer(accuracy_score),
           'kappa': make_scorer(cohen_kappa_score)}

def update_progress(progress):
    barLength = 100 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), round(progress*100, 2), status)
    sys.stdout.write(text)
    sys.stdout.flush()  


def mlTest(data, target, name_dataset):
    # Nested Cross Validation:
    inner_loop_cv = 5   
    outer_loop_cv = 5
    
    # Number of random trials:
    NUM_TRIALS = 50
    
    # Grid of Parameters:
    K_grid = {"clf__n_neighbors": [3, 5, 7, 9, 11, 13, 15]}
    C_grid = {"clf__C": [0.001, 0.01, 0.1, 1, 10, 100, 1000]}
    est_grid = {"clf__n_estimators": [2, 4, 8, 16, 32, 64]}
    MLP_grid = {"clf__alpha": [0.001, 0.01, 0.1, 1, 10], "clf__learning_rate_init": [0.001, 0.01, 0.1, 1],
        "clf__hidden_layer_sizes": [10, 20]}
    SVC_grid = {"clf__gamma": [0.0001, 0.001, 0.01, 0.1], "clf__C": [0.001, 0.01, 0.1, 1, 10, 100, 1000]}
    DT_grid = {"clf__max_depth": [10, 20, 30, 50, 100]}
    XGBoost_grid = {"clf__n_estimators": [2, 4, 8, 16, 32, 64], "clf__learning_rate": [0.001, 0.01, 0.1, 1]}
        
    # Classifiers:
    names = ["KNN", "Logistic Regression", "RBF SVM", "Extra Trees", "Random Forest", "Neural Net"]

    classifiers = [
        KNeighborsClassifier(),
        LogisticRegression(),
        SVC(),
        ExtraTreesClassifier(),
        RandomForestClassifier(),
        MLPClassifier()]
    
    # names = ["KNN", "RBF SVM", "Extra Trees", "Random Forest", "AdaBoost", "XGBoost"]

    # classifiers = [
    #     KNeighborsClassifier(),
    #     SVC(),
    #     ExtraTreesClassifier(),
    #     RandomForestClassifier(),
    #     AdaBoostClassifier(),
    #     GradientBoostingClassifier()]
            
    # # Check minimum number of samples:
    # count_class = Counter(target)
    # print(count_class)
    # if count_class[0] < 12 or count_class[1] < 12:
    #     continue

    # Remove low variance:
    print("Before removing low variance:{}".format(data.shape))
    selector = VarianceThreshold(threshold=0)
    selector.fit_transform(data)
    cols=selector.get_support(indices=True)
    data = data.iloc[:,cols]
    # features_anti = features_name[cols]
    # n_features = len(features_anti)
    print("After removing low variance:{}".format(data.shape))
    
    # sm = SMOTE(random_state=42)
    # data, target = sm.fit_resample(data, target)
    # print('Resampled dataset shape %s' % Counter(target))
    
    # Preprocessing - Feature Selection
    std_scaler = MinMaxScaler()
    data = std_scaler.fit_transform(data)
    
    _, pvalue = chi2(data, target)

    threshold = 0.05
    cols_model = np.where(pvalue < threshold)[0]
    # coef_sel = pvalue[cols_model]
    print(len(cols_model))
    
    # if len(cols_model) == 0:
    #     continue
        
    # features_anti = features_anti[cols_model]
    # n_features = len(features_anti)
    data = data[:, cols_model]
    
    # print("After select from model:{}".format(data.shape))
    
    # Initialize Variables:
    scores_auc = np.zeros([NUM_TRIALS,len(classifiers)])
    scores_acc = np.zeros([NUM_TRIALS,len(classifiers)])
    scores_sens = np.zeros([NUM_TRIALS,len(classifiers)])
    scores_spec = np.zeros([NUM_TRIALS,len(classifiers)])
    scores_kappa = np.zeros([NUM_TRIALS,len(classifiers)])
    scores_prec = np.zeros([NUM_TRIALS,len(classifiers)])

    # Loop for each trial
    update_progress(0)
    for i in range(NUM_TRIALS):
        #print("Trial = {}".format(i))
    
        inner_cv = StratifiedKFold(n_splits=inner_loop_cv, shuffle=True, random_state=i)
        outer_cv = StratifiedKFold(n_splits=outer_loop_cv, shuffle=True, random_state=i)
    
        k = 0
    
        for name, clf in zip(names, classifiers):
            model = Pipeline([('clf', clf)])
            
            if name == "KNN":
                grid = K_grid
            elif name == "RBF SVM":
                grid = SVC_grid              
            elif name == "Random Forest" or name == "AdaBoost" or name == "Extra Trees":
                grid = est_grid
            elif name == "Neural Net":
                grid = MLP_grid
            elif name == "Linear SVM":
                grid = C_grid
            elif name == "Decision Tree":
                grid = DT_grid
            elif name == "XGBoost":
                grid = XGBoost_grid
            else:
                grid = C_grid

            # Inner Search
            classif = GridSearchCV(estimator=model, param_grid=grid, cv=inner_cv)
            classif.fit(data, target)
            print(name + ": ")
            print(classif.best_params_)
        
            # Outer Search
            cv_results = cross_validate(classif, data, target, scoring=scoring, cv=outer_cv, return_estimator=True)

            tp = cv_results['test_tp']
            tn = cv_results['test_tn']
            fp = cv_results['test_fp']
            fn = cv_results['test_fn']
            
            sens = np.zeros(outer_loop_cv)
            spec = np.zeros(outer_loop_cv)
            prec = np.zeros(outer_loop_cv)
            
            for j in range(outer_loop_cv):
                TP = tp[j]
                TN = tn[j]
                FP = fp[j]
                FN = fn[j]
                
                # Sensitivity, hit rate, recall, or true positive rate
                sens[j] = TP/(TP+FN)
                
                # Fall out or false positive rate
                FPR = FP/(FP+TN)
                spec[j] = 1 - FPR
                if TP + FP > 0:
                    prec[j] = TP / (TP + FP)

            scores_sens[i,k] = sens.mean()
            scores_spec[i,k] = spec.mean()
            scores_prec[i,k] = prec.mean()
            scores_auc[i,k] = cv_results['test_auc'].mean()
            scores_acc[i,k] = cv_results['test_acc'].mean()
            scores_kappa[i,k] = cv_results['test_kappa'].mean()
            
            k = k + 1
            
        update_progress((i+1)/NUM_TRIALS)

    results = np.zeros((12,len(classifiers)))
    scores = [scores_auc, scores_acc, scores_sens, scores_spec, scores_kappa, scores_prec]
    for counter_scr, scr in enumerate(scores):
        results[2*counter_scr,:] = np.mean(scr,axis=0)
        results[2*counter_scr + 1,:] = np.std(scr,axis=0)
        
    names_scr = ["AUC_Mean", "AUC_Std", "Acc_Mean", "Acc_Std", 
        "Sens_Mean", "Sens_Std", "Spec_Mean", "Spec_Std", 
        "Kappa_Mean", "Kappa_Std", "Prec_Mean", "Prec_Std"]
        
    directory = "./results/"+name_dataset
    if not os.path.exists(directory):
        os.makedirs(directory)

    results_df=pd.DataFrame(results, columns=names, index=names_scr)
    results_df.to_csv(directory+"/SMOTE_results_"+name_dataset+".csv")

    df_auc = pd.DataFrame(scores_auc, columns=names)
    df_auc.to_csv(directory+"/SMOTE_"+name_dataset+"_auc.csv")
    
    df_acc = pd.DataFrame(scores_acc, columns=names)
    df_acc.to_csv(directory+"/SMOTE_"+name_dataset+"_acc.csv")
    
    df_sens = pd.DataFrame(scores_sens, columns=names)
    df_sens.to_csv(directory+"/SMOTE_"+name_dataset+"_sens.csv")
    
    df_spec = pd.DataFrame(scores_spec, columns=names)
    df_spec.to_csv(directory+"/SMOTE_"+name_dataset+"_spec.csv")
    
    df_kappa = pd.DataFrame(scores_kappa, columns=names)
    df_kappa.to_csv(directory+"/SMOTE_"+name_dataset+"_kappa.csv")
    
    df_prec = pd.DataFrame(scores_prec, columns=names)
    df_prec.to_csv(directory+"/SMOTE_"+name_dataset+"_prec.csv")
    
    # df_features = pd.DataFrame(coef_sel, columns = ["p-value chi2"], index=features_anti)
    # df_features.to_csv(directory+"/features_"+name_dataset+".csv")




if __name__ == "__main__":   
    
    # ignore all future warnings
    simplefilter(action='ignore', category= FutureWarning)
    simplefilter(action='ignore', category= UserWarning)
    simplefilter(action='ignore', category= DeprecationWarning)

    # read raw data and clean the data by replacing NAN with 0
    data_raw = pd.read_csv('../raw_data/final_inputed_results_withChemicals.csv', header=[0])
    data_raw = data_raw.fillna(0)

    # read group information from processed raw data
    group = data_raw["group"]

    # read demography information from processed raw data
    demography = data_raw[["gender", "age", "BMI", "education"]]

    # read dietary and behavioral habits information from processed raw data
    dietary_and_behavioral_habits = data_raw[["smoking_history", "drinking_habit", "green_tea_habit", "regular_staple_food", 
                                        "frenquency_on_pickled_food", "frenquency_on_fresh_fruit_vegetable", 
                                        "preference_of_hot_food_and_water", "preference_of_hard_food", "daily_meal_times"]]

    # read clinic biochemical index information from processed raw data
    clinic_biochemical_index = data_raw[["LYM_pcent", "MON_pcent", "NEU_pcent", "GLB", "ALT_/AST", "ALP", "γ_GGT", "UREA", "CRE"]]

    # read chemical exposome from processed raw data
    chemical_exposome = data_raw.iloc[:,28:] # from column AC (namely Dodecanoic acid) to the end

    # read medication history from processed raw data
    medication_history = data_raw[["medication_history"]]
    
    # try ml models
    mlTest(chemical_exposome, group, "chemical_exposome")
