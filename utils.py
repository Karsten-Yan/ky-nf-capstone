# model selector class to choose models

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import math as math
import scikitplot as skplt

from sklearn.model_selection import cross_val_score, cross_validate
from sklearn.metrics import classification_report
from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.neural_network import MLPClassifier
from xgboost import XGBClassifier

from imblearn.pipeline import Pipeline
from sklearn.utils.class_weight import compute_class_weight

seed = 7

class ModelSelector():

    def __init__(self, preprocessor,cv_state =True, sampler=None, instance_name=None):
        # initialise class with predetrmined models and tests, maybe add functions to add/change them afterwards?
        # add option to use a sampler (default = none)
        self.results = []
        self.name = []
        self.preprocessor = preprocessor
        self.sampler = sampler
        self.cv = []
        self.labels = []
        self.cv_state = cv_state
        self.scoring = {'precision':'precision',
                        'recall': 'recall',
                        'f1score':'f1',
                        'roc auc':'roc_auc'}
        self.models = [
                        ('LR', LogisticRegression(random_state=seed, max_iter = 1000, n_jobs=-1)),
                        ('KNN', KNeighborsClassifier(n_jobs = -1)),
                        ('RF', RandomForestClassifier(random_state=seed, n_jobs = -1)),
                        ('ADAB', AdaBoostClassifier(random_state=seed)),
                        ('XGB',  XGBClassifier(random_state=seed, n_jobs = -1)),
                        ('SVC', SVC(random_state=seed)),
                        ('GNB', GaussianNB()),
                        ('MLP', MLPClassifier(random_state=seed))
                    ]
        self.instance_name = instance_name
        self.tests = ["test_"+elm for elm in self.scoring.keys()]
        self.best_results= []


    # function for determining the cvs, if True do cv for each ref_num
    def __cv__(self, X):
        if self.cv_state == True:
            self.labels = X.ref_number.values
            self.cv = [(np.where(self.labels != label)[0], np.where(self.labels == label)[0]) for label in np.unique(self.labels)]
        else:
            self.cv = 5


    # fit each model and do cv and record cvs and names
    def select_model(self, X, y, state):
        self.cv_state = state
        self.__cv__(X)
        for name, model in self.models:
            model_pipe = Pipeline([(self.preprocessor.__class__.__name__,self.preprocessor),
                                   (self.sampler.__class__.__name__,self.sampler),("name",model)])
            cv_results = cross_validate(model_pipe, X, y, cv=self.cv, scoring=self.scoring, n_jobs=-1, return_train_score=True)
            self.results.append(cv_results)
            self.name.append(name)

    def plot_selection(self):
        # plot the test scores for each test and each model
        plt.figure(figsize = (15,8))
        plt.suptitle('Algorithm Comparison for ' + self.instance_name)
        for idx, test in enumerate(self.tests):
            temp_results = [self.results[i][test] for i in range(len(self.results))]
            plt.subplot(1,len(self.tests),idx+1)
            sns.boxplot(y = temp_results, x = self.name)
            plt.title(test)
            plt.xticks(rotation = 90)
            plt.ylim(0,1.05)
        plt.tight_layout()
        plt.subplots_adjust(top=0.85)

    def get_scores(self, top = 3):
        # get test scores for top 3 tests for each measure and return them as a dataframe
        for idx, test in enumerate(self.tests):
            temp_results = [np.mean(self.results[i][test]) for i in range(len(self.results))]
            test_name = [test for i in range(len(self.name))]
            instance_name = [self.instance_name for i in range (len(self.name))]
            temp_best = sorted(zip(temp_results,self.name,test_name, instance_name), reverse = True)[:top]
            self.best_results.append(temp_best)

        test= pd.DataFrame()
        for i in range(np.array(self.best_results).shape[0]):
            test= pd.concat([test, pd.DataFrame(np.array(self.best_results)[i])], axis=0)
        test.columns = ["test_result","model","test", "instance_name"]
        test.test_result = test.test_result.astype(float)
        test.set_index(["instance_name","test","model"], inplace = True)
        return(test)


def modified_ratio(X, y, threshold = 0):
    # for every position in the predictions calculate a ratio of how myn positions were modified
    # print out those over threshold and return them as a dataframe
    modified_pos = []
    modified_pos_ratio = []
    y_pred_ref_pos = X.ref_pos.values
    ratio_df = pd.DataFrame(pd.concat([pd.Series(y),pd.Series(y_pred_ref_pos)],axis=1))
    ratio_df.columns =["predicted","ref_pos"]
    for elm in ratio_df.ref_pos.unique():
        msk = ratio_df[ratio_df.ref_pos == elm]
        ratio = round((len(msk[msk.predicted == 1])/len(msk))*100,2)
        if ratio > threshold:
            print("Ratio of modified Reads for ref_pos " + str(elm) + " is :" + str(ratio))
            modified_pos.append(elm)
            modified_pos_ratio.append(ratio)
    df = pd.DataFrame({"ratio":modified_pos_ratio}, index = modified_pos )
    df.index.rename("ref_pos",inplace=True)
    return (df)


def feature_importances (model, preprocessor, X, y, debug = False, custom_feature_state = False, custom_feature_list=None):
    # plot feature importances with predefined model and preprocessor, get feature names either from input dataframe
    # if all are used or define a custom one if subset is used
    if custom_feature_state:
        feature_names = custom_feature_list
    else:
        cat_features = list(X.columns[X.dtypes ==  "category"])
        numeric_features = list(X.columns[X.dtypes ==  "float"])
        cat = None
        if "base_1" in X.columns:
            bases = ["_A","_C","_T","_G"]
            cat = [elm + base for elm in cat_features for base in bases ]
        feature_names = numeric_features + cat

    X_trans = preprocessor.fit_transform(X)

    if debug:
        print(feature_names)
        print(len(feature_names))
    model.fit(X_trans, y)

    skplt.estimators.plot_feature_importances(model, feature_names=feature_names, max_num_features=10)
    plt.xticks(rotation=90);



def evaluation (model,X, y, thresh = 30):
    # for model produce confusion matrix, classification report andmodified ratio
    y_pred=model.predict(X)
    skplt.metrics.plot_confusion_matrix(y, y_pred, figsize=(10,10), text_fontsize=20)
    print("-"*30)
    print(classification_report(y,y_pred))
    print("-"*30)
    ratio = modified_ratio(X,y_pred, thresh)
    print("-"*30)
    return ratio


def yeast_prediction(preprocessor, X, y, yeast_df, model, model_name):
    #use test data to fit inputted model and predict the yeast data, bar graph for ratio of wt and ko in predicted
    # states, to see if either condition has more modified than the other
    X_trans = preprocessor.fit_transform(X)
    X_yeast_trans = preprocessor.fit_transform(yeast_df)
    model.fit(X_trans, y)
    y_yeast_pred = model.predict(X_yeast_trans)
    temp_df = pd.DataFrame(yeast_df["file_type"])
    temp_df["modified_status"] = y_yeast_pred
    yeast_group = temp_df.file_type.groupby(temp_df["modified_status"]
                                           ).value_counts(normalize=True).rename("wt_ko_ratio").reset_index()

    sns.barplot(data = yeast_group, x = "modified_status",y = "wt_ko_ratio", hue = "file_type"
            ,edgecolor="grey" , linewidth = 2.5 );
    yeast_group.sort_values(["modified_status","file_type"],inplace=True)
    yeast_group.set_index(["modified_status","file_type"], inplace=True)
    yeast_group = pd.concat([yeast_group], keys=[model_name], names=['Model'])
    return yeast_group
