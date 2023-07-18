import warnings 
warnings.filterwarnings('ignore')
# Libraries
#Basic libraries
import os
import numpy as np
import pandas as pd
from pathlib import Path

from sklearn.model_selection import GroupShuffleSplit, RepeatedStratifiedKFold, LeaveOneOut, RepeatedKFold
from sklearn.base import clone
from scipy.stats import mannwhitneyu, spearmanr
from stabl.visualization import scatterplot_features, boxplot_features
from stabl.stabl import Stabl, save_stabl_results

from sklearn.linear_model import Lasso, LassoCV, LogisticRegressionCV, LogisticRegression, LinearRegression, ElasticNetCV, Lasso

#STABL pipelines
from stabl.multi_omic_pipelines import multi_omic_stabl, multi_omic_stabl_cv
from stabl.single_omic_pipelines import single_omic_stabl, single_omic_stabl_cv

#Preprocessing functions
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from stabl.preprocessing import LowInfoFilter, remove_low_info_samples

# Import Data
X_noEGA_pen = pd.read_csv('../Onset of Labor csv/immunome_noEGA_pen_OOL.csv',index_col="ID")
X_noEGA = pd.read_csv('../Onset of Labor csv/immunome_noEGA_OOL.csv',index_col="ID")

X = X_noEGA_pen

EGA_error = pd.read_csv('../Onset of Labor csv/EGA_error.csv',index_col="ID").iloc[:,0]
# Preprocessing
remove_low_info_samples(X, threshold=1.)
preprocessing = Pipeline(
	steps=[
		('lif', LowInfoFilter(0.2)),
		('variance', VarianceThreshold(0.0)),
		('impute', SimpleImputer(strategy='median'))
	])

X = pd.DataFrame(
	data=preprocessing.fit_transform(X),
	index=X.index,
	columns=preprocessing.get_feature_names_out()
)

# Models
## single_omic_stabl_pipeline personalized

from sklearn.svm import l1_min_c

from stabl.metrics import jaccard_matrix
from stabl.utils import compute_CI, permutation_test_between_clfs
from stabl.preprocessing import LowInfoFilter, remove_low_info_samples

from scipy import stats
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score, average_precision_score, r2_score, mean_squared_error, mean_absolute_error

from stabl.pipelines_utils import save_plots, compute_scores_table

lasso = Lasso(max_iter=int(1e6))
lasso_cv = LassoCV(n_alphas=50, max_iter=int(1e6), n_jobs=-1)
en_cv = ElasticNetCV(n_alphas=50, max_iter=int(1e6), n_jobs=-1, l1_ratio=.5)

logit_lasso_cv = LogisticRegressionCV(penalty="l1", solver="liblinear", Cs=np.logspace(-2, 2, 50),
                                      max_iter=int(1e6), class_weight="balanced", scoring="roc_auc",
                                      n_jobs=-1
                                      )

logit_en_cv = LogisticRegressionCV(penalty="elasticnet", solver="saga", Cs=np.logspace(-2, 2, 50),
                                   max_iter=int(1e6), class_weight="balanced", scoring="roc_auc",
                                   n_jobs=-1, l1_ratios=[.5]
                                   )

logit = LogisticRegression(penalty=None, class_weight="balanced", max_iter=int(1e6))
linreg = LinearRegression()

def compute_scores_table_without_STABL(
        predictions_dict,
        y,
        task_type="binary",
        selected_features_dict=None
):
    scores_columns = []
    if selected_features_dict is not None:
        if task_type == "binary":
            scores_columns = ["ROC AUC", "Average Precision", "N features", "CVS"]

        elif task_type == "regression":
            scores_columns = ["R2", "RMSE", "MAE", "N features", "CVS"]
            
    else:
        if task_type == "binary":
            scores_columns = ["ROC AUC", "Average Precision"]

        elif task_type == "regression":
            scores_columns = ["R2", "RMSE", "MAE"]

    table_of_scores = pd.DataFrame(data=None, columns=scores_columns)

    for model, preds in predictions_dict.items():

        if task_type == "binary":
            model_roc = roc_auc_score(y, preds)
            model_roc_CI = compute_CI(y, preds, scoring="roc_auc")
            cell_value = f"{model_roc:.3f} [{model_roc_CI[0]:.3f}, {model_roc_CI[1]:.3f}]"
            table_of_scores.loc[model, "ROC AUC"] = cell_value

            model_ap = average_precision_score(y, preds)
            model_ap_CI = compute_CI(y, preds, scoring="average_precision")
            cell_value = f"{model_ap:.3f} [{model_ap_CI[0]:.3f}, {model_ap_CI[1]:.3f}]"
            table_of_scores.loc[model, "Average Precision"] = cell_value

        elif task_type == "regression":
            model_r2 = r2_score(y, preds)
            model_r2_CI = compute_CI(y, preds, scoring="r2")
            table_of_scores.loc[model, "R2"] = f"{model_r2:.3f} [{model_r2_CI[0]:.3f}, {model_r2_CI[1]:.3f}]"

            model_rmse = np.sqrt(mean_squared_error(y, preds))
            model_rmse_CI = compute_CI(y, preds, scoring="rmse")
            table_of_scores.loc[model, "RMSE"] = f"{model_rmse:.3f} [{model_rmse_CI[0]:.3f}, {model_rmse_CI[1]:.3f}]"

            model_mae = mean_absolute_error(y, preds)
            model_mae_CI = compute_CI(y, preds, scoring="mae")
            table_of_scores.loc[model, "MAE"] = f"{model_mae:.3f} [{model_mae_CI[0]:.3f}, {model_mae_CI[1]:.3f}]"

        if selected_features_dict is not None:

            if model != "STABL":
                sel_features = selected_features_dict[model]["Fold nb of features"]
                jaccard_mat = jaccard_matrix(selected_features_dict[model]["Fold selected features"], remove_diag=False)
                jaccard_val = jaccard_mat[np.triu_indices_from(jaccard_mat, k=1)]

                median_features = np.median(sel_features)
                iqr_features = np.quantile(sel_features, [.25, .75])
                cell_value = f"{median_features:.3f} [{iqr_features[0]:.3f}, {iqr_features[1]:.3f}]"
                table_of_scores.loc[model, "N features"] = cell_value

                jaccard_median = np.median(jaccard_val)
                jaccard_iqr = np.quantile(jaccard_val, [.25, .75])
                cell_value = f"{jaccard_median:.3f} [{jaccard_iqr[0]:.3f}, {jaccard_iqr[1]:.3f}]"
                table_of_scores.loc[model, "CVS"] = cell_value

    return table_of_scores

def single_omic_cv(
        X,
        y,
        outer_splitter,
        task_type,
        save_path,
        outer_groups=None
):
    
    models = ["Lasso", "Lasso 1SE", "ElasticNet"]

    os.makedirs(Path(save_path, "Training CV"), exist_ok=True)
    os.makedirs(Path(save_path, "Summary"), exist_ok=True)

    # Initializing the df containing the data of all omics
    predictions_dict = dict()
    selected_features_dict = dict()
    coefficients_dict = dict()

    for model in models:
        predictions_dict[model] = pd.DataFrame(data=None, index=y.index)
        selected_features_dict[model] = []
        coefficients_dict[model] = pd.DataFrame(data=None, index=X.columns)

    i = 1
    for train, test in outer_splitter.split(X, y, groups=outer_groups):
        # Jonas additional code in case outer_splitter is LeaveOneOut
        if isinstance(outer_splitter, LeaveOneOut):
            print(f" Iteration {i} over {X.shape[0]} ".center(80, '*'), "\n")
        else:
            print(f" Iteration {i} over {outer_splitter.get_n_splits()} ".center(80, '*'), "\n")
        # end additional code
        train_idx, test_idx = y.iloc[train].index, y.iloc[test].index

        fold_selected_features = dict()
        for model in models:
            fold_selected_features[model] = []

        print(f"{len(train_idx)} train samples, {len(test_idx)} test samples")

        # __other models__
        y_train, y_test = y.loc[train_idx], y.loc[test_idx]
        X_train = X.loc[train_idx]
        X_test = X.loc[test_idx]

        # __Lasso__
        if task_type == "binary":
            inner_splitter = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=42)
            model = clone(logit_lasso_cv).set_params(cv=inner_splitter)
            predictions = model.fit(X_train, y_train).predict_proba(X_test)[:, 1]
        else:
            inner_splitter = RepeatedKFold(n_splits=5, n_repeats=5, random_state=42)
            model = clone(lasso_cv).set_params(cv=inner_splitter)
            predictions = model.fit(X_train, y_train).predict(X_test)

        selected_features_dict["Lasso"].append(list(X_train.columns[np.where(model.coef_.flatten())]))
        predictions_dict["Lasso"].loc[test_idx, f"Fold n°{i}"] = predictions
        coefficients_dict["Lasso"][f"Fold n°{i}"] = model.coef_

        # __Lasso 1SE__
        if task_type == "binary":
            # Jonas additional code
            new_best_c_corr = model.C_[0] - model.scores_[True].std() / np.sqrt(inner_splitter.get_n_splits())
            if new_best_c_corr < 0:
                best_c_corr = abs(model.C_[0])
            else:
                best_c_corr = new_best_c_corr
            # end of new code
            model = LogisticRegression(penalty='l1', solver='liblinear', C=best_c_corr, class_weight='balanced',
                                       max_iter=2_000_000)
            predictions = model.fit(X_train, y_train).predict_proba(X_test)[:, 1]

        selected_features_dict["Lasso 1SE"].append(list(X_train.columns[np.where(model.coef_.flatten())]))
        predictions_dict["Lasso 1SE"].loc[test_idx, f"Fold n°{i}"] = predictions
        coefficients_dict["Lasso 1SE"][f"Fold n°{i}"] = model.coef_

        # __EN__
        if task_type == "binary":
            model = clone(logit_en_cv).set_params(cv=inner_splitter)
            predictions = model.fit(X_train, y_train).predict_proba(X_test)[:, 1]

        else:
            model = clone(en_cv).set_params(cv=inner_splitter)
            predictions = model.fit(X_train, y_train).predict(X_test)

        selected_features_dict["ElasticNet"].append(list(X_train.columns[np.where(model.coef_.flatten())]))
        predictions_dict["ElasticNet"].loc[test_idx, f"Fold n°{i}"] = predictions
        coefficients_dict["ElasticNet"][f"Fold n°{i}"] = model.coef_

        i += 1

    # __SAVING_RESULTS__

    if y.name is None:
        y.name = "outcome"

    summary_res_path = Path(save_path, "Summary")
    cv_res_path = Path(save_path, "Training CV")

    jaccard_matrix_dict = dict()
    formatted_features_dict = dict()

    for model in models:

        jaccard_matrix_dict[model] = jaccard_matrix(selected_features_dict[model])
        
        # Jonas additional code in case outer_splitter is LeaveOneOut
        if isinstance(outer_splitter, LeaveOneOut):
            index=[f"Fold {i}" for i in range(X.shape[0])]
        else:
            index=[f"Fold {i}" for i in range(outer_splitter.get_n_splits())]
        # end additional code

        formatted_features_dict[model] = pd.DataFrame(
            data={
                "Fold selected features": selected_features_dict[model],
                "Fold nb of features": [len(el) for el in selected_features_dict[model]]
            },
            index=index # Jonas'additional code linked to this parameter
        )
        formatted_features_dict[model].to_csv(Path(cv_res_path, f"Selected Features {model}.csv"))
        coefficients_dict[model].to_csv(Path(cv_res_path, f"{model} coefficients.csv"))

    predictions_dict = {model: predictions_dict[model].median(axis=1) for model in predictions_dict.keys()}

    table_of_scores = compute_scores_table_without_STABL(
        predictions_dict=predictions_dict,
        y=y,
        task_type=task_type,
        selected_features_dict=formatted_features_dict
    )

    table_of_scores.to_csv(Path(summary_res_path, "Scores training CV.csv"))
    table_of_scores.to_csv(Path(cv_res_path, "Scores training CV.csv"))

    save_plots(
        predictions_dict=predictions_dict,
        y=y,
        task_type=task_type,
        save_path=cv_res_path
    )

    return predictions_dict

single_omic_cv(
        X=X,
        y=EGA_error,
        outer_splitter=RepeatedStratifiedKFold(n_repeats=10, n_splits=5, random_state=42),
        task_type='regression',
        save_path=f"../Results_EGA_correction/immunome_noEGA_pen_OOL/Other_models",
        outer_groups=None
)