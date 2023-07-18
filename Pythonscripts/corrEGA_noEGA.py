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

from sklearn.linear_model import LassoCV, LogisticRegressionCV, LogisticRegression, LinearRegression, ElasticNetCV, Lasso

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

X = X_noEGA
data_name = "immunome_noEGA_OOL"

EGA_error = pd.read_csv('../Onset of Labor csv/EGA_error.csv',index_col="ID").iloc[:,0]
# Preprocessing
remove_low_info_samples(X)
preprocessing = Pipeline(
	steps=[
		('lif', LowInfoFilter(0.2)),
		('variance', VarianceThreshold(0.0)),
		('impute', SimpleImputer(strategy='median')),
		('std', StandardScaler())
	])
# Training CV
run_name = "LassoKF_L0.2_V0.1_B0.5"
logit_en = LogisticRegression(penalty='l1', max_iter=int(1e6), solver='liblinear', class_weight='balanced')

stabl = Stabl(base_estimator=clone(logit_en),
	lambda_name='C',
	lambda_grid=list(np.linspace(0.01, 1, 30)),
	n_bootstraps=1000,
	artificial_type='knockoff',
	artificial_proportion=1.,
	sample_fraction=0.5,
	replace=False,
	fdr_threshold_range=list(np.arange(0., 1., .01)),
	sample_weight_bootstrap=None,
	bootstrap_threshold=1e-5,
	backend_multi='threading',
	verbose=0,
	n_jobs=-1,
	random_state=42)

stability_selection = clone(stabl).set_params(artificial_type=None, hard_threshold=0.3)
outer_splitter = RepeatedKFold(n_splits=5, n_repeats=10, random_state=42)

single_omic_stabl_cv(
	X=X,
	y=EGA_error,
	outer_splitter=outer_splitter,
	stabl=stabl,
	stability_selection=stability_selection,
	task_type='regression',
	save_path=f"../Results_EGA_correction/{data_name}/{run_name}"
)
# Univariate
os.makedirs(f"../Results_EGA_correction/{data_name}/{run_name}" + '/Univariate', exist_ok=True)
Spearmancorr = {}
features = X.columns
for feature in features:
	corr, pval = spearmanr(X[feature], EGA_error)
	Spearmancorr[feature] = [corr, pval]

SpearmanPvalue = pd.DataFrame(Spearmancorr).T
SpearmanPvalue.columns = ['Spearman corr', 'pvalue']
SpearmanPvalue.sort_values('pvalue', inplace=True)
SpearmanPvalue.to_csv(f"../Results_EGA_correction/{data_name}/{run_name}"+'/Univariate/SpearmanCorrelationsPval.csv', index=True)

scatterplot_features(
	SpearmanPvalue[:10].index,
	X,
	EGA_error,
	show_fig=False,
	export_file=True,
	path=f"../Results_EGA_correction/{data_name}/{run_name}/Univariate/")
# Final STABL
X_STD = pd.DataFrame(
	data=preprocessing.fit_transform(X),
	index=X.index,
	columns=preprocessing.get_feature_names_out()
)

finalstabl = clone(stabl)
finalstabl.fit(X_STD,EGA_error)

save_stabl_results(finalstabl,f"../Results_EGA_correction/{data_name}/{run_name}"+'/FinalSTABL/',X_STD,EGA_error,task_type='regression')