import os
import numpy as np
import seaborn
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from sklearn.model_selection import train_test_split, GridSearchCV, learning_curve, RepeatedKFold
from sklearn.metrics import mean_squared_error
from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor, ExtraTreesRegressor

from matplotlib import rcParams
plt.rcParams["font.family"] = "arial"
plt.rcParams["font.size"] = 12
plt.rcParams["axes.labelsize"] = 14
plt.rcParams["axes.linewidth"] = 1
plt.rcParams["legend.frameon"] = False
plt.rcParams["legend.framealpha"] = 1.0

def make_dataframe_from_json(jsonfile=None):
	import json
	df = json.load(open(jsonfile, "r"))
	df = pd.DataFrame(df["_default"])
	df = df.T
	df = df.set_index("system")
	df = pd.DataFrame(df, dtype=float)
	return df

def make_dataframe_form_csv(csvfile=None):
	df = pd.read_csv(csvfile)
	if "system" in df.columns:
		df = df.drop("system", axis=1)
	if "Unnamed: 0" in df.columns:
		df = df.drop("Unnamed: 0", axis=1)
	return df

def remove_irregular_samples(df=None):
	before = len(df)
	# remove positive adssorption energy
	df = df[df["E_ads"] < 0.0]

	# remove zero or negative height and width
	for i in df.columns:
		if ("height" in i) or ("width" in i):
			df = df[df[i] > 0.0]

	after = len(df)
	print("removing outliear: {0:d} --> {1:d}".format(before, after))

	return df

def make_shap_plot(model=None, X=None, outdir=None):
	"""
	SHAP value plot.

	Args:
		model:
		X:
		outdir:
	Returns:
	"""
	import shap

	shap_values = shap.TreeExplainer(model).shap_values(X)

	plt.tight_layout()
	shap.summary_plot(shap_values, X_train, plot_type="bar")
	plt.savefig(outdir + "/" + "shap_values.png")
	plt.show()
	plt.close()
	#fig.savefig(outdir + "/" + "shap_bar.png", bbox_inches="tight")

	plt.tight_layout()
	shap.summary_plot(shap_values, X_train)
	plt.savefig(outdir + "/" + "shap_summary.png")
	plt.show()
	plt.close()
	#fig.savefig(outdir + "/" + "shap_values.png", bbox_inches="tight")

def test(df=None, model=None, model_type="lasso"):
	from sklearn.model_selection import cross_validate, RepeatedKFold

	feature_names = df.columns
	cv_model = cross_validate(model, X, y, cv=RepeatedKFold(n_splits=5, n_repeats=4), return_estimator=True)
	if model_type == "lasso":
		coefs = pd.DataFrame([est.coef_ for est in cv_model["estimator"]], columns=feature_names[1:])
	else:  # tree regressor
		coefs = pd.DataFrame([est.feature_importances_ for est in cv_model["estimator"]], columns=feature_names[1:])

	seaborn.boxplot(data=coefs, orient="h", saturation=0.5)
	plt.grid()
	plt.show()

def plot_feature_importance(rf=None, X=None, outdir=None):
	#std = np.std([tree.feature_importances_ for tree in rf.estimators_], axis=0)
	feature_imp = pd.DataFrame({"name": X.columns, "Coef": rf.feature_importances_})

	_, ax = plt.subplots(figsize=(10, 10))
	ax.barh(feature_imp["name"].iloc[::-1], feature_imp["Coef"].iloc[::-1], height=0.6, color="limegreen")

	ax.set_xlabel("Feature importance")
	ax.axvline(x=0, color="black", linewidth=0.5)
	plt.tight_layout()
	plt.savefig(outdir + "/" + "feature_importance.png")
	plt.show()
	plt.close()

# main
outdir = "regression_results"
os.makedirs(outdir, exist_ok=True)
os.system("rm {}/*".format(outdir))

df = make_dataframe_from_json(jsonfile="sample.json")
df = remove_irregular_samples(df)

X = df.drop("E_ads", axis=1)
y = -df["E_ads"]

#cv = 10
#test_size = 1.0/cv
cv = RepeatedKFold(n_splits=5, n_repeats=4)
test_size = 0.2

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=1)

scaler = StandardScaler()
#scaler = MinMaxScaler()
#scaler = RobustScaler()

methods = [Ridge(), Lasso()]
names   = ["ridge", "lasso"]

print("-------- Linear Regression ---------")
lr = LinearRegression()
lr.fit(X_train, y_train)

print(pd.DataFrame({"name": X.columns, "Coef": lr.coef_}).sort_values(by="Coef"))
print("Training set score: {:.3f}".format(lr.score(X_train, y_train)))
print("Test set score: {:.3f}".format(lr.score(X_test, y_test)))
print("RMSE: {:.3f}".format(np.sqrt(mean_squared_error(y_test, lr.predict(X_test)))))

for name, method in zip(names, methods):
	print("-------- %s ---------" % name)

	pipe = make_pipeline(scaler, method)
	param_grid = {name + "__alpha": list(10**np.arange(-3, 3, 0.1))}

	grid = GridSearchCV(pipe, param_grid=param_grid, cv=cv)
	grid.fit(X_train, y_train)

	print(pd.DataFrame({"name": X.columns,
						"Coef": grid.best_estimator_.named_steps[name].coef_}).sort_values(by="Coef"))
	print("Best parameters: {}".format(grid.best_params_))
	print("Training set score: {:.3f}".format(grid.score(X_train, y_train)))
	print("Test set score: {:.3f}".format(grid.score(X_test, y_test)))
	print("RMSE: {:.3f}".format(np.sqrt(mean_squared_error(y_test, grid.predict(X_test)))))

	seaborn.regplot(y=y.values, x=grid.predict(X))
	plt.title("%s regression" % name)
	plt.show()

# plot coefficient of LASSO
lasso_coef = pd.DataFrame({"name": X.columns, "Coef": grid.best_estimator_.named_steps["lasso"].coef_})
fig, ax = plt.subplots(figsize=(10, 10))
ax.barh(lasso_coef["name"].iloc[::-1], lasso_coef["Coef"].iloc[::-1], height=0.6, color="royalblue")
ax.set_xlabel("Coefficient")

ax.axvline(x=0, color="black", linewidth=0.5)
plt.tight_layout()
plt.savefig(outdir + "/" + "lasso_coef.png")
plt.show()
plt.close()

# learning_curve for LASSO
best_param = list(grid.best_params_.values())[0]
pipe = Pipeline([("scl", scaler), ("lasso", Lasso(alpha=best_param))])
train_sizes, train_scores, test_scores = learning_curve(estimator=pipe, X=X_train, y=y_train,
														train_sizes=np.linspace(0.2, 1.0, 10), cv=cv)
train_mean = np.mean(train_scores, axis=1)
train_std  = np.std(train_scores,  axis=1)
test_mean  = np.mean(test_scores,  axis=1)
test_std   = np.std(test_scores,   axis=1)

plt.plot(train_sizes, train_mean, color="blue", marker="o", markersize=5, label="training accuracy")
plt.fill_between(train_sizes, train_mean+train_std, train_mean-train_std, alpha=0.15, color="blue")

plt.plot(train_sizes, test_mean, color="green", marker="s", markersize=5, label="validation accuracy")
plt.fill_between(train_sizes, test_mean+test_std, test_mean-test_std, alpha=0.15, color="green")

plt.xticks()
plt.yticks()
plt.xlabel("Number of training samples")
plt.ylabel("Accuracy")
#plt.legend(loc="lower right", fontsize=14)
plt.ylim([0.0, 1.0])
plt.savefig(outdir + "/" + "learning_curve_lasso.png")
plt.show()
plt.close()

# Tree regression
print("-------- Tree Regression ---------")
#rf = RandomForestRegressor()
rf = GradientBoostingRegressor()
#rf = ExtraTreesRegressor()

# do grid search to find hyper parameters
param_grid = {"n_estimators": [50, 100, 150], "max_depth": [5, 10, 15]}
grid_tree = GridSearchCV(rf, param_grid=param_grid, cv=cv)
grid_tree.fit(X_train, y_train)

rf = grid_tree.best_estimator_

print("Best parameters: {}".format(grid_tree.best_params_))
print("Training set score: {:.3f}".format(rf.score(X_train, y_train)))
print("Test set score: {:.3f}".format(rf.score(X_test, y_test)))
print("RMSE : {:.3f}".format(np.sqrt(mean_squared_error(y_test, rf.predict(X_test)))))

seaborn.regplot(y=y.values, x=rf.predict(X))
plt.title("Random forest")
plt.show()

# learning_curve for Random forest
train_sizes, train_scores, test_scores = learning_curve(estimator=rf, X=X_train, y=y_train,
														train_sizes=np.linspace(0.2, 1.0, 10), cv=cv)
train_mean = np.mean(train_scores, axis=1)
train_std  = np.std(train_scores, axis=1)
test_mean  = np.mean(test_scores, axis=1)
test_std   = np.std(test_scores, axis=1)

plt.plot(train_sizes, train_mean, color="blue", marker="o", markersize=5, label="training accuracy")
plt.fill_between(train_sizes, train_mean+train_std, train_mean-train_std, alpha=0.15, color="blue")

plt.plot(train_sizes, test_mean, color="green", marker="s", markersize=5, label="validation accuracy")
plt.fill_between(train_sizes, test_mean+test_std, test_mean-test_std, alpha=0.15, color="green")

plt.xticks()
plt.yticks()
plt.xlabel("Number of training samples")
plt.ylabel("Accuracy")
#
#plt.legend(loc="lower right", fontsize=14)
plt.ylim([0.0, 1.0])
plt.savefig(outdir + "/" + "learning_curve_RF.png")
plt.show()
plt.close()

# plot correlation matrix
corr = df.corr()
_, ax = plt.subplots(figsize=(12, 12))
seaborn.heatmap(corr, vmax=1, vmin=-1, center=0, annot=False, annot_kws={"size": 10},
				cbar=True, cmap="RdBu_r", square=True, fmt=".1f", ax=ax)
plt.savefig(outdir + "/" + "correlation.png")
plt.show()
plt.close()

make_shap_plot(model=rf, X=X_train, outdir=outdir)

test(df=df, model=grid.best_estimator_.named_steps["lasso"], model_type="lasso")
test(df=df, model=rf, model_type="rf")
