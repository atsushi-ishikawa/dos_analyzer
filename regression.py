import os
import numpy as np
import seaborn
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import Ridge, Lasso
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from sklearn.model_selection import train_test_split, GridSearchCV, learning_curve, RepeatedKFold
from sklearn.metrics import mean_squared_error
from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor, ExtraTreesRegressor
from xgboost import XGBRegressor
from matplotlib import rcParams
from sklearn.model_selection import cross_val_score

plt.rcParams["font.family"] = "arial"
plt.rcParams["font.size"] = 16
plt.rcParams["axes.labelsize"] = 18
plt.rcParams["axes.linewidth"] = 1
plt.rcParams["legend.frameon"] = False
plt.rcParams["legend.framealpha"] = 1.0
plt.rcParams["axes.axisbelow"] = True

n_splits  = 5  # K-fold CV (default: 5 -> 10)
n_repeats = 1  # number of repeats (default: 10)
random_state = 0
cv = RepeatedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=random_state)

# whether to show figure
showfigure = False

def make_dataframe_from_json(jsonfile=None):
    """
    Get descriptors from json file.

    Args:
        jsonfile:
    Returns:
        df:
    """
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


def remove_irregular_samples(df=None, key="E_ads"):
    from sklearn.impute import SimpleImputer

    # remove positive adssorption energy
    before = len(df)
    df = df[df[key] < 0.0]
    after = len(df)
    print("removing positove adsorption energy: {0:d} --> {1:d}".format(before, after))

    # remove negative height and width
    before = len(df)
    for i in df.columns:
        if ("height" in i) or ("width" in i):
            df = df[df[i] >= 0.0]
    after = len(df)
    print("removing strange peak: {0:d} --> {1:d}".format(before, after))

    imp_mean = SimpleImputer(missing_values=0.0, strategy="median")   # median is better
    df2 = pd.DataFrame(imp_mean.fit_transform(df))
    df2.index = df.index
    df2.columns = df.columns

    return df2


def make_shap_plot(model=None, model_name=None, X=None, outdir=None):
    """
    SHAP (Shapley additive explanation) value plot.

    Args:
        model: regression model
        model_name: name of regression model
        X:
        outdir: directory to save file
    Returns:
    """
    import shap

    shap_values = shap.TreeExplainer(model).shap_values(X)

    # plotting SHAP value
    fig = plt.figure()
    fig.tight_layout()
    shap.summary_plot(shap_values, X_train, plot_type="bar", plot_size=(10, 10), show=showfigure)
    fig.savefig(outdir + "/" + "shap_values_" + model_name + ".png", bbox_inches="tight")
    plt.clf()
    plt.close()

    # SHAP summary plot
    fig = plt.figure()
    fig.tight_layout()
    shap.summary_plot(shap_values, X_train, plot_size=(6, 10), show=showfigure)
    fig.savefig(outdir + "/" + "shap_summary_" + model_name + ".png", bbox_inches="tight")
    plt.clf()
    plt.close()


def plot_variability_of_coefficients(df=None, model=None, model_name="lasso", outdir=None):
    """
    Plot the variability of LASSO or tree-regression coefficients.
    See https://scikit-learn.org/stable/auto_examples/inspection/plot_linear_model_coefficient_interpretation.html

    Args:
        df:
        model: regression model
        model_name: name of regression model
        outdir: directory to save file
    Returns:
        None
    """
    from sklearn.model_selection import cross_validate

    feature_names = df.columns
    if "E_ads" in feature_names:
        feature_names.remove("E_ads")

    if model_name is None:
        raise ValueError("input model name")

    cv_model = cross_validate(model, X, y, cv=cv, return_estimator=True)

    if model_name == "lasso":
        coefs = pd.DataFrame([est.coef_ for est in cv_model["estimator"]], columns=feature_names)
    else:
        coefs = pd.DataFrame([est.feature_importances_ for est in cv_model["estimator"]], columns=feature_names)

    fig, ax = plt.subplots(figsize=(8, 10))
    seaborn.boxplot(data=coefs, orient="h", saturation=0.5, color="cyan", linewidth=1.0)
    ax.set_xlabel("Coefficient")
    ax.axvline(x=0, color="black", linewidth=0.5)
    ax.yaxis.grid(color="gray", linewidth=0.25)
    fig.tight_layout()
    plt.savefig(outdir + "/" + "coef_variability_" + model_name + ".png")
    if showfigure:
        plt.show()
    plt.clf()
    plt.close()


def plot_feature_importance(model=None, X=None, outdir=None):
    """
    Plot feature importance of tree regressors.

    Args:
        model: regression model
        X:
        outdir: directory to save file
    Returns:
        None
    """
    feature_imp = pd.DataFrame({"name": X.columns, "Coef": model.feature_importances_})

    _, ax = plt.subplots(figsize=(10, 10))
    ax.barh(feature_imp["name"].iloc[::-1], feature_imp["Coef"].iloc[::-1], height=0.6, color="limegreen")
    ax.set_xlabel("Feature importance")
    ax.axvline(x=0, color="black", linewidth=0.5)
    plt.tight_layout()
    plt.savefig(outdir + "/" + "feature_importance.png")
    if showfigure:
        plt.show()
    plt.clf()
    plt.close()


def plot_correlation_matrix(df=None, outdir=None):
    """
    Plot correlation matrix among descriptors.

    Args:
        df: DataFrame
        outdir: directory to save file
    Returns:
        None
    """
    corr = df.corr()
    _, ax = plt.subplots(figsize=(14, 14))
    seaborn.heatmap(corr, vmax=1, vmin=-1, center=0, annot=False,
                    cbar=True, cmap="RdBu_r", square=True, fmt=".1f", ax=ax)
    plt.tight_layout()
    plt.savefig(outdir + "/" + "correlation.png")
    if showfigure:
        plt.show()
    plt.clf()
    plt.close()


def plot_scatter_and_line(x=None, y=None, model_name=None, outdir=None):
    """
    Plot the scattered points and regression line.

    Args:
        x: X
        y: y
        model_name:
        outdir: directory to save file
    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=(6, 6))
    seaborn.regplot(y=y, x=x, scatter_kws={"color": "navy", 'alpha': 0.3}, line_kws={"color": "navy"})
    ax.set_xlabel("Predicted value")
    ax.set_ylabel("True value")
    fig.tight_layout()
    plt.savefig(outdir + "/" + "regplot_" + model_name + ".png")
    if showfigure:
        plt.show()
    plt.clf()
    plt.close()


# main starts
outdir = "regression_results"
os.makedirs(outdir, exist_ok=True)
os.system("rm {}/*".format(outdir))

# setup dataframe
df = make_dataframe_from_json(jsonfile="surf_x.json")
X  = df

# adsorption energy
df_Eads = make_dataframe_from_json(jsonfile="Eads.json")

df_all = pd.DataFrame()

adsorbates = ["CO", "CH3", "NO"]
adsorbate_desciptor = "p_center"
#adsorbate_desciptor = "e_fermi"

for adsorbate in adsorbates:
    df_ads = make_dataframe_from_json(jsonfile=adsorbate + "_x" + ".json")
    for i, name in enumerate(df.index):
        label = name + "_" + adsorbate

        ads_value = df_ads.loc[adsorbate][adsorbate_desciptor]
        E_ads = df_Eads.loc[label]["E_ads"]

        tmp = pd.Series(X.loc[name])
        series_ads   = pd.Series([ads_value], index=["ads_descriptor"])
        series_name  = pd.Series([label], index=["system"])
        series_Eads  = pd.Series([E_ads], index=["E_ads"])
        tmp = pd.concat([tmp, series_ads])
        tmp = pd.concat([tmp, series_name])
        tmp = pd.concat([tmp, series_Eads])

        df_all = pd.concat([df_all, tmp], axis=1)

df_all = df_all.T
df_all = df_all.set_index("system")
df_all = remove_irregular_samples(df_all)

X = df_all.copy()
y = df_all["E_ads"]  # more negative = stronger adsorption
#y = -df_all["E_ads"]  # more positive = stronger adsorption

del X["E_ads"]

# plot correlation matrix
plot_correlation_matrix(df=df_all, outdir=outdir)

# train-test split
test_size = 1.0 / n_splits
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=random_state)

scaler = StandardScaler()
#scaler = MinMaxScaler()
#scaler = RobustScaler()

#
# ridge and lasso
#
max_iter = 3000
methods = [Ridge(max_iter=max_iter), Lasso(max_iter=max_iter)]
names = ["ridge", "lasso"]
for name, method in zip(names, methods):
    print("----- %s -----" % name)

    pipe = make_pipeline(scaler, method)
    param_grid = {name + "__alpha": list(10**np.arange(-3, 3, 0.1))}

    grid = GridSearchCV(pipe, param_grid=param_grid, cv=cv)
    grid.fit(X_train, y_train)
    cv_score = cross_val_score(grid, X_train, y_train, cv=cv, scoring="neg_root_mean_squared_error")

    print("Best parameters: {}".format(grid.best_params_))
    print("Score(train): {:.3f}".format(grid.score(X_train, y_train)))
    print("Score(test) : {:.3f}".format(grid.score(X_test, y_test)))
    print("RMSE(train) : {:.3f}".format(np.sqrt(mean_squared_error(y_train, grid.predict(X_train)))))
    print("RMSE(CV)    : {0:.3f}, std: {1:.3f}".format(-cv_score.mean(), cv_score.std()))
    print("RMSE(test)  : {:.3f}".format(np.sqrt(mean_squared_error(y_test, grid.predict(X_test)))))

    plot_scatter_and_line(x=grid.predict(X), y=y.values, model_name=name, outdir=outdir)

# plot coefficient of LASSO
lasso_coef = pd.DataFrame({"name": X.columns, "Coef": grid.best_estimator_.named_steps["lasso"].coef_})

fig, ax = plt.subplots(figsize=(8, 10))
ax.barh(lasso_coef["name"].iloc[::-1], lasso_coef["Coef"].iloc[::-1], height=0.6, color="royalblue")
ax.set_xlabel("Coefficient")
ax.axvline(x=0, color="black", linewidth=0.5)
ax.yaxis.grid(color="gray", linewidth=0.25)
plt.tight_layout()
plt.savefig(outdir + "/" + "lasso_coef.png")
if showfigure:
    plt.show()
plt.clf()
plt.close()

lasso = grid.best_estimator_.named_steps["lasso"]
plot_variability_of_coefficients(df=X, model=lasso, model_name="lasso", outdir=outdir)

# learning_curve for LASSO
best_param = list(grid.best_params_.values())[0]
pipe = Pipeline([("scl", scaler), ("lasso", Lasso(alpha=best_param))])
train_sizes, train_scores, val_scores = learning_curve(estimator=pipe, X=X_train, y=y_train, scoring="r2",
                                                       train_sizes=np.linspace(0.2, 1.0, 10), cv=cv)
train_mean = np.mean(train_scores, axis=1)
train_std  = np.std(train_scores, axis=1)
val_mean   = np.mean(val_scores, axis=1)
val_std    = np.std(val_scores, axis=1)

plt.plot(train_sizes, train_mean, color="blue", marker="o", markersize=5, label="training accuracy")
plt.fill_between(train_sizes, train_mean+train_std, train_mean-train_std, alpha=0.15, color="blue")
plt.plot(train_sizes, val_mean, color="green", marker="s", markersize=5, label="validation accuracy")
plt.fill_between(train_sizes, val_mean+val_std, val_mean-val_std, alpha=0.15, color="green")

plt.xticks()
plt.yticks()
plt.xlabel("Number of training samples")
plt.ylabel("Accuracy")
plt.ylim([0.0, 1.0])
plt.tight_layout()
plt.savefig(outdir + "/" + "learning_curve_lasso.png")
if showfigure:
    plt.show()
plt.clf()
plt.close()
#
# Tree regression
#
print("===== Tree Regression =====")

#names = ["randomforest", "gradientboosting", "extratree", "xgb"]
#methods = [RandomForestRegressor(), GradientBoostingRegressor(), ExtraTreesRegressor(), XGBRegressor()]
names = ["randomforest", "gradientboosting", "extratree"]
methods = [RandomForestRegressor(), GradientBoostingRegressor(), ExtraTreesRegressor()]
for name, method in zip(names, methods):
    print("----- %s -----" % name)

    # do grid search to find hyper parameters
    param_grid = {"n_estimators": [50, 100, 150], "max_depth": [1, 2, 3, 4, 5]}
    grid_tree = GridSearchCV(method, param_grid=param_grid, cv=cv)
    grid_tree.fit(X_train, y_train)

    method = grid_tree.best_estimator_

    cv_score = cross_val_score(method, X_train, y_train, cv=cv, scoring="neg_root_mean_squared_error")
    print("Best parameters: {}".format(grid_tree.best_params_))
    print("Score(train): {:.3f}".format(method.score(X_train, y_train)))
    print("Score(test) : {:.3f}".format(method.score(X_test, y_test)))
    print("RMSE(train) : {:.3f}".format(np.sqrt(mean_squared_error(y_train, method.predict(X_train)))))
    print("RMSE(CV)    : {0:.3f}, std: {1:.3f}".format(-cv_score.mean(), cv_score.std()))
    print("RMSE(test)  : {:.3f}".format(np.sqrt(mean_squared_error(y_test, method.predict(X_test)))))

    # make regplot
    plot_scatter_and_line(x=method.predict(X), y=y.values, model_name=name, outdir=outdir)

    # learning_curve for Random forest
    train_sizes, train_scores, val_scores = learning_curve(estimator=method, X=X_train, y=y_train, scoring="r2",
                                                           train_sizes=np.linspace(0.2, 1.0, 10), cv=cv)
    train_mean = np.mean(train_scores, axis=1)
    train_std = np.std(train_scores, axis=1)
    val_mean = np.mean(val_scores, axis=1)
    val_std = np.std(val_scores, axis=1)

    plt.plot(train_sizes, train_mean, color="blue", marker="o", markersize=5, label="training accuracy")
    plt.fill_between(train_sizes, train_mean+train_std, train_mean-train_std, alpha=0.15, color="blue")
    plt.plot(train_sizes, val_mean, color="green", marker="s", markersize=5, label="validation accuracy")
    plt.fill_between(train_sizes, val_mean+val_std, val_mean-val_std, alpha=0.15, color="green")

    plt.xticks()
    plt.yticks()
    plt.xlabel("Number of training samples")
    plt.ylabel("Accuracy")
    plt.ylim([0.0, 1.0])
    plt.savefig(outdir + "/" + "learning_curve_" + name + ".png")
    if showfigure:
        plt.show()
    plt.clf()
    plt.close()

    plot_variability_of_coefficients(df=X, model=method, model_name=name, outdir=outdir)

    # shap plot
    make_shap_plot(model=method, model_name=name, X=X_train, outdir=outdir)
