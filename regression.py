import seaborn
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import mean_squared_error
from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor, ExtraTreesRegressor

df = pd.read_csv("tmpout.csv")
if "system" in df.columns:
	df = df.drop("system", axis=1)
if "Unnamed: 0" in df.columns:
	df = df.drop("Unnamed: 0", axis=1)

#vars = ["E_ads","efermi","d-position1","p-position1","s-position1"]
#vars = ["E_ads","efermi","d-position1","d-position2","p-position1","s-position1"]
#vars = ["E_ads","efermi","d-position1","d-position2","d-position3","p-position1","s-position1"]
#seaborn.pairplot(df, plot_kws={"s":10}, vars=vars, height=1.8)
#seaborn.pairplot(df, plot_kws={"s":10}, height=1.8)
#seaborn.pairplot(df, height=1.8, kind="reg")
#plt.savefig("pairplot.png")
#plt.close()

X  = df.drop("E_ads",axis=1)
y  = df["E_ads"]

cv = 4
test_size = 1.0/cv

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=0)

scaler = StandardScaler()
#scaler = MinMaxScaler()
#scaler = RobustScaler()

methods = [Ridge(), Lasso()]
names   = ["ridge", "lasso"]

print("-------- Linear Regression ---------")
lr = LinearRegression()
lr.fit(X_train, y_train)

print(pd.DataFrame({"name":X.columns, "Coef":lr.coef_}).sort_values(by="Coef"))
print("Training set score: {:.3f}".format(lr.score(X_train, y_train)))
print("Test set score: {:.3f}".format(lr.score(X_test, y_test)))
print("RMSE: {:.3f}".format( np.sqrt(mean_squared_error(y_test, lr.predict(X_test)))))

for name, method in zip(names, methods):
	print("-------- %s ---------" % name)

	pipe = make_pipeline(scaler, method)
	param_grid = {name + "__alpha":list(10**np.arange(-5, 5, 0.1))}

	grid = GridSearchCV(pipe, param_grid=param_grid, cv=cv)
	grid.fit(X_train, y_train)

	print(pd.DataFrame({"name":X.columns, "Coef":grid.best_estimator_.named_steps[name].coef_}).sort_values(by="Coef"))
	print("Best parameters: {}".format(grid.best_params_))
	print("Training set score: {:.3f}".format(grid.score(X_train, y_train)))
	print("Test set score: {:.3f}".format(grid.score(X_test, y_test)))
	print("RMSE: {:.3f}".format( np.sqrt(mean_squared_error(y_test, grid.predict(X_test)))))

	#seaborn.scatterplot(y.values, grid.predict(X))
	seaborn.regplot(y.values, grid.predict(X))
	plt.show()
#
# learning_curve
#
from sklearn.model_selection import learning_curve
best_param = list(grid.best_params_.values())[0]
pipe = Pipeline([("scl",StandardScaler()), ("lasso",Lasso(alpha=best_param))])
#pipe = Pipeline([("scl",StandardScaler()), ("lasso",Ridge(alpha=8))])
train_sizes, train_scores,test_scores = learning_curve(estimator=pipe, X=X_train, y=y_train, train_sizes=np.linspace(0.1, 1.0, 10), cv=cv)

train_mean = np.mean(train_scores, axis=1)
train_std  = np.std(train_scores,  axis=1)
test_mean  = np.mean(test_scores,  axis=1)
test_std   = np.std(test_scores,   axis=1)

plt.plot(train_sizes, train_mean, color="blue", marker="o", markersize=5, label="training accuracy")
plt.fill_between(train_sizes, train_mean+train_std, train_mean-train_std, alpha=0.15, color="blue")

plt.plot(train_sizes, test_mean, color="green", marker="s", markersize=5, label="validation accuracy")
plt.fill_between(train_sizes, test_mean+test_std, test_mean-test_std, alpha=0.15, color="green")

font = {'family': 'arial', 'weight': 'normal', 'size': 16 }

plt.xticks(fontname="arial", fontsize=12)
plt.yticks(fontname="arial", fontsize=12)
plt.xlabel("Number of training samples", fontdict=font)
plt.ylabel("Accuracy", fontdict=font)

from matplotlib import rcParams
plt.rcParams['font.family'] = "arial"
plt.rcParams['font.size']   = 10
plt.rcParams["axes.labelsize"] = 12

plt.legend(loc="lower right", fontsize=14)
plt.ylim([0.0, 1.0])
#plt.grid()
plt.show()
plt.close()
#
# Random forest or Gradient boosting forest
#
n_estimators=100
#n_components = 1

print("-------- Random Forest Regression ---------")
lr = RandomForestRegressor(n_estimators=n_estimators)
#lr = GradientBoostingRegressor(n_estimators=n_estimators)
#lr = ExtraTreesRegressor(n_estimators=n_estimators)
#lr = PLSRegression(n_components=n_components)

lr.fit(X_train, y_train)
print("Training set score: {:.3f}".format(lr.score(X_train, y_train)))
print("Test set score: {:.3f}".format(lr.score(X_test, y_test)))
print("RMSE : {}".format( np.sqrt(mean_squared_error(y_test, lr.predict(X_test)))))

plt.barh(range(lr.n_features_), lr.feature_importances_, align="center")
plt.yticks(np.arange(lr.n_features_), X.columns)
plt.show()
#
# plot correlation matrix
#
corr = df.corr()
fig, ax = plt.subplots(figsize=(12,12))
seaborn.heatmap(corr, vmax=1, vmin=-1, center=0, annot=False, annot_kws={"size":10}, cbar=True, cmap="RdBu_r", square=True, fmt=".1f", ax=ax)
plt.savefig("correlation.png")
plt.show()
plt.close()
#
# bolasso (in R)
#
import pyper
r = pyper.R(use_pandas="True")
r("source(file='do_bolasso.r')")

freq = pd.read_csv("bolasso_frequency.csv")
freq.index = freq.iloc[:,0]
freq.index.name = "Descriptors"
freq = freq.drop("Unnamed: 0", axis=1)
freq = freq.drop("Intercept", axis=0)

col = freq.columns.astype(float)
col = col.map("{:.3f}".format)
freq.columns = col

fig, ax = plt.subplots(figsize=(10,10))
seaborn.heatmap(freq, vmax=1, vmin=0, cbar=True, cmap="Blues", square=False, ax=ax)
plt.savefig("bolasso_frequency.png")
#plt.grid(color="lightgray", ls=":", linewidth=1)
plt.show()
plt.close()
