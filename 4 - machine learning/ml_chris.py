"""Build some classifiers."""

from sklearn.linear_model import Lasso
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
import pandas as pd
import numpy as np

data = pd.read_csv("ml_data.csv", index_col=0)
status = data.status
binstatus = (status > 1) + 0
print("Class balance: \n", binstatus.value_counts())
data = data.drop("status", axis=1)
log_data = np.log2(data + 0.5)

rf = RandomForestClassifier(500)
lassie = Lasso(0.1)

cv_rf = cross_val_score(rf, data, binstatus, scoring="roc_auc", cv=10)
cv_lasso = cross_val_score(lassie, log_data, binstatus, scoring="roc_auc",
                           cv=10)

print("\nrandom forest AUC: %f +- %f" % (cv_rf.mean(), cv_rf.std()))
print("Lasso AUC: %f +- %f" % (cv_lasso.mean(), cv_lasso.std()))

imps = pd.Series(rf.fit(data, binstatus).feature_importances_, data.columns)
coefs = pd.Series(lassie.fit(log_data, binstatus).coef_, data.columns)

print("\nTop10 RF importances:\n", imps.sort_values(ascending=False)[0:10])
print("\nTop10 Lasso coefs:\n", coefs.abs().sort_values(ascending=False)[0:10])
