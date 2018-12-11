import pandas as pd
import numpy as np
from scipy import interp
import matplotlib.pyplot as plt

from sklearn.metrics import roc_curve, auc
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold

data = pd.read_csv("ml_data.csv", index_col=0)
status = data.status

#bin patients into two categories: no diabetes or have some degree of diabetes
binstatus = (status > 1) + 0
data = data.drop("status", axis=1)
#take the log of abundances
log_data = np.log2(data + 0.5)
#initiate RF classifier
rf = RandomForestClassifier(500)

cv_rf = cross_val_score(rf, data, binstatus, scoring="roc_auc", cv=10)

imps = pd.Series(rf.fit(data, binstatus).feature_importances_, data.columns)

cv = StratifiedKFold(n_splits=10)

tprs = [] #true positive rates
aucs = [] 
mean_fpr = np.linspace(0, 1, 100)

with open("rf_classification_results.txt", "w+") as f:
    f.write("Class balance: \n")
    f.write(str(binstatus.value_counts()))
    f.write(str("\nrandom forest AUC: %f +- %f" % (cv_rf.mean(), cv_rf.std())))
    f.write("\nTop10 RF importances:\n")
    f.write(str(imps.sort_values(ascending=False)[0:10]))
i = 1
for train, test in cv.split(data, binstatus):

    probas_ = rf.fit(data.iloc[train], binstatus.iloc[train]).predict_proba(data.iloc[test])
    # Compute ROC curve and area under the curve
    fpr, tpr, thresholds = roc_curve(binstatus.iloc[test], probas_[:, 1])
    tprs.append(interp(mean_fpr, fpr, tpr))
    tprs[-1][0] = 0.0
    roc_auc = auc(fpr, tpr)
    aucs.append(roc_auc)
    plt.plot(fpr, tpr, lw=1, alpha=0.3,
             )

    i += 1
plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
         label='AUC = 0.5', alpha=.8)

mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
plt.plot(mean_fpr, mean_tpr, color='b',
         label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
         lw=2, alpha=.8)

std_tpr = np.std(tprs, axis=0)
tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
tprs_lower = np.maximum(mean_tpr - std_tpr, 0)

#plotting
plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                 label=r'$\pm$ 1 std. dev.')

plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Diabetic Status RF Classifier ROC')
plt.legend(loc="lower right")
plt.savefig("rf_roc.png", format="png")
plt.savefig("rf_roc.svg", format="svg")
