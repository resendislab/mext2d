## Classifier Performance

I took the average AUC over 100 runs with kfold cross-validation (#folds = 6). Classes were balanced by random subsampling from the larger class.

Genera:
- T2D vs NGT: 0.5104166666666666
- T2D vs IGT: 0.4791666666666667
- IGT vs NGT: 0.5346200980392155

Family:
- T2D vs NGT: 0.5848484848484848
- T2D vs IGT: 0.3628787878787879
- IGT vs NGT: 0.6021241830065359

SparCC Clusters:
- T2D vs NGT: 0.575
- T2D vs IGT: 0.4469696969696971
- IGT vs NGT: 0.5334967320261437

## Important Features

For each run, I recorded the five most important features in the decision tree. I then counted the number of times each feature was one of the five most important features.

#### Genera

**T2D vs NGT**

|Genera   | Occurrence out of 100 runs|
|----------------------|---|
| Escherichia/Shigella    | 81 |
| Ruminiclostridium_5     | 57 |
| Lactobacillus           | 50 |
| Collinsella             | 49 |
| Anaerostipes            | 39 |
| Ruminococcus_2          | 37 |
| Lachnospiraceae_UCG-004 | 19 |
| Subdoligranulum         | 18 |

**IGT vs NGT**

| Genera   | Occurrence out of 100 runs|
|----------------------|---|
| Escherichia/Shigella         | 90 |
| Anaerostipes                 | 79 |
| Ruminococcus_2               | 71 |
| Streptococcus                | 45 |
| Intestinibacter              | 28 |
| Lachnospiraceae_ND3007_group | 26 |
| Fusicatenibacter             | 25 |


#### Family:

**T2D vs NGT**

| Family   | Occurrence out of 100 runs|
|-----------------------|----|
| Enterobacteriaceae    | 98 |
| Peptostreptococcaceae | 72 |
| Veillonellaceae       | 44 |
| Coriobacteriaceae     | 34 |
| Acidaminococcaceae    | 27 |
| Bifidobacteriaceae    | 26 |
| Erysipelotrichaceae   | 24 |
| Clostridiaceae_1      | 21 |

**IGT vs NGT**

| Family   | Occurrence out of 100 runs|
|-----------------------|----|
| Enterobacteriaceae    | 89 |
| Peptostreptococcaceae | 76 |
| Coriobacteriaceae     | 70 |
| Veillonellaceae       | 64 |
| Christensenellaceae   | 45 |
| Enterococcaceae       | 27 |
| Porphyromonadaceae    | 26 |
