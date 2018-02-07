# Cooccurence groups in the Mexican type II diabetes cohort

We'll use [SparCC](https://bitbucket.org/yonatanf/sparcc) which is from [this paper from the Alm lab](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002687) to identify cooccurence networks amongst taxonomic groups of bacteria from this cohort.

For now, there is a single tar-compressed directory sparcc_fam_gen.tar.gz. This in includes the following files:

1. var_pass1pct_family_sparcc_corr_JA112018.txt
2. var_pass1pct_family_sparcc_corr_JA112018_pvals_2side.txt
3. var_pass1pct_family_sparcc_cov_JA112018.txt
4. var_pass1pct_genus_sparcc_corr_JA112018.txt
5. var_pass1pct_genus_sparcc_corr_JA112018_pvals_2side.txt
6. var_pass1pct_genus_sparcc_cov_JA112018.txt

The description each file is as follows:

1. A tab-delimited table of correlation values between different families.
2. A tab-delimited table of 2-sided *p*-values based on random data label shuffling (1000 permutations) for the family correlations.
3. A tab-delimited table of covariance values between different families. Not as informative or trustworthy as correlations.
4. A tab-delimited table of correlation values between different genera.
5. A tab-delimited table of 2-sided *p*-values based on random data label shuffling (1000 permutations) for the genera correlations.
6. A tab-delimited table of covariance values between different genera. Not as informative or trustworthy as correlations.
