## Differential analysis (association tests)

**source file(s)**: `testing.R`

Tests associations between genera and clinical variables.

**Parameters**:

All of those are used for filtering out specific taxa.

- maximum percentage of samples where the taxa was not observed: 10%
- minimum normalized average abundance across all samples: 10 reads
- acceptable diabetes status codes: 1-5

**Comments**:

- some samples are duplicated and have to be omitted since those are samples taken
  at different time-points (unfortunately the info which time point each samples
  corresponds has been lost)
- "poscount" size factor estimation was recommended by Mike Love and is based on
  one of their preprints where this gave the best results for microbiome data
- currently uses only gender as confounder, others are possible but not be missing
  for the tested target variable
- p-values are adjusted using independent hypothesis weighting (https://www.nature.com/articles/nmeth.3885).
  There is no discernible enrichment of low p-values for low abundance genera
  (rather the reverse is true)
- genera with a `baseMean` smaller 10 (few counts in average) are removed
  since they result in a bimodal p-val histogram (no FDR adjustment methods can deal with that)
  ~~this parameter could be increased to remove low abundance genera from the analysis
  (maybe 10 would be a good cutoff)~~

