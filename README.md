## Analysis scripts fro 16S

Analysis scripts and discussion for the samples from Guanajuato. Also see
the [Trello board](https://trello.com/b/rHtrpyiz/microbiome) for more info.

Makes use of helper functions from our [general pipeline](https://github.com/resendislab/microbiome).

Still in a very rudimentary stage. After discussion this should be converted
to R notebooks.

### Trimming and filtering

This was done using the `filterAndTrim` method from dada2. > 75% of reads pass
the filters.

**Parameters**:

- trimLeft: 10bp (standard for Illumina)
- truncLen: 240 and 200 (forward and reverse)
- maxEE: 2 (max expected errors in read using the Illumina error model)

### Sequence variants

**source file(s)**: `build_taxonomy.R`

Basically goes all the way from filtered FastQ to sequence variants + taxonomy.

**Comments:**

- error models are fitted to each run individually as recommended by Ben
- `anntotate_files` will probably become obsolete since the only thing you really
  need is the ID ("folio")
- there are many bimeras, however most of the reads are kept after bimera-removal
  (> 70%)
- taxonomy assignment uses RDP with the SILVA DB v128
- still missing the dada2 species assignment by exact match
- final product is a [phyloseq](https://joey711.github.io/phyloseq/) object

### Differential analysis (association tests)

**source file(s)**: `testing.R`

Tests associations between genera and clinical variables

**Comments**:

- some samples are duplicated and have to be omitted since those are samples taken
  at different time-points (unfortunately the info which time point each samples
  corresponds has been lost)
- "poscount" size factor estimation was recommended by Mike Love and is based on
  one of their preprints where this gave the best results for microbiome data
- currently uses only gender as confounder, others are possible but not be missing
  for the tested target variable
- p-values are adjusted using [independent-hypothesis weighting](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4930141/)
  using the mean normalized count across all samples. This is supposed to remove
  the dependence of testing power on the actual abundance. There is no discernible
  enrichment of low p-values for low abundance genera (rather the reverse is true)
- genera with a `baseMean` smaller 1 (not even a single count in average) are removed
  since they result in a bimodal p-val histogram (no FDR adjustment methods can deal with that)
  this parameter could be increased to remove low abundance genera from the analysis
  (maybe 10 would be a good cutoff)

