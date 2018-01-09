## Analysis scripts for 16S

Analysis scripts and discussion for the samples from Guanajuato. Also see
the [Trello board](https://trello.com/b/rHtrpyiz/microbiome) for more info.

Makes use of helper functions from our [general pipeline](https://github.com/resendislab/microbiome).

### Getting set up

To make installation easier all functionality required is implemented is
delivered with a [specialized R package](https://github.com/resendislab/microbiome).
You can use the provided docker container from there or install it yourself.

You will need [R (>=3.4)](https://r-project.org). To install the package you
will need `devtools` and bioconductor. Open an R session and run the following:

```R
install.packages("devtools")
source("https://bioconductor.org/biocLite.R")
biocLite("BiocInstaller")
devtools::install_github("cdiener/microbiome/mbtools")
```

Now you're ready to go :)

All individual steps are contained in its own sub-folder along with
documentation.

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
- genera with a `baseMean` smaller 8 (few counts in average) are removed
  since they result in a bimodal p-val histogram (no FDR adjustment methods can deal with that)
  ~~this parameter could be increased to remove low abundance genera from the analysis
  (maybe 10 would be a good cutoff)~~

