## Analysis scripts for 16S

Analysis scripts and discussion for the MEXT2D project.

Makes use of helper functions from our [general pipeline](https://github.com/resendislab/mbtools).

### Getting set up

To make installation easier all functionality required is implemented is
delivered with a [specialized R package](https://github.com/resendislab/mbtools).
You can use the provided docker container from there or install it yourself.

You will need [R (>=3.4)](https://r-project.org). To install the package you
will need `devtools` and bioconductor. Open an R session and run the following:

```R
install.packages("devtools")
source("https://bioconductor.org/biocLite.R")
biocLite("BiocInstaller")
devtools::install_github("resendislab/mbtools")
```

Now you're ready to go :)

All individual steps are contained in its own sub-folder along with
documentation.

Input and intermediate data files can be found in the `data` folder and some
plots can be found in the `figures` folder.
