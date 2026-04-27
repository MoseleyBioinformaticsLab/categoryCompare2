# categoryCompare2

A Bioconductor package for meta analysis of high-throughput datasets
using enriched feature annotations instead of just the features
themselves.

Note that this is a rewrite of the categoryCompare package. For more
information about how things have changed and why, please see [this
poster](https://figshare.com/articles/categoryCompare_v2_0/1427435).

## API

This version is mostly complete and works for **many** full analyses.
The user facing API is expected to be close to fixed, especially for the
core functionality. Some methods need actual `S4` or `R6` based objects
and methods, but we expect that the function calls themselves will
remain the same.

## To Do

Things that are still needed include:

- Wrapper for initial analysis when given feature data and annotations
- Example for importing the users own annotations
- Better exploration of features that are linked to specific
  annotations, including the original data associated with the features

## Documentation

The
[Vignette](https://github.com/MoseleyBioinformaticsLab/categoryCompare2/blob/master/vignettes/categoryCompare_vignette_v2_visnetwork.Rmd)
provides a description of the thinking behind this package as well as a
toy example for demonstration purposes.

## Installation

Installation of this package from Github requires the
[remotes](https://cran.r-project.org/web/packages/remotes/ "remotes")
and
[BiocManager](https://moseleybioinformaticslab.github.io/categoryCompare2/biocmanagerLink)
packages.

``` r

install.packages("BiocManager")
BiocManager::install("Biobase")
BiocManager::install("remotes")
BiocManager::install("moseleybioinformaticslab/categoryCompare2", dependencies = c("Depends", "Imports", "Suggests"))
```

### Mac Installation

There is one issue with installation on MacOS, and that is the {Cairo}
package. This is actually due to having to have `xquartz` installed. The
easiest way I can find to install it is using homebrew. In a terminal,
you can do:

    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    brew install --cask xquartz

    R
    install.packages("Cairo") # skip if already installed
    library(Cairo)

Hopefully that worked fine, and now you should be able to use all the
functionality in {categoryCompare2}.

## Citation

Flight RM, Harrison BJ, Mohammad F, Bunge MB, Moon LDF, Petruska JC and
Rouchka EC (2014). .CATEGORYCOMPARE, an analytical tool based on feature
annotations. *Frontiers in Genetics*.
[link](http://dx.doi.org/10.3389/fgene.2014.00098)
