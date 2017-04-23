# categoryCompare2

A Bioconductor package for meta analysis of high-throughput datasets using 
enriched feature annotations instead of just the features themselves.

Note that this is a rewrite of the categoryCompare package. For more information about how things have changed and why, please see [this poster](https://figshare.com/articles/categoryCompare_v2_0/1427435).

## Documentation

The [Vignette](https://github.com/rmflight/categoryCompare2/blob/master/vignettes/categoryCompare_vignette_v2_visnetwork.Rmd) provides a description of the thinking behind this package as well as a toy example for demonstration purposes.

## Installation

Installation of this package from Github requires the [devtools][devtoolsLink]
package.

```r
install.packages("devtools")
library(devtools)
install_github("categoryCompare2", "rmflight")
```

[devtoolsLink]: https://github.com/hadley/devtools "devtools"

## Citation

Flight RM, Harrison BJ, Mohammad F, Bunge MB, Moon LDF, Petruska JC and Rouchka EC (2014). .CATEGORYCOMPARE, an analytical tool based on feature annotations.
_Frontiers in Genetics_. [link](http://dx.doi.org/10.3389/fgene.2014.00098)
