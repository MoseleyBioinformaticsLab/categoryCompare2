[![Bioconductor Time](http://bioconductor.org/shields/years-in-bioc/categoryCompare2.svg)](http://bioconductor.org/packages/release/bioc/html/categoryCompare2.html "Bioconductor status")

[![Bioconductor Availability](http://bioconductor.org/shields/availability/release/categoryCompare2.svg)](http://bioconductor.org/packages/release/bioc/html/categoryCompare2.html#archives "Platform availability") 
[![Bioconductor Downloads](http://bioconductor.org/shields/downloads/categoryCompare2.svg)](http://bioconductor.org/packages/stats/bioc/categoryCompare2.html "Percentile downloads")
[![Bioconductor Commits](http://bioconductor.org/shields/commits/bioc/categoryCompare2.svg)](http://bioconductor.org/packages/release/bioc/html/categoryCompare2.html#svn_source "svn commits")
[![Support posts](http://bioconductor.org/shields/posts/categoryCompare2.svg)](https://support.bioconductor.org/t/categorycompare2/ "Bioconductor support posts")

[![Build Status](https://travis-ci.org/rmflight/categoryCompare2.svg?branch=master)](https://travis-ci.org/rmflight/categoryCompare2 "travis build status") [![Bioconductor Release Build](http://bioconductor.org/shields/build/release/bioc/categoryCompare.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/categoryCompare2/ "Bioconductor release build") [![Bioconductor Devel Build](http://bioconductor.org/shields/build/devel/bioc/categoryCompare2.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/categoryCompare2/ "Bioconductor devel build")

# categoryCompare2

A Bioconductor package for meta analysis of high-throughput datasets using 
enriched feature annotations instead of just the features themselves.

Note that this is a rewrite of the categoryCompare package. As of May 28, it is still only available from Github, but is expected to be part of the Fall 2015 Bioconductor (3.2) release.

See the "Description" file for additional requirements.

## Documentation

The [Vignette][vignLink] provides a description of the thinking behind
this package as well as a toy example for demonstration purposes.

## Installation

Installation of this package from Github requires the [devtools][devtoolsLink]
package.

```r
install.packages("devtools")
library(devtools)
install_github("categoryCompare2", "rmflight")
```

[vignLink]: http://rmflight.github.io/categoryCompare/index.html "categoryCompare Vignette"
[devtoolsLink]: https://github.com/hadley/devtools "devtools"

## Citation

Flight RM, Harrison BJ, Mohammad F, Bunge MB, Moon LDF, Petruska JC and Rouchka EC (2014). .CATEGORYCOMPARE, an analytical tool based on feature annotations.
_Frontiers in Genetics_. [link](http://dx.doi.org/10.3389/fgene.2014.00098)
