# categoryCompare2 0.100.31

- changes to `gocats_2_annotation`, checking for GO terms with no description, the actual base GO terms for each portion, and adding min and max arguments for the number of features annotated to a term.

# categoryCompare2 0.100.30

- adds documentation in the CLI vignette to using GSEA on the command line.
- fixes `filter_and_group.R` to allow filtering of CLI GSEA results.
- adds documentation in the CLI vignette about the JSON structure of each of the input files, both for hypergeometric testing and GSEA.

# categoryCompare2 0.100.28

- adds gene-set enrichment analysis (GSEA) to the command line interface. Note that this is only available by running `run_enrichment.R` directly. Currently, there is no way to modify the GSEA options, such as `min_n` and `max_n`.

# categoryCompare2 0.100.25

- adds gene-set enrichment analysis (GSEA) using the `fgsea` package. See `gsea_features` and `gsea_feature_enrichment`.

# categoryCompare2 0.100.22

- Fixed a bug where if nothing was significant, then all annotations were pulled into the graph.

# categoryCompare2 0.100.21

- Added more output in the binomial testing.

# categoryCompare2 0.100.16

- Fixing bugs around binomial testing.

# categoryCompare2 0.100.9

- Better importing of gocats objects.

# categoryCompare2 0.100.8

- Made it possible to do iterative filtering and grouping of the graph in the CLI.

# categoryCompare2 0.100.7

- Fixes to the CLI.

# categoryCompare2 0.100.1

- P-value adjustment is done by default.

# categoryCompare2 0.99.161

- Binomial enrichment.

# categoryCompare2 0.99.160

- Document finding communities in the vignettes.

# categoryCompare2 0.99.159

- CLI executables are available.

# categoryCompare2 0.99.129

- Really became the try inheritor to the original categoryCompare.
