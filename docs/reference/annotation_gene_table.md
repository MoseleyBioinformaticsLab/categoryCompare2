# annotation to genes

Creates a tabular output of annotations to genes providing lookup of
which genes are contributing to a particular annotation.

## Usage

``` r
annotation_gene_table(
  combined_enrichment,
  annotations = NULL,
  use_db = NULL,
  input_type = "ENTREZID",
  gene_info = c("SYMBOL", "GENENAME")
)
```

## Arguments

- combined_enrichment:

  combined enrichment object

- annotations:

  which annotations to grab features from

- use_db:

  the annotation database

- input_type:

  what type of gene id was it?

- gene_info:

  what type of info to return for each gene

## Value

data.frame
