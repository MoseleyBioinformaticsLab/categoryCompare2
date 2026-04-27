# do GSEA

Performs gene-set enrichment analysis using the \`fgsea\` package.

## Usage

``` r
gsea_feature_enrichment(
  gsea_features,
  min_features = 15,
  max_features = 500,
  return_type = "cc2",
  ...
)
```

## Arguments

- gsea_features:

  a GSEA features object

- min_features:

  the minimum number of features for an annotation (default = 15)

- max_features:

  the maximum number of features for an annotation (default = 500)

- return_type:

  what type of object should be returned? ("cc2" or "fgsea")

- ...:

  other \`fgsea\` options

## Value

enriched_result

## Details

The runtime is dependent on the maximum size of the provided annotation,
so the authors of \`fgsea\` recommend a maximum size of 500. In
addition, to calculate statistics, a minimum size of annotated features
are required. Going below 15 may not be advised. If you want to use
other \`fgsea\` functions, it is recommended to set \`return_type =
"fgsea"\`. Otherwise, you should keep the default of "cc2".

## See also

fgsea::fgsea
