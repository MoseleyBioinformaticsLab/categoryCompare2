# do hypergeometric enrichment

do hypergeometric enrichment

## Usage

``` r
hypergeometric_feature_enrichment(
  hypergeometric_features,
  direction = "over",
  p_adjust = "BH",
  min_features = 1
)
```

## Arguments

- hypergeometric_features:

  a hypergeometric_features object

- direction:

  which direction to do the enrichment (over or under)

- p_adjust:

  how to correct the p-values (default is "BH")

- min_features:

  how many features should be annotated before testing it?

## Value

enriched_result

## Details

The `min_features` argument here applies to the minumum number of
features an annotation has from the universe of features supplied,
**not** the minumum number of features from the differential list. For
more about the p-value adjustment, see
[`stats::p.adjust`](https://rdrr.io/r/stats/p.adjust.html)
