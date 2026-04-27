# do binomial testing

do binomial testing

## Usage

``` r
binomial_feature_enrichment(
  binomial_features,
  p_expected = 0.5,
  direction = "two.sided",
  p_adjust = "BH",
  conf_level = 0.95,
  min_features = 6
)
```

## Arguments

- binomial_features:

  a binomial_features object

- p_expected:

  the expected probability (default 0.5)

- direction:

  which direction to do the enrichment (two.sided, less, greater)

- p_adjust:

  how to correct the p-values (default is "BH")

- conf_level:

  the confidence level for the confidence interval (default is 0.95)

- min_features:

  a minimum number of features that are annotated to each annotation
  (default is 6)

## Value

enriched_result
