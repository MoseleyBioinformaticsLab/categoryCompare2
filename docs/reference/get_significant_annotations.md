# get significant annotations

given a
[`statistical_results`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/statistical_results-class.md)
object and some conditional expressions, return the significant
annotations

In the case where we have a
[`combined_enrichment`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_enrichment.md)
and we want to get all of the significant annotations from each of them,
and put them together so we can start doing real meta-analysis.

## Usage

``` r
get_significant_annotations(in_results, ...)

# S4 method for class 'statistical_results'
get_significant_annotations(in_results, ...)

# S4 method for class 'combined_enrichment'
get_significant_annotations(in_results, ...)
```

## Arguments

- in_results:

  a
  [`combined_enrichment`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_enrichment.md)
  object

- ...:

  conditional expressions

## Value

vector of significant annotation_id's

[`combined_enrichment`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_enrichment.md)
object

## Details

Note that this function returns the original
[`combined_enrichment`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_enrichment.md)
object with a modified
[`combined_statistics`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_statistics.md)
slot where the significant annotations have been added in.

## Examples

``` r

test_stat <- new("statistical_results",
                 annotation_id = c("a1", "a2", "a3"),
                 statistic_data = list(pvalues = c(a1 = 0.01, a2 = 0.5, a3 = 0.0001),
                   counts = c(a1 = 5, a2 = 10, a3 = 1),
                   odds = c(a1 = 20, a2 = 100, a3 = 0)))
get_significant_annotations(test_stat, pvalues < 0.05)
#> [1] "a1" "a3"
get_significant_annotations(test_stat, odds > 10)
#> [1] "a1" "a2"
get_significant_annotations(test_stat, pvalues < 0.05, counts >= 1)
#> [1] "a1" "a3"
```
