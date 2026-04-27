# extract statistics

extract all statistics from a
[`combined_enrichment`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_enrichment.md)
object and create a
[`combined_statistics`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_statistics.md)
where each statistic from the underlying
[`statistical_results`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/statistical_results-class.md)
object in each of the enrichments is named according to which enrichment
it was in and what statistic it was.

## Usage

``` r
# S4 method for class 'combined_enrichment'
extract_statistics(in_results)
```

## Arguments

- in_results:

  the
  [`combined_enrichment`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_enrichment.md)
  object

## Value

combined_statistics
