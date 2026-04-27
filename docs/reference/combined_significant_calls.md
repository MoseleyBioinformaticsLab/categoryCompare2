# get significant annotations calls

In the case where we have a
[`combined_enrichment`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_enrichment.md)
and we want to get all of the significant annotations from each of them,
and put them together so we can start doing real meta-analysis.

## Usage

``` r
combined_significant_calls(in_results, queries)
```

## Arguments

- in_results:

  a
  [`combined_enrichment`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_enrichment.md)
  object

- queries:

  a list of queries that can form a call object

## Value

[`combined_enrichment`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_enrichment.md)
object

## Details

Note that this function returns the original
[`combined_enrichment`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_enrichment.md)
object with a modified
[`combined_statistics`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_statistics.md)
slot where the significant annotations have been added in.
