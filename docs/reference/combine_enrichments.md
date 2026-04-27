# combine enrichments

This is one of the primary workhorse functions behind categoryCompare2.
The primary function of `categoryCompare` is to enable *comparisons* of
different enrichment analyses. To facilitate that, we must first
**combine** one (really, we can do this with a single) or more
[`enriched_result`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/enriched_result.md).

## Usage

``` r
combine_enrichments(...)

# S4 method for class 'enriched_result'
combine_enrichments(...)

# S4 method for class 'list'
combine_enrichments(...)
```

## Arguments

- ...:

  list of enriched_result

## Value

[`combined_enrichment`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_enrichment.md)
