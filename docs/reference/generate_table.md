# generate statistical table

given a
[`combined_enrichment`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_enrichment.md)
object, get out the data.frame either for investigation or to add data
to the
[`cc_graph`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/cc_graph.md).

## Usage

``` r
generate_table(comb_enrichment, entries = "signficant", link_type = "explicit")

# S4 method for class 'combined_enrichment'
generate_table(
  comb_enrichment,
  entries = "significant",
  link_type = "explicit"
)
```

## Arguments

- comb_enrichment:

  the
  [`combined_enrichment`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_enrichment.md)
  object

- entries:

  which entries to return, "significant" or "all"

- link_type:

  should their be an "explicit" link (see details)

## Value

data.frame

## Details

the `link_type` controls whether to create an "explicit" link that is
actually a column in the data.frame, or create an "implicit" html link
that is part of the `@name` column in the returned data.frame. Useful if
you are embedding the data.frame in an html report.
