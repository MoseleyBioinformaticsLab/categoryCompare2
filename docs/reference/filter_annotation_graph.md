# filter graph by significant entries

If a graph has already been generated, it may be faster to filter a
previously generated one than generate a new one from significant data.

## Usage

``` r
filter_annotation_graph(in_graph, comb_enrich)
```

## Arguments

- in_graph:

  the
  [`cc_graph`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/cc_graph.md)
  previously generated

- comb_enrich:

  the
  [`combined_enrichment`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_enrichment.md)
  that you want to use to filter with

## Value

cc_graph
