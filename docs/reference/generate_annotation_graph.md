# generate annotation graph

given a
[`combined_enrichment`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_enrichment.md),
generate the annotation similarity graph

## Usage

``` r
generate_annotation_graph(
  comb_enrichment,
  annotation_similarity = "combined",
  low_cut = 5,
  hi_cut = 500
)

# S4 method for class 'combined_enrichment'
generate_annotation_graph(
  comb_enrichment,
  annotation_similarity = "combined",
  low_cut = 5,
  hi_cut = 500
)
```

## Arguments

- comb_enrichment:

  the combined_enrichment object

- annotation_similarity:

  which similarity measure to use

- low_cut:

  keep only those annotations in the graph with at least this many
  annotated features

- hi_cut:

  keep only those annotations with less than this many annotated
  features

## Value

[`cc_graph`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/cc_graph.md)
