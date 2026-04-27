# annotation similarity graph

given an annotation-feature list, generate a similarity graph between
all of the annotations

## Usage

``` r
generate_annotation_similarity_graph(
  annotation_features,
  similarity_type = "combined"
)
```

## Arguments

- annotation_features:

  list where each entry is a set of features to that annotation

- similarity_type:

  which type of overlap coefficient to report

## Value

cc_graph
