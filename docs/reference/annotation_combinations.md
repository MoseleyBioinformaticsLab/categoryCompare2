# unique annotation combinations

determine the unique combinations of annotations that exist in the
significant matrix of the
[`cc_graph`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/cc_graph.md)
and assign each node in the graph to a group.

determine the unique combinations of annotations that exist in the
significant matrix of the
[`combined_statistics`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_statistics.md)
and assign each annotation to a group.

## Usage

``` r
annotation_combinations(object)

# S4 method for class 'cc_graph'
annotation_combinations(object)

# S4 method for class 'significant_annotations'
annotation_combinations(object)
```

## Arguments

- object:

  the
  [`combined_statistics`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_statistics.md)
  to work on

## Value

node_assignment

node_assignment
