# cc_graph

A `cc_graph` class is a `graphNEL` with the added slot of `significant`,
a matrix of rows (nodes / annotations) and whether they were found to be
significant in a given enrichment (columns). This matrix is used for
classifying the annotations into different groups, and generating either
pie-charts or coloring the nodes in a visualization.

constructs a *cc_graph* given a
[`graphNEL`](https://rdrr.io/pkg/graph/man/graphNEL-class.html) and a
*significant* matrix.

## Usage

``` r
cc_graph(graph, significant)

# S4 method for class 'cc_graph'
show(object)

cc_graph(graph, significant)
```

## Arguments

- graph:

  the [`graphNEL`](https://rdrr.io/pkg/graph/man/graphNEL-class.html)

- significant:

  a matrix indicating which nodes are significant in which experiment

- object:

  the cc_graph to show

## Slots

- `significant`:

  numeric matrix of ones and zeros
