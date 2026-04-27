# add tooltip

before passing to Cytoscape, add a tooltip attribute to the graph

## Usage

``` r
add_tooltip(
  in_graph,
  node_data = c("name", "description"),
  description,
  separator = "\n"
)
```

## Arguments

- in_graph:

  the graph to work with

- node_data:

  which pieces of node data to use

- description:

  other descriptive text to use

- separator:

  what separator to use for the tooltip

## Value

the graph with a new nodeData member "tooltip"
