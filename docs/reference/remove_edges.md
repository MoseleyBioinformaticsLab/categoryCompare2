# remove edges

given a RCy3 network connection, remove edges according to provided
values.

## Usage

``` r
remove_edges(edge_obj, cutoff, edge_attr = "weight", value_direction = "under")

# S4 method for class 'character,numeric'
remove_edges(edge_obj, cutoff, edge_attr = "weight", value_direction = "under")

# S4 method for class 'cc_graph,numeric'
remove_edges(edge_obj, cutoff, edge_attr = "weight", value_direction = "under")
```

## Arguments

- edge_obj:

  cc_graph

- cutoff:

  the cutoff to use

- edge_attr:

  which attribute to use

- value_direction:

  remove edges with value under or over

## Value

nothing

cc_graph
