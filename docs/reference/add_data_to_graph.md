# add table data to graph

given the annotation_graph and a data.frame, add all of the data in the
data.frame to the graph so it is available elsewhere. Note that for NA
integer and numerics, the value is modified to -100, and for infinite
values, it is modified to 1e100.

## Usage

``` r
add_data_to_graph(graph, data)
```

## Arguments

- graph:

  the graph to work on

- data:

  the data to add to it

## Value

graphNEL
