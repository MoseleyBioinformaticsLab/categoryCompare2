# table from graph

Creates a table from the annotation graph, and if provided, adds the
community information to the table.

## Usage

``` r
table_from_graph(in_graph, in_assign = NULL, community_info = NULL)
```

## Arguments

- in_graph:

  the `cc_graph` object

- in_assign:

  the `node_assign` object

- community_info:

  the `community_info` object

## Value

data.frame
