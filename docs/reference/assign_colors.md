# assign colors

given a
[`node_assign`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/node_assign-class.md),
assign colors to either the independent groups of unique annotations, or
to each of the experiments independently.

## Usage

``` r
assign_colors(in_assign, type = "experiment")
```

## Arguments

- in_assign:

  the
  [`node_assign`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/node_assign-class.md)
  object generated from a
  [`cc_graph`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/cc_graph.md)

- type:

  either "group" or "experiment"

## Value

node_assign with colors
