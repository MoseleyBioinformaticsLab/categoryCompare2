# combine annotations

Takes multiple
[`annotation`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/annotation.md)
objects and combines them so that there is a consistent sole set for
creating the
[`cc_graph`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/cc_graph.md)
and providing other information about each annotation entry.

## Usage

``` r
combine_annotations(annotation_list)

# S4 method for class 'list'
combine_annotations(annotation_list)
```

## Arguments

- annotation_list:

  one or more
  [`annotation`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/annotation.md)

## Value

[`annotation`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/annotation.md)
