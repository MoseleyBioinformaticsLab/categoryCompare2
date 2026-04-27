# get significant annotations calls

In the case where we have a
[`statistical_results`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/statistical_results-class.md)
and we want to get all of the significant annotations from it

## Usage

``` r
get_significant_annotations_calls(in_results, queries)
```

## Arguments

- in_results:

  a
  [`statistical_results`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/statistical_results-class.md)
  object

- queries:

  a list of queries that can form a call object

## Value

vector of significant annotation_id's
