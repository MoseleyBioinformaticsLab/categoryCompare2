# the enriched results class

given all the slots for an `enriched_result`, checks that all the data
is self-consistent, and creates the `enriched_result` object.

## Usage

``` r
hypergeometric_result(features, universe, annotation, statistics)

# S4 method for class 'enriched_result'
show(object)

hypergeometric_result(features, universe, annotation, statistics)

# S4 method for class 'binomial_result'
show(object)
```

## Arguments

- features:

  the features that were differentially expressed (see details)

- universe:

  all of the features that were measured

- annotation:

  an
  [`annotation`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/annotation.md)
  object

- statistics:

  a
  [`statistical_results`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/statistical_results-class.md)
  object

- object:

  the binomial_result object to show

## Value

hypergeometric_result

## Slots

- `annotation`:

  list giving the annotation to feature relationship

- `statistics`:

  a
  [`statistical_results`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/statistical_results-class.md)
  object

- `features`:

  the significant features

- `universe`:

  the universe of features

- `ranks`:

  the ranked features

- `positivefc`:

  the positive log-fold-changed genes, a vector of class ANY

- `negativefc`:

  the negative log-fold-changed genes
