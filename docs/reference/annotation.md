# annotation class

This class holds an annotation object that defines how annotations
relate to features, as well as various pieces about each annotation

Does sensical checks when creating an `annotation` object.

## Usage

``` r
annotation(
  annotation_features,
  annotation_type = NULL,
  description = character(0),
  links = character(0),
  feature_type = NULL
)

# S4 method for class 'annotation'
show(object)

annotation(
  annotation_features,
  annotation_type = NULL,
  description = character(0),
  links = character(0),
  feature_type = NULL
)
```

## Arguments

- annotation_features:

  list of annotation to feature relationships

- annotation_type:

  a simple one word description of the annotations

- description:

  character vector providing descriptive text about the annotation

- links:

  character vector defining html links for each annotation (may be
  empty)

- feature_type:

  one word description of the feature type

- object:

  the annotation object

## Details

These objects may be created by hand, or may result from specific
functions to create them. Most notably, this package provides functions
for creating a Gene Ontology annotation.

See the `annotation`, each slot is a parameter.

## Slots

- `annotation_features`:

  list of annotation to feature relationships

- `description`:

  character vector providing descriptive text about the annotation

- `counts`:

  numeric vector of how many features are in each annotation

- `links`:

  character vector defining html links for each annotation (may be
  empty)

- `annotation_type`:

  a one word short description of the "type" of annotation

- `feature_type`:

  a one word short description of the "type" of features
