# gocats to annnotations

Transforms a gocats ancestors JSON list to a GO annotation object.

## Usage

``` r
gocats_to_annotation(
  ancestors_file = "ancestors.json",
  namespace_file = "namespace.json",
  annotation_type = "gocatsGO",
  feature_type = "Uniprot",
  feature_translation = NULL,
  feature_min = 5,
  feature_max = 5000,
  add_description = "yes"
)
```

## Arguments

- ancestors_file:

  the ancestors.json file from gocats (required)

- namespace_file:

  the namespace.json file from gocats (optional)

- annotation_type:

  what annotations are we making? (gocatsGO by default)

- feature_type:

  what type of features are we using (assume Uniprot)

- feature_translation:

  a data.frame used to convert the feature IDs

- feature_min:

  minimum number of features annotated (default: 5)

- feature_max:

  maximum number of features annotated (default: 5000)

- add_description:

  should GO.db be queried for annotation descriptions ("yes")

## Value

annotation object
