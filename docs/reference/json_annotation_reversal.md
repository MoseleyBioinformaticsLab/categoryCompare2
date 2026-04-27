# annotation reversal

Given a JSON file of features to annotations, reverse to turn it into
annotations to features, and optionally add some meta-information about
them.

## Usage

``` r
json_annotation_reversal(
  json_file,
  out_file = "annotations.json",
  feature_type = NULL,
  annotation_type = NULL
)
```

## Arguments

- json_file:

  the json file to use

- out_file:

  the json file to write out to

- feature_type:

  the type of features

- annotation_type:

  the type of annotations

## Value

the json object, invisibly
