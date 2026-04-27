# significant annotations

The `significant_annotations` class holds which annotations from which
enrichment were both **measured** and **significant**. Each of these
slots is a *logical matrix* with rows named by *annotation_id* and
columns named by the names of the
[`enriched_result`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/enriched_result.md)
that was combined.

Makes a new significant_annotation while checking that everything is
valid.

## Usage

``` r
significant_annotations(significant, measured, sig_calls = NULL)

significant_annotations(significant, measured, sig_calls = NULL)
```

## Arguments

- significant:

  logical matrix of annotations (rows) and experiments (columns)

- measured:

  logical matrix of annotations (rows) and experiments (columns)

- sig_calls:

  character vector of deparsed calls that resulted in signficant and
  measured

## Slots

- `significant`:

  logical matrix

- `measured`:

  logical matrix

- `sig_calls`:

  character representations of calls used to filter the data
