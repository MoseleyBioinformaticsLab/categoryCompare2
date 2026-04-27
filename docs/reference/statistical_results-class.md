# statistical results class

This class holds the part of an enrichment that is the statistical
results. It has two pieces, a `list` of `statistics` that is a named
list with the actual numerical results of applying the statistics. The
other piece is the `annotation_id` vector defining which entry in each
vector of the `statistics` is.

## Slots

- `statistic_data`:

  list of numerical statistics

- `annotation_id`:

  vector of ids

- `method`:

  how the statistics were calculated
