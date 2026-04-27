# combined enrichments

The `combined_enrichment` class holds the results of combining several
[`enriched_result`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/enriched_result.md)s
together, which includes the original
[`enriched_result`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/enriched_result.md)s,
as well as the
[`cc_graph`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/cc_graph.md)
and combined
[`annotation`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/annotation.md)
objects.

## Slots

- `enriched`:

  list of enriched objects

- `annotation`:

  [`annotation`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/annotation.md)
  where the annotation_features have been combined across the
  [`enriched_result`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/enriched_result.md)

- `statistics`:

  [`combined_statistics`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_statistics.md)
  of both
