# combined statistics

holds the results of extracting a bunch of statistics from a
[`combined_enrichment`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/combined_enrichment.md)
into one entity. This is useful because we want to enable multiple data
representations and simple filtering on the actual `data.frame` of
statistics, and this provides flexibility to enable that.

constructor function for the combined_statistics object, makes sure that
empty things get initialized correctly

## Usage

``` r
combined_statistics(
  statistic_data,
  which_enrichment,
  which_statistic,
  annotation_id,
  significant = NULL,
  measured = NULL,
  use_names = NULL
)

combined_statistics(
  statistic_data,
  which_enrichment,
  which_statistic,
  annotation_id,
  significant = NULL,
  measured = NULL,
  use_names = NULL
)
```

## Arguments

- statistic_data:

  the data.frame of statistics

- which_enrichment:

  which enrichment gave the results

- which_statistic:

  which statistics were calculated in each case

- annotation_id:

  the annotations for which we are returning statistics

- significant:

  the significant annotations

- measured:

  the measured annotations

- use_names:

  the order of naming

## Value

combined_statistics

## Slots

- `statistic_data`:

  a `data.frame` of all of the statistics from all of the enrichments

- `significant`:

  a
  [`significant_annotations`](https://moseleybioinformaticslab.github.io/categoryCompare2/reference/significant_annotations.md)
  object, that may be empty

- `which_enrichment`:

  a `vector` giving which enrichment each column of the statistics came
  from

- `which_statistic`:

  a `vector` providing which statistic each column contains
