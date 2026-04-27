# node_assign

The `node_assign` class holds the unique annotation combinations and the
assignment of the nodes to those combinations for use in visualization.

## Slots

- `groups`:

  the unique groups, as a logical matrix

- `assignments`:

  named character vector providing association with groups

- `description`:

  named character vector providing a description to group

- `colors`:

  named character vector of hex colors for groups or experiments

- `color_type`:

  whether doing group or experiment based colors

- `pie_locs`:

  if doing experiment colors, then pie graphs were generated here
