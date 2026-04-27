# generate colors

given a bunch of items, generate a set of colors for either single node
colorings or pie-chart annotations. Colors are generated using the *hcl*
colorspace, and for `n_color >= 5`, the colors are re-ordered in an
attempt to create the largest contrasts between colors, as they result
from being picked on a circle in *hcl* space.

## Usage

``` r
generate_colors(n_color)
```

## Arguments

- n_color:

  how many colors to generate
