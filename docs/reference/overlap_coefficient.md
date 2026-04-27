# overlap coefficient

calculates the similarity using the "overlap" coefficient, which is

## Usage

``` r
overlap_coefficient(n1, n2)
```

## Arguments

- n1:

  group 1 of objects

- n2:

  group 2 of objects

## Value

double

## Details

length(intersect(n1, n2)) / length(union(n1, n2))
