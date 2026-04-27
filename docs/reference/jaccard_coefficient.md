# jaccard coefficient

calculates similarity of two groups of objects using "jaccard"
coefficient, defined as:

## Usage

``` r
jaccard_coefficient(n1, n2)
```

## Arguments

- n1:

  group 1

- n2:

  group 2

## Value

double

## Details

length(intersect(n1, n2)) / min(c(length(n1), length(n2)))
