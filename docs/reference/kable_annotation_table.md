# print table kable

print the annotation gene table in
[`knitr::kable`](https://rdrr.io/pkg/knitr/man/kable.html) format

## Usage

``` r
kable_annotation_table(annotation_gene_table, header_level = 3, cat = TRUE)
```

## Arguments

- annotation_gene_table:

  list of tables

- header_level:

  what header level should the labels be done at?

- cat:

  whether to write it directly, or just return the table for later

## Value

character
