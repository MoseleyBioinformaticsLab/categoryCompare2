# index a list

Provided a list, and a condition, returns the logical indices into the
named part of the list provided. Uses `subset` like non-standard
evaluation so that we can define appropriate expressions.

## Usage

``` r
multi_query_list(list_to_query, queries)
```

## Arguments

- list_to_query:

  the list to run the query on

- queries:

  the expressions that do the queries (as rlang::enquos)

## Value

logical "&" of all queries
