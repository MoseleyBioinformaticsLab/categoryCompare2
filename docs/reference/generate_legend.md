# generate a legend

it often helps to have a legend displayed for reference.

## Usage

``` r
generate_legend(
  in_assign,
  upper_names = TRUE,
  img = FALSE,
  width = 800,
  height = 400,
  pointsize = 70,
  ...
)
```

## Arguments

- in_assign:

  the assign object from `annotation_combinations`

- upper_names:

  whether to make names uppercase for easier viewing

- img:

  should a base64 encoded data uri be returned for embedding?

- width:

  how wide should the image be if saving to an image

- height:

  how high should it be

- pointsize:

  the pointsize parameter for ragg, determines textsize in the image

- ...:

  any other parameter to `pie`
