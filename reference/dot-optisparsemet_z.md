# Internal z-score normalisation helper

Computes the z-score of a numeric vector. If the standard deviation is
zero or non-finite, returns a zero vector of the same length.

## Usage

``` r
.optisparsemet_z(x)
```

## Arguments

- x:

  Numeric vector.

## Value

Numeric vector of z-scores, or zeros if `sd(x) == 0`.
