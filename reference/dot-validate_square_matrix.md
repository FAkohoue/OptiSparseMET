# Internal square matrix validator

Checks that an object is a base R matrix with the expected square
dimensions. Used internally to validate relationship matrices (GRM, A,
K).

## Usage

``` r
.validate_square_matrix(M, p, nm = "matrix")
```

## Arguments

- M:

  Object to validate.

- p:

  Expected dimension. Must have `nrow(M) == p` and `ncol(M) == p`.

- nm:

  Character scalar. Object name used in error messages.

## Value

Invisibly returns `TRUE` if validation passes.
