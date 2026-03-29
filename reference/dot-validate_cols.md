# Internal data frame column validator

Checks that a data frame contains all required column names. Used
internally across allocation and seed-planning functions.

## Usage

``` r
.validate_cols(dt, cols, nm = "data")
```

## Arguments

- dt:

  A data frame.

- cols:

  Character vector of required column names.

- nm:

  Character scalar. Object name used in error messages.

## Value

Invisibly returns `TRUE` if validation passes.
