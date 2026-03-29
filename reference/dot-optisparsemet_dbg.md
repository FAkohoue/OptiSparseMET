# Internal conditional debug message helper

Emits a formatted message via
[`base::message()`](https://rdrr.io/r/base/message.html) when
`debug = TRUE`. Used internally to provide optional verbose output
during design construction and optimisation.

## Usage

``` r
.optisparsemet_dbg(debug, ...)
```

## Arguments

- debug:

  Logical. If `TRUE`, the message is emitted.

- ...:

  Arguments passed to
  [`base::sprintf()`](https://rdrr.io/r/base/sprintf.html).

## Value

Invisibly returns `NULL`.
