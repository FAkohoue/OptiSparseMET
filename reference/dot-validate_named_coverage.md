# Internal named vector coverage validator

Checks that a named vector or list contains all required names. Used
internally to validate that `env_design_specs` covers all environments.

## Usage

``` r
.validate_named_coverage(x, required, nm = "object")
```

## Arguments

- x:

  Named vector or list.

- required:

  Character vector of required names.

- nm:

  Character scalar. Object name used in error messages.

## Value

Invisibly returns `TRUE` if validation passes.
