# build_inclusion_mat

Supporting function for functions full_subsets_gam and fit_model_set.
Not called directly.

## Usage

``` r
build_inclusion_mat(included.vars, formula.list)
```

## Arguments

- included.vars:

  A character vector of variables included in the model set

- formula.list:

  A list of model formula, as obtained through generate_model_set

## Value

A matrix of variables included in the model set

## Details

Builds var.inclusion matrix based on the included variables and set of
model formula

## Examples

``` r
included.vars <- c("depth", "complexity")
formula.list <- list(depth = ~1, "depth+complexity" = ~1)
build_inclusion_mat(included.vars, formula.list)
#>      depth complexity
#> [1,]     1          0
#> [2,]     1          1
```
