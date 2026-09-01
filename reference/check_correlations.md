# check_correlations

generates a correlation matrix among all columns of a data.frame

## Usage

``` r
check_correlations(dat, parallel = FALSE, n.cores = 4)
```

## Arguments

- dat:

  the data.frame containing the columns for which a correlation matrix
  is sought.

- parallel:

  a logical indicating if calaculation of the correlation matrix should
  be done in parallel. Defaults to FALSE.

- n.cores:

  a numeric value indicating the number of cores to utilise if parallel
  is TRUE.

## Value

a correlation matrix

## Details

The function uses cor to calculate the Pearson correlation coefficient
among continuous variables, lm to approximate the correlation
coefficient among a continuous variable and a factor variable through
the call lm(continuous~factor), and nnet to approximate the correlation
among factor variables using a multinomial model fit.

## Examples

``` r
data(case_study1)
check_correlations(case_study1[, c("depth", "complexity", "ZONE")])
#>                 depth complexity       ZONE
#> depth      1.00000000 0.33364348 0.08234682
#> complexity 0.33364348 1.00000000 0.01410832
#> ZONE       0.08234682 0.01410832 1.00000000
```
