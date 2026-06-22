# check.non.linear.correlations

generates a correlation matrix among all columns of a data.frame

## Usage

``` r
check.non.linear.correlations(dat)
```

## Arguments

- dat:

  the data.frame containing the columns for which a correlation matrix
  is sought.

## Value

an approximate correlation matrix

## Details

The function uses gam to estimate a correlation coefficient among
continuous variables (continuous~s(continuous), lm to approximate the
correlation coefficient between a continuous variable (as response) and
a factor variable (as a predictor) through the call
lm(continuous~factor), and nnet to apporoximate the correlation for
factor variables as responses using a multnomial model fit through a
call to multinom(factor~factor) or (factor~continuous).

## Note

The resulting "correlation" matrix is assymetric as the row variable is
used as the "response" and the column variable is used as the
"predictor". The use of gam may be slightly oversensitive for
continuous-continuous correlations and users may wish to increase
cor.cutoff. Inspect individual replationships manually. Values are only
approximate "correlations" and are in fact the sqrt of the R-square
values reported for each of the fitted relationships. Note that the
function assumes a gaussian distribution for continuous response
variables. Substantial deviations from this assumption will yield
spurious "correlation" estimates.

## Examples

``` r
data(case_study1)
check.non.linear.correlations(case_study1[, c("depth", "complexity", "ZONE")])
#>                 depth complexity       ZONE
#> depth      1.00000000  0.3918893 0.08234682
#> complexity 0.38159050  1.0000000 0.01410832
#> ZONE       0.07004624  0.0119875 1.00000000
```
