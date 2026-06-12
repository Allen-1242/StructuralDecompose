# Minimum level length checks

Minimum level length checks

## Usage

``` r
LevelCheck(timeseries, level_length = 10, breaks)
```

## Arguments

- timeseries:

  Given time series

- level_length:

  Mean distance between two levels

- breaks:

  breakpoints returned

## Value

The series cleaned with the minimum level check

## Examples

``` r
LevelCheck(timeseries = StructuralDecompose::Nile_dataset[,1], breaks = c(1,4,5))
#> [1]   0 100

LevelCheck(timeseries = runif(n = 50, min = 1, max = 10), breaks = c(1,4,5))
#> [1]  0 50
```
