# Median level checks

Median level checks

## Usage

``` r
MedianCleaning(timeseries, median_level = 0.5, breaks, frequency = 52)
```

## Arguments

- timeseries:

  Given time series

- median_level:

  Median distance between two levels

- breaks:

  Breaks identified

- frequency:

  Timeseries frequency, defaults to 12 points

## Value

The series cleaned with the median check

## Examples

``` r
MedianCleaning(timeseries = StructuralDecompose::Nile_dataset[,1], breaks = c(1,4,5))
#> [1] 1 5

MedianCleaning(timeseries = runif(n = 50, min = 1, max = 10), breaks = c(1,4,5))
#> [1] 1 5
```
