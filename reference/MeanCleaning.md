# Mean level checks

Mean level checks

## Usage

``` r
MeanCleaning(timeseries, mean_level = 0.5, breaks, frequency = 52)
```

## Arguments

- timeseries:

  Given time series

- mean_level:

  Mean distance between two levels

- breaks:

  breakpoints returned

- frequency:

  Timeseries frequency, defaults to 12 points

## Value

The series cleaned with the mean check

## Examples

``` r
MeanCleaning(timeseries = StructuralDecompose::Nile_dataset[,1], breaks = c(1,4,5), frequency = 1)
#> Error in MeanCleaning(timeseries = StructuralDecompose::Nile_dataset[,     1], breaks = c(1, 4, 5), frequency = 1): unused argument (frequency = 1)

MeanCleaning(timeseries = runif(n = 50, min = 1, max = 10), breaks = c(1,4,5), frequency = 12)
#> Error in MeanCleaning(timeseries = runif(n = 50, min = 1, max = 10), breaks = c(1,     4, 5), frequency = 12): unused argument (frequency = 12)
```
