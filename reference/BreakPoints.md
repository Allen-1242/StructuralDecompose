# Generation of breakpoints

Generation of breakpoints

## Usage

``` r
BreakPoints(
  timeseries,
  frequency = 52,
  break_algorithm = "strucchange",
  break_level = 0.05
)
```

## Arguments

- timeseries:

  Given time series

- frequency:

  Timeseries frequency, defaults to 12 points

- break_algorithm:

  Breakpoint algorithm to be used. Defaults to strucchange

- break_level:

  Additional parameters for breakpoint algorithm

## Value

A list of breakpoints

## Examples

``` r
BreakPoints(timeseries = seq(100), frequency = 52, break_level = 0.05)
#>  [1]   0   5  11  17  22  28  33  38  43  48  53  59  64  69  74  79  84  89  94
#> [20] 100
BreakPoints(timeseries = StructuralDecompose::Nile_dataset[,1], frequency = 52)
#> [1]   0 100
```
