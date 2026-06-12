# Decomposition

``` r

library(StructuralDecompose)
```

## Introduction

All univariate time series data can be broken down into three major
components. A cyclic repeated motion, A general growth pattern and the
remaining unexplained portion. Popularly they are Seasonality/Cyclicity,
Trend and Residuals respectively. Each component can be treated
separately, and then summed together to get the final time series.

## Trend

Trend is the underlying motion of the time series and can have a variety
of forms. It can be increasing, decreasing or it can have highly random
movements. In StructuralDecompose, we assume that the trend can be
broken up into many different parts based upon breakpoints that we have
identified in the series.

``` r

data <- StructuralDecompose::Nile_dataset[,1]

x <- StructuralDecompose(data)
x$trend 
```

Here we have a few different versions of Trend

## Seasonality

Seasonality is a repeated motion of the time series. Seasonality can
very broadly be either Additive or multiplicative. StructuralDecompose
automatically identifies the seasonality type and treats it accordingly.

``` r

x <- StructuralDecompose(data)
x$seasonality 
```

Here we have two different types of seasonality
