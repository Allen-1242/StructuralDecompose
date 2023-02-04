
# StructuralDecompose

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/StructuralDecompose)](https://CRAN.R-project.org/package=StructuralDecompose)
<!-- badges: end -->

Note: This package is under construction, please use the current version of R and Python here. The website will be updated in the coming days

StructuralDecompose is an algorithm suited to the decomposition of a time series into it's component terms of trend, seasonality and residuals. It is well suited to decompose a series in the presence of significant level shifts. 

The algorithm outputs the decomposed trend, seasonality, residuals as well as anomalies detected. 

Please note that the package is currently submitted to CRAN. CRAN submissions take a long time. Please use the github download to access the code. 

## Installation

You can install the development version of StructuralDecompose like so:

``` r
package(StructuralDecompose)
install_github("StructuralDecompose/StructuralDecompose")
```

## Example

An example code:

Note that we can specify the break algorithm and the smoothing algorithm as well. If the setting is set to 'auto', it will optimize which algorithm to use. 

It is best to keep the default algorithms. 

``` r
library(StructuralDecompose)

StructuralDecompose <- function(Data, frequency = 12, break_algorithm = 'strucchange', smoothening_algorithm = 'lowess', break_level = 0.05, median_level = 0.5, mean_level = 0.5, level_length = 0.5, conf_level = 0.5)

```

