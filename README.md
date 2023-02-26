
# StructuralDecompose

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/StructuralDecompose)](https://CRAN.R-project.org/package=StructuralDecompose)
<!-- badges: end -->

Note: Website documentation is still under construction

StructuralDecompose is an algorithm suited to the decomposition of a time series into it's component terms of trend, seasonality and residuals. It is well suited to decompose a series in the presence of significant level shifts. 

The algorithm outputs the decomposed trend, seasonality, residuals as well as anomalies detected. 


## Installation

The StructuralDecompose Package is now available on CRAN

You can install the development version of TangledFeatures like so:

| Type        | Source     | Command                                                                       |
|-------------|------------|-------------------------------------------------------------------------------|
| Release     | CRAN       | `install.packages("StructuralDecompose")`                                             |
| Development | Github     | `install_github("StructuralDecompose/StructuralDecompose")` |

Once you have downloaded the package, you can then load it using:

``` r
library(StructuralDecompose)
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

