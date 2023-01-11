# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#StructuralDecompose

#Novel Method to decompose a levelShifted time series
#Temp change

#library import
library(data.table)
library(dplyr)

library(changepoint)
library(segmented)
library(strucchange)

#library(smooth)

#Testing with a toy dataset
#data1 <- read.csv("C:\\Users\\sunny\\OneDrive\\Documents\\GitHub\\StructuralDecompose\\Nile_DataSet.csv", header=TRUE, stringsAsFactors=FALSE)
#Data <- data1$x

#Replace this with Twitters anomaly detection method, advanded methods

#' Generation of breakpoints
#'
#' @param timeseries Given time series
#' @param frequency Tiemseries frequency, defaults to 12 points
#' @param break_algorithm Breakpoint algorithm to be used. Defaults to strcchange
#' @param break_level Additional parameters for breakpoint algorithm
#'
#' @return A list of breakpoints
#' @export
#'
#' @examples
BreakPoints <- function(timeseries, frequency = 52, break_algorithm = 'strucchange', break_level = 0.05)
{
  #Strucchange algorithm
  if(break_algorithm == 'strucchange')
  {
    mvalue <- NA
    bp <- NA

    tryCatch(
      {
        mvalue = strucchange::breakpoints(timeseries ~ 1, h = break_level)
      }, error = function(e){bp <<- NA}
    )


    if(!(is.na(mvalue$breakpoints)))
    {
      bp = mvalue$breakpoints
    }

    if(is.na(bp)){print('Change break value , min segment size must be larger than the number of regressors')}

  }

  #Changpoint
  if(break_algorithm == 'changepoint')
  {
    changepoints <- changepoint::cpt.mean(timeseries, method="BinSeg")

    bp =  cpts(changepoints)

    if(is.na(bp)){print('Change break value , min segment size must be larger than the number of regressors')}

  }

  #Segmented
  if(break_algorithm == 'segmented')
  {
    changepoints <- changepoint::cpt.mean(timeseries, method="BinSeg")

    bp =  segmented(timeseries)$psi

    if(is.na(bp)){print('Change break value , min segment size must be larger than the number of regressors')}

  }

  #Writing the breakpoints
  breaks <- vector()

  if(is.na(bp))
  {
    breaks <- c(0, length(timeseries))
  }else
  {
    breaks <- c(0, bp, length(timeseries))
  }


  return(breaks)

}

#' Median level checks
#'
#' @param timeseries Given time series
#' @param median_level Median distance between two levels
#' @param Breakpoints Breakpoints identified
#' @param breaks Breaks identified
#' @param frequency Tiemseries frequency, defaults to 12 points
#'
#' @return The series cleaned with the median check
#' @export
#'
#' @examples
MedianCleaning <- function(timeseries, median_level = 0.5, Breakpoints = c(), breaks, frequency = 52)
{
  maj_filt_rast<- focal(input_rast, fun=modal, w=matrix(1,nrow=3,ncol=3))
}


#' Automatic Anomaly detection
#'
#' @param timeseries Given time series
#' @param frequency Tiemseries frequency, defaults to 12 points
#' @param conf_level Confidence level for Anomaly detection
#' @param breaks breakpoints identified
#' @param window_len Window length for anomaly detection
#'
#' @return the list of anomalies in the time series, along with the time series plot
#' @export
#'
#' @examples
AnomalyDetection <- function(timeseries, frequency = 52, conf_level = 0.05, breaks, window_len = 14)
{    if(is.null(num_obs_per_period)){
  stop("must supply period length for time series decomposition")
}

  num_obs <- nrow(data)

  # Check to make sure we have at least two periods worth of data for anomaly context
  if(num_obs < num_obs_per_period * 2){
    stop("Anom detection needs at least 2 periods worth of data")
  }

  # Check if our timestamps are posix
  posix_timestamp <- if (class(data[[1L]])[1L] == "POSIXlt") TRUE else FALSE

  # Handle NAs
  if (length(rle(is.na(c(NA,data[[2L]],NA)))$values)>3){
    stop("Data contains non-leading NAs. We suggest replacing NAs with interpolated values (see na.approx in Zoo package).")
  } else {
    data <- na.omit(data)
  }

  # -- Step 1: Decompose data. This returns a univarite remainder which will be used for anomaly detection. Optionally, we might NOT decompose.
  data_decomp <- stl(ts(data[[2L]], frequency = num_obs_per_period),
                     s.window = "periodic", robust = TRUE)

  # Remove the seasonal component, and the median of the data to create the univariate remainder
  data <- data.frame(timestamp = data[[1L]], count = (data[[2L]]-data_decomp$time.series[,"seasonal"]-median(data[[2L]])))

  # Store the smoothed seasonal component, plus the trend component for use in determining the "expected values" option
  data_decomp <- data.frame(timestamp=data[[1L]], count=(as.numeric(trunc(data_decomp$time.series[,"trend"]+data_decomp$time.series[,"seasonal"]))))

  if(posix_timestamp){
    data_decomp <- format_timestamp(data_decomp)
  }
  # Maximum number of outliers that S-H-ESD can detect (e.g. 49% of data)
  max_outliers <- trunc(num_obs*k)

  if(max_outliers == 0){
    stop(paste0("With longterm=TRUE, AnomalyDetection splits the data into 2 week periods by default. You have ", num_obs, " observations in a period, which is too few. Set a higher piecewise_median_period_weeks."))
  }

  func_ma <- match.fun(median)
  func_sigma <- match.fun(mad)

  ## Define values and vectors.
  n <- length(data[[2L]])
  if (posix_timestamp){
    R_idx <- as.POSIXlt(data[[1L]][1L:max_outliers], tz = "UTC")
  } else {
    R_idx <- 1L:max_outliers
  }

  num_anoms <- 0L

  # Compute test statistic until r=max_outliers values have been
  # removed from the sample.
  for (i in 1L:max_outliers){
    if(verbose) message(paste(i,"/", max_outliers,"completed"))

    if(one_tail){
      if(upper_tail){
        ares <- data[[2L]] - func_ma(data[[2L]])
      } else {
        ares <- func_ma(data[[2L]]) - data[[2L]]
      }
    } else {
      ares = abs(data[[2L]] - func_ma(data[[2L]]))
    }

    # protect against constant time series
    data_sigma <- func_sigma(data[[2L]])
    if(data_sigma == 0)
      break

    ares <- ares/data_sigma
    R <- max(ares)

    temp_max_idx <- which(ares == R)[1L]

    R_idx[i] <- data[[1L]][temp_max_idx]

    data <- data[-which(data[[1L]] == R_idx[i]), ]

    ## Compute critical value.
    if(one_tail){
      p <- 1 - alpha/(n-i+1)
    } else {
      p <- 1 - alpha/(2*(n-i+1))
    }

    t <- qt(p,(n-i-1L))
    lam <- t*(n-i) / sqrt((n-i-1+t**2)*(n-i+1))

    if(R > lam)
      num_anoms <- i
  }

  if(num_anoms > 0) {
    R_idx <- R_idx[1L:num_anoms]
  } else {
    R_idx = NULL
  }}

#' Smoothening of the time series
#'
#' @param timeseries Given time series
#' @param frequency Tiemseries frequency, defaults to 12 points
#' @param smoothening_algorithm Smoothening algorithm required
#' @param breaks Breakpoints identified
#'
#' @return The smoothened time series
#' @export
#'
#' @examples
Smoothing <- function(timeseries, frequency = 52, smoothening_algorithm = 'lowess', breaks = Break_points)
{
  #Smoothening the series
  k <- vector()
  k <- breaks
  k <- sort(k)

  trend_line <- vector()
  if(length(k) == 0)
  {
    smooth_series = lowess(timeseries)
    trend_line = c(trend_line, smooth_series$y)
  }else
  {
    trend_line <- vector()
    pointer <- vector()

    for(i in seq(from = 0, to = c(length(k) - 2), by = 1))
    {
      if(i == 0)
      {
        pointer <- 0
      }else
      {
        pointer <- 1
      }

      v <- timeseries[(c(k[i+1]+pointer) : k[i+2])]

      if(smoothening_algorithm == 'lowess')
      {
        smooth_series = lowess(v)$y
      }

      if(smoothening_algorithm == 'SMA')
      {
        #smooth_series = smooth::sma(v)
      }

      trend_line = c(trend_line, smooth_series)
    }
  }

  return(trend_line)

}

#' Main decomposition algorithm
#'
#' @param Data Time series required
#' @param frequency Frequency of the tine series
#' @param break_algorithm breakpoints algorithm used. Defaults to strucchange
#' @param smoothening_algorithm Smoothing algorithm used. Defaults to lowess
#' @param break_level Break level for the breakpoints algorithm
#' @param median_level Average median distance between two level
#' @param mean_level Average mean distance between a group of points near breakpoints
#' @param level_length Minimum number of points required to determine a level
#' @param conf_level Confidence interval used for Anomaly detection
#'
#' @return The decomposed time series along with a host of other metrics
#' @export
#'
#' @examples
StructuralDecompose <- function(Data, frequency = 12, break_algorithm = 'strucchange', smoothening_algorithm = 'lowess', break_level = 0.05, median_level = 0.5, mean_level = 0.5, level_length = 0.5, conf_level = 0.5)
{
  #Strucchange algorithm

  #Initial Sanity checks
  if(!is.numeric(frequency)  || !is.numeric(break_level) || !is.numeric(mean_check) || !is.numeric(median_level) || !is.numeric(level_length) || !is.numeric(conf_level) || !is.numeric(window_len))
  {
    stop(print('Value needs to be numeric'))
  }

  if(!is.logical(plot))
  {
    stop(print('Value needs to be boolean'))
  }

  #If the time series contains NA
  if(any(is.na(y)))
  {
    stop(print('Interpolation of time series needed , recommend the zoo package'))
  }


  #Calling the main break-point algorithm
  Break_points <- BreakPoints(timeseries = Data, frequency = frequency, break_algorithm = break_algorithm, break_level = break_level)

  #Median Cleaning
  Break_points <- MedianCleaning(timeseries = Data, breaks = Break_points, frequency = frequency, median_level = 0.5)

  #Mean Cleaning
  Break_points <- MeanCleaning(timeseries = Data, breaks = Break_points, frequency = frequency, mean_level = 0.5)

  #Level check
  Break_points <- LevelCheck(timeseries = Data, breaks = Break_points, level_length = 10)


  #Anomaly Detection
  Cleanseries <- AnomalyDetection(timeseries = Data, frequency = frequency, breaks = Break_points, conf_level = 0.5)

  #Smoothing of the time series
  Decomposedtrend <- Smoothing(timeseries = Cleanseries, breaks = Break_points)

  #Transforming into a series, final decomp
  ts.QTY1 = ts(data = as.vector(t(Decomposedtrend)), frequency = frequency)

  decomposed <- NA
  tryCatch(
    {
      decomposed <- stl(ts.QTY1, s.window = 'periodic')
    }, error = function(e){ts_QTY1 <<- NULL}
  )

  if(length(decomposed) != 1)
  {
    seasonal <- decomposed$time.series[,1]
    trend <- decomposed$time.series[,2]
    remainder <- decomposed$time.series[,3]

    #Removing seasonality
    ts_QTY1 <- ts.QTY1 - seasonal
  }else
  {
    seasonal <- NULL
    trend <- NULL
    remainder <- NULL
  }

  newList <- list('anomalies' = anomalies, 'trend_line' = trend_line, 'ds_series' = ts_QTY1,
                  'breakpoints' = k)

  return(newList)
}
