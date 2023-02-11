#StructuralDecompose
#Novel Method to decompose a levelShifted time series


#' Nile River Dataset
#' @name Nile_dataset
#' @keywords datasets
NULL

#' Generation of breakpoints
#'
#' @param timeseries Given time series
#' @param frequency Timeseries frequency, defaults to 12 points
#' @param break_algorithm Breakpoint algorithm to be used. Defaults to strucchange
#' @param break_level Additional parameters for breakpoint algorithm
#'
#' @return A list of breakpoints
#' @examples
#' BreakPoints(timeseries = seq(100), frequency = 52, break_level = 0.05)
#' BreakPoints(timeseries = StructuralDecompose::Nile_dataset[,1], frequency = 52)
#' @importFrom changepoint cpts
#' @importFrom segmented segmented


#' @export
#'
BreakPoints <- function(timeseries, frequency = 52, break_algorithm = 'strucchange', break_level = 0.05)
{
  #Strucchange algorithm
  if(break_algorithm == 'strucchange')
  {
    mvalue <- NA
    bp <- NA

    tryCatch(
      {
        mvalue <- strucchange::breakpoints(timeseries ~ 1, h = break_level)
      }, error = function(e){bp <<- NA}
    )
    bp <- mvalue$breakpoints
  }

  #Changepoint
  if(break_algorithm == 'changepoint')
  {
    changepoints <- changepoint::cpt.mean(timeseries, method="BinSeg")

    bp <-  cpts(changepoints)

    if(is.na(bp)){warning('Change break value , min segment size must be larger than the number of regressors')}

  }

  #Segmented
  if(break_algorithm == 'segmented')
  {
    changepoints <- segmented::segmented(timeseries)$psi

    bp <-  changepoints

    if(is.na(bp)){warning('Change break value , min segment size must be larger than the number of regressors')}

  }

  #Writing the breakpoints
  breaks <- vector()

  if(any(is.na(bp)))
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
#' @param breaks Breaks identified
#' @param frequency Timeseries frequency, defaults to 12 points
#' @examples
#' MedianCleaning(timeseries = StructuralDecompose::Nile_dataset[,1], breaks = c(1,4,5))
#'
#' MedianCleaning(timeseries = runif(n = 50, min = 1, max = 10), breaks = c(1,4,5))
#' @return The series cleaned with the median check
#' @export
#'
MedianCleaning <- function(timeseries, median_level = 0.5, breaks, frequency = 52)
{

  #Writing the breakpoints
  t <- vector()

  if(any(is.na(breaks)))
  {
    t <- c(0, length(timeseries))
  }else
  {
    t <- breaks
  }

  p1 <- vector()
  p2 <- vector()
  #Seasonal check for level changes
  if(frequency != 1)
  {
    difference_mat <- outer(t,t, "-")
    difference_table <- which(difference_mat == frequency , arr.ind = TRUE)
    difference_table <- data.frame(difference_table)

    if(dim(difference_table)[1] != 0)
    {
      for(l in seq (from = 1 , to = c(dim(difference_table)[1]), by = 1))
      {
        p1 <- c(t[difference_table['row'][l,]], t[difference_table['col'][l,]])
        p2 <- c(p2, p1)
      }
    }

  }else
  {
    p2 <- t
  }

  t <- t[!(t %in% p2)]

  #Median cleaning of breakpoints
  med_flag <- FALSE
  n <- length(timeseries)
  k <- vector()
  i <- 1
  j <- 3

  while(i < n)
  {
    while(j < n+1)
    {
      L <- t[i]
      R <- t[j]
      Lr <- t[i+1]

      if((is.na(R)))
      {
        med_flag <- TRUE
        k <- c(0, length(timeseries))

        return(k)
      }

      mid1 <- median(timeseries[c(L+1) : c(Lr)])
      mid2 <- median(timeseries[c(Lr + 1)] : (R))

      if(mid2 == 0)
      {
        mid2 <- 0.0000001
      }

      if((abs(mid1/mid2) > median_level) || abs(mid2/mid1 > median_level))
      {
        k <- c(k , Lr)

        if(is.na(k))
        {
          k <- c(0, length(timeseries))
        }else
        {
          k <- c(0, k, length(timeseries))
        }

        return(k)

      }else
      {
        t <- t[t!=Lr]
      }
    }

    i <- i + 1
    j <- j + 1
    L <- t[i]

    if(med_flag == TRUE)
    {
      if(is.na(k))
      {
        k <- c(0, length(timeseries))
      }else
      {
        k <- c(0, k, length(timeseries))
      }

      return(k)
    }

  }

  #Writing the breakpoints
  if(any(is.na(k)))
  {
    k <- c(0, length(timeseries))
  }else
  {
    k <- c(0, k, length(timeseries))
  }

  return(k)


}


#' Mean level checks
#'
#' @param timeseries Given time series
#' @param mean_level Mean distance between two levels
#' @param breaks breakpoints returned
#' @param frequency Timeseries frequency, defaults to 12 points
#' @examples
#' MeanCleaning(timeseries = StructuralDecompose::Nile_dataset[,1], breaks = c(1,4,5), frequency = 1)
#'
#' MeanCleaning(timeseries = runif(n = 50, min = 1, max = 10), breaks = c(1,4,5), frequency = 12)
#' @return The series cleaned with the mean check
#' @export
#'
MeanCleaning <- function(timeseries, mean_level = 0.5, breaks, frequency = 52)
{
  breaks_new <- vector()
  for(x in seq(0, c(length(breaks)), 1))
  {
    if(x == c(length(breaks) - 1))
    {
      breaks_new <- unique(c(0, breaks_new, length(timeseries)))
      return(breaks_new)
    }else
    {
      before <- timeseries[c(breaks[[x + 1]] - 4) : c(c(breaks[[x + 1]]) - 1)]
      rownames(before) <- NULL

      after <- timeseries[c(breaks[[x + 1]] + 1): c(c(breaks[[x + 1]]) + 4)]
      rownames(after) <- NULL

      #Getting the mean values
      M1 <- mean(before)
      M2 <- mean(after)

      #Checking the mean ratio
      if((abs(M1/M2)) > mean_level || ((abs(M2/M1)) > mean_level))
      {
        breaks_new <- c(breaks, breaks[x + 1])
      }

    }
  }

  return(breaks_new)
}

#' Minimum level length checks
#'
#' @param timeseries Given time series
#' @param level_length Mean distance between two levels
#' @param breaks breakpoints returned
#' @examples
#' LevelCheck(timeseries = StructuralDecompose::Nile_dataset[,1], breaks = c(1,4,5))
#'
#' LevelCheck(timeseries = runif(n = 50, min = 1, max = 10), breaks = c(1,4,5))
#' @return The series cleaned with the minimum level check
#' @importFrom utils tail
#' @importFrom stats median

#' @export
#'
LevelCheck <- function(timeseries, level_length = 10, breaks)
{
  breaks_new <- vector()

  if(length(breaks) != 2)
  {

    for(i in seq(from = 1, to = c(length(breaks) - 1), by = 1))
    {
      if(((breaks[i+1] - breaks[i] >= level_length)))
      {
        breaks_new <- c(breaks[i+1], breaks_new)
      }
    }

    breaks_new <- c(0, breaks_new, length(timeseries))
    breaks_new <- unique(sort(breaks_new))


    if(length(breaks) == 2)
    {
      return(breaks_new)

    }else if(c(length(timeseries) - tail(breaks, 2)[1] <= level_length))
    {
      breaks_new <- breaks_new[-match(tail(breaks_new, 2)[1], breaks_new)]
    }

  }else
  {
    breaks_new <- c(0, breaks_new, length(timeseries))
    breaks_new <- unique(sort(breaks_new))
  }

  return(breaks_new)
}

#' Automatic Anomaly detection
#'
#' @param timeseries Given time series
#' @param frequency Timeseries frequency, defaults to 12 points
#' @param conf_level Confidence level for Anomaly detection
#' @param breaks breakpoints identified
#' @param window_len Window length for anomaly detection
#' @param window_len Window length for anomaly detection
#' @examples
#' AnomalyDetection(timeseries = StructuralDecompose::Nile_dataset[,1], breaks = c(4, 50, 80))
#'
#' AnomalyDetection(timeseries = runif(n = 50, min = 1, max = 10),  breaks = c(4, 20, 30))
#' @return the list of anomalies in the time series, along with the time series plot
#' @importFrom stats mad
#' @importFrom stats median

#' @export
#'

AnomalyDetection <- function(timeseries, frequency = 52, conf_level = 1.5, breaks, window_len = 14)
{
  #Initialization
  y_med <- timeseries
  temp1 <- timeseries

  window_medians <- vector()


  #Updating the breakpoints to include the first and last digits
  breaks <- unique(sort(c(breaks, 0, length(timeseries))))

  #Taking the median series
  for(i in seq(from = 1, to = c(length(breaks) - 1), by = 1))
  {
    num_windows <- ceiling(length(timeseries[c(breaks[i]+1) : c(breaks[i+1])])/window_len)

    for(w in seq(num_windows))
    {
      if(c(breaks[i] + c(window_len * c(w))) >= breaks[i+1])
      {
        window_len_new <- breaks[i+1] - c(breaks[i] + c(window_len * c(w-1)))

        med_1 <- median(y_med[c(c(breaks[i] + 1) + c(window_len * c(w-1))) : c(breaks[i+1])])
        y_med[c(c(breaks[i] + 1) + c(window_len * c(w-1))) : c(breaks[i+1])] <- med_1

        #Storing the medians
        window_medians <- c(window_medians, rep(med_1, window_len_new))
        window_len_new <- c()
      }else
      {
        med_1 <- median(y_med[c(c(breaks[i] + 1) + c(window_len * c(w-1))) : c(c(breaks[i] + 1) + c(window_len * c(w)))])
        y_med[c(c(breaks[i] + 1) + c(window_len * c(w-1))) : c(c(breaks[i] + 1) + c(window_len * c(w)))] <- med_1

        #Storing the medians
        window_medians <- c(window_medians, rep(med_1, window_len))
      }
    }
  }

  #Subtracting the median values
  temp1 <- timeseries - y_med

  #Extracting Anomalies
  Anomalies <- vector()
  for(i in seq(from = 1, to = c(length(breaks) - 1), by = 1))
  {
    med_1 <- mad(temp1[c(breaks[i] + 1) : c(breaks[i+1])], constant = conf_level)
    median_val <- median(c(temp1[c(breaks[i] + 1) : c(breaks[i+1])]))

    if(breaks[i+1] == length(temp1))
    {
      next
    }else
    {
      for(p in seq(1, (length(temp1[c(breaks[i] + 1) : c(breaks[i+1])]))))
      {
        Upper_bound <- median_val + (med_1 * conf_level)
        Lower_bound <- median_val - (med_1 * conf_level)

        #Upper confidence bound
        if(temp1[c(breaks[i] + 1) : c(breaks[i+1])][p] > Upper_bound)
        {
          temp1[c(breaks[i] + 1) : c(breaks[i+1])][p] <- Upper_bound
          Anomalies <- c(Anomalies, c(breaks[i] + p))
        }

        #Lower confidence bound
        if(temp1[c(breaks[i] + 1) : c(breaks[i+1])][p] < Lower_bound)
        {
          temp1[c(breaks[i] + 1) : c(breaks[i+1])][p] <- Lower_bound
          Anomalies <- c(Anomalies, c(breaks[i] + p))

        }
      }
    }

  }

  #Re-adding the median values
  timeseries <- (temp1 + window_medians)
  Anomalies <- sort(unique(Anomalies))

  Results <- list('DeAnomalized_series' = timeseries, 'Anomalies' = Anomalies)

  return(Results)

  }


#' Smoothening of the time series
#'
#' @param timeseries Given time series
#' @param frequency Timeseries frequency, defaults to 12 points
#' @param smoothening_algorithm Smoothening algorithm required
#' @param breaks Breakpoints identified by the previous algorithm
#' @param lowess Lowess smoothener
#' @examples
#' Smoothing(timeseries = StructuralDecompose::Nile_dataset[,1], breaks = c(4, 50, 80))
#'
#' Smoothing(timeseries = runif(n = 50, min = 1, max = 10), breaks = c(4, 20, 30))
#' @return The smoothened time series
#' @importFrom utils tail
#' @importFrom stats lowess
#' @importFrom utils tail

#' @export
#'
Smoothing <- function(timeseries, frequency = 52, smoothening_algorithm = 'lowess', breaks)
{
  #Smoothening the series
  k <- vector()
  k <- breaks
  k <- sort(k)

  trend_line <- vector()
  if(length(k) == 0)
  {
    smooth_series <- lowess(timeseries)
    trend_line <- c(trend_line, smooth_series$y)
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
        smooth_series <- lowess(v)$y
      }

      if(smoothening_algorithm == 'SMA')
      {
        #smooth_series = smooth::sma(v)
      }

      trend_line <- c(trend_line, smooth_series)
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
#' @param conf_level Confidence level for Anomaly detection, best to keep this a static value
#' @param window_len Length of the Moving window for Anomaly Detection
#' @param plot True of False indicating if you want the internal plots to be generated
#' @examples
#' StructuralDecompose(Data = StructuralDecompose::Nile_dataset[,1])
#'
#' StructuralDecompose(Data = runif(n = 50, min = 1, max = 10))
#' @return The decomposed time series along with a host of other metrics
#' @importFrom stats ts
#' @importFrom stats stl
#' @export
#'
StructuralDecompose <- function(Data, frequency = 12, break_algorithm = 'strucchange', smoothening_algorithm = 'lowess', break_level = 0.05, median_level = 0.5, mean_level = 0.5, level_length = 12, conf_level = 0.5, window_len = 12, plot = FALSE)
{

  #Initial Sanity checks
  if(!is.numeric(frequency)  || !is.numeric(break_level) || !is.numeric(mean_level) || !is.numeric(median_level) || !is.numeric(level_length) || !is.numeric(conf_level) || !is.numeric(window_len))
  {
    stop('Value needs to be numeric')
  }

  if(!is.logical(plot))
  {
    stop('Value needs to be boolean')
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
  Anom_output <- AnomalyDetection(timeseries = Data, frequency = frequency, breaks = Break_points, conf_level = 0.5)
  Cleanseries <- Anom_output$DeAnomalized_series
  Anomalies <- Anom_output$Anomalies

  #Smoothing of the time series
  Decomposedtrend <- Smoothing(timeseries = Cleanseries, breaks = Break_points)

  #Detrending the series
  Detrended_Data <- c(Data - Decomposedtrend)


  #Transforming into a series, final decomp
  Detrended_Data <- ts(data = as.vector(t(Detrended_Data)), frequency = frequency)


  decomposed <- NA
  tryCatch(
    {
      decomposed <- stl(Detrended_Data, s.window = 'periodic')
    }, error = function(e){decomposed <<- NULL}
  )

  if(length(decomposed) != 1)
  {
    seasonal <- decomposed$time.series[,1]
    trend <- decomposed$time.series[,2]
    remainder <- decomposed$time.series[,3]

    #Removing seasonality
    Deseasonalized <- Detrended_Data - seasonal
  }else
  {
    seasonal <- NULL
    trend <- NULL
    remainder <- NULL
  }

  newList <- list('anomalies' = Anomalies, 'trend_line' = Decomposedtrend, 'Deseaonalized_Series' = Deseasonalized,
                  'breakpoints' = Break_points, 'trend' = trend, 'seasonality' = seasonal, 'remainder' = remainder)



  return(newList)
}
