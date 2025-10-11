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
BreakPoints <- function(timeseries, frequency = 52, break_algorithm = 'strucchange', break_level = 0.01)
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

#' Seasonal Breaks level checks
#'
#' @param timeseries Given time series
#' @param tolerance Tolerence of seasonality values
#' @param min_repetitions Minimun repititions for detecting seasonality
#' @param frequency Timeseries frequency, defaults to 12 points
#' @examples
#' MedianCleaning(timeseries = StructuralDecompose::Nile_dataset[,1], breaks = c(1,4,5))
#'
#' MedianCleaning(timeseries = runif(n = 50, min = 1, max = 10), breaks = c(1,4,5))
#' @return The series cleaned with the median check
#' @export
#'
SeasonalBreaks <- function(breaks, frequency, tolerance = 0.10, min_repetitions = 3) 
{
  t <- sort(unique(as.numeric(breaks)))
  n <- length(t)
  
  if (n < min_repetitions)
  {
    return(t)
  }
    
  
  # --- Step 2: Define tolerance band (relative) ---
  freq_min <- frequency * (1 - tolerance)
  freq_max <- frequency * (1 + tolerance)
  
  # --- Step 3: Build adjacency matrix of "near-frequency" links ---
  diff_mat <- abs(outer(t, t, "-"))
  adj <- (diff_mat >= freq_min) & (diff_mat <= freq_max)
  diag(adj) <- FALSE   # no self-links
  
  visited <- rep(FALSE, n)
  groups <- list()
  
  for (i in seq_len(n)) 
  {
    if (!visited[i]) 
    {
      stack <- c(i)
      component <- c()
      
      while (length(stack) > 0) 
      {
        node <- stack[1]
        stack <- stack[-1]
        if (!visited[node]) 
        {
          visited[node] <- TRUE
          component <- c(component, node)
          neighbors <- which(adj[node, ])
          stack <- c(stack, neighbors[!visited[neighbors]])
        }
      }
      groups <- c(groups, list(component))
    }
  }
  
  # --- Step 5: Remove groups that have enough seasonal hits ---
  to_remove <- rep(FALSE, n)
  
  for (g in groups) 
  {
    if (length(g) >= min_repetitions) 
    {
      to_remove[g] <- TRUE
    }
  }
  
  # --- Step 6: Return filtered set ---
  t_filtered <- t[!to_remove]
  if (length(t_filtered) < 2)
  {
    t_filtered <- c(min(t), max(t))
  }
  
  return(t_filtered)
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

MedianCleaning <- function(timeseries, median_level = 0.1, breaks, frequency = 52, max_iter = 10000) 
{
  
  # Initialize breakpoints
  if (any(is.na(breaks))) 
  {
    t <- c(0, length(timeseries))
  } else 
  {
    t <- sort(unique(breaks))
  }
  
  # --- Seasonal repetition pruning ---
  if (frequency != 1) 
  {
    t <- SeasonalBreaks(breaks = t, frequency = frequency,  tolerance = 0.10, min_repetitions = 3)
  }
  
  # --- Recursive median-based merging ---
  iter <- 0L
  repeat 
  {
    iter <- iter + 1L
    if (iter > max_iter) stop("Max iterations reached — possible infinite loop")
    if (length(t) <= 2) break

    merged <- FALSE

    # scan triplets left->right and drop at most one boundary per pass
    for (i in seq_len(length(t) - 2L)) 
    {
      L  <- t[i]
      Lr <- t[i + 1L]
      R  <- t[i + 2L]

      # guard invalid/empty segments
      if (Lr <= L || R <= Lr)
      {
        next
      }
      
      seg1 <- timeseries[(L + 1L):Lr]
      seg2 <- timeseries[(Lr + 1L):R]
      
      if (length(seg1) == 0L || length(seg2) == 0L) 
      {
        next
      }
      
      mid1 <- median(seg1, na.rm = TRUE)
      mid2 <- median(seg2, na.rm = TRUE)

      dist <- abs(mid1 - mid2) / max(abs(mid1), abs(mid2), 1e-7) #0 is almost same and 1 is very different

      if (dist < median_level) #Higher the median, more merging happens i.e we drop more breaks
      {
        # LEFT MERGE: drop the start of the right segment
        t <- t[-(i + 1L)]
        merged <- TRUE
        break
      }
    }
    
    if (!merged) 
    {
      break
    }
  }

  return(t)
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

MeanCleaning <- function(timeseries, mean_level = 0.25, breaks,  window = 4) 
{
  n <- length(timeseries)
  
  if (length(breaks) <= 2)
  {
    return(breaks)  # no internal breaks to clean
  }
  
  breaks_new <- c(0)
  
  for (b in breaks[-c(1, length(breaks))]) 
  {  # exclude first (0) and last (n)
    
    # Define local window safely within bounds
    start_before <- max(1, b - window)
    end_before   <- max(1, b - 1)
    
    start_after  <- min(n, b + 1)
    end_after    <- min(n, b + window)
    
    before <- timeseries[start_before:end_before]
    after  <- timeseries[start_after:end_after]
    
    M1 <- mean(before, na.rm = TRUE)
    M2 <- mean(after, na.rm = TRUE)
    
    M1 <- ifelse(M1 == 0, 1e-7, M1)
    M2 <- ifelse(M2 == 0, 1e-7, M2)
    
    dist <- abs(M1 - M2) / max(abs(M1), abs(M2), 1e-7)
    
    if (dist > mean_level) 
    {
      # Keep this breakpoint
      breaks_new <- c(breaks_new, b)
    }
  }
  
  breaks_new <- c(breaks_new, n)
  return(sort(unique(breaks_new)))
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
LevelCheck <- function(timeseries, breaks, level_length = 10, max_iter = 100000)
{
  n <- length(timeseries)
  t <- sort(unique(c(0, breaks, n)))  # ensure boundaries
  
  iter <- 0L
  repeat 
  {
    iter <- iter + 1L
    if (iter > max_iter) stop("Max iterations reached — possible infinite loop")
    
    seg_lengths <- diff(t)
    k <- which(seg_lengths < level_length)[1]  # leftmost short segment
    if (is.na(k)) break  # all segments long enough
    
    # merge direction: right by default; left if it's the last segment
    drop_idx <- if (k == length(seg_lengths)) k else (k + 1L)
    
    # safety: only drop internal boundaries (not 0 or n)
    if (drop_idx <= 1L || drop_idx >= length(t)) break
    
    t <- t[-drop_idx]
  }
  
  return(t)
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

AnomalyDetection <- function(timeseries, breaks, frequency = 52, conf_level = 1.5, window_len = 14) 
{
  n <- length(timeseries)
  breaks <- sort(unique(c(0, breaks, n)))
  
  y_med <- numeric(n)
  window_medians <- numeric(n)
  
  # --- Build window-median baseline ---
  for (i in seq_len(length(breaks) - 1)) 
  {
    seg_start <- breaks[i] + 1
    seg_end   <- breaks[i + 1]
    seg <- timeseries[seg_start:seg_end]
    
    n_seg <- length(seg)
    num_windows <- ceiling(n_seg / window_len)
    
    for (w in seq_len(num_windows)) 
    {
      w_start <- seg_start + (w - 1) * window_len
      w_end   <- min(seg_start + w * window_len - 1, seg_end)
      
      med_val <- median(timeseries[w_start:w_end], na.rm = TRUE)
      y_med[w_start:w_end] <- med_val
      window_medians[w_start:w_end] <- med_val
    }
  }
  
  # --- Compute residuals ---
  resid <- timeseries - y_med
  Anomalies <- integer(0)
  
  # --- MAD-based anomaly clipping ---
  for (i in seq_len(length(breaks) - 1)) 
  {
    seg_idx <- (breaks[i] + 1):breaks[i + 1]
    seg_resid <- resid[seg_idx]
    
    seg_median <- median(seg_resid, na.rm = TRUE)
    seg_mad <- mad(seg_resid, constant = 1.4826, na.rm = TRUE)
    
    upper <- seg_median + conf_level * seg_mad
    lower <- seg_median - conf_level * seg_mad
    
    high_out <- which(seg_resid > upper)
    low_out  <- which(seg_resid < lower)
    out_idx  <- c(high_out, low_out)
    
    if (length(out_idx) > 0) 
    {
      resid[seg_idx[out_idx]] <- pmax(pmin(seg_resid[out_idx], upper), lower)
      Anomalies <- c(Anomalies, seg_idx[out_idx])
    }
  }
  
  # --- Recombine de-anomalized series ---
  de_anom <- resid + window_medians
  
  return(list(
    DeAnomalized_series = de_anom,
    Anomalies = sort(unique(Anomalies))
  ))
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
Smoothing <- function(timeseries, breaks, frequency = 52, smoothing_algorithm = "lowess", span = 0.3) 
{
  n <- length(timeseries)
  breaks <- sort(unique(c(0, breaks, n)))
  trend_line <- numeric(0)
  
  for (i in seq_len(length(breaks) - 1)) 
  {
    seg_start <- breaks[i] + 1
    seg_end   <- breaks[i + 1]
    seg <- timeseries[seg_start:seg_end]
    
    # Handle very short segments gracefully
    if (length(seg) < 3) 
    {
      smooth_seg <- rep(mean(seg, na.rm = TRUE), length(seg))
    } else if (tolower(smoothing_algorithm) == "lowess") 
    {
      smooth_seg <- lowess(seq_along(seg), seg, f = span)$y
    } else if (tolower(smoothing_algorithm) == "sma") 
    {
      if (!requireNamespace("TTR", quietly = TRUE)) 
      {
        stop("Package 'TTR' required for SMA smoothing.")
      }
      smooth_seg <- TTR::SMA(seg, n = floor(length(seg) / 5))
      smooth_seg[is.na(smooth_seg)] <- seg[is.na(smooth_seg)] # fill edges
    } else if (tolower(smoothing_algorithm) == "loess") 
    {
      smooth_seg <- predict(loess(seg ~ seq_along(seg), span = span))
    } else 
    {
      stop("Unknown smoothing algorithm. Use 'lowess', 'SMA', or 'loess'.")
    }
    
    trend_line <- c(trend_line, smooth_seg)
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
StructuralDecompose <- function(Data, frequency = 12, break_algorithm = 'strucchange', smoothening_algorithm = 'lowess', break_level = 0.01, median_level = 0.1, mean_level = 0.25, level_length = 12, conf_level = 2, window_len = 12, plot = FALSE)
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
  Break_points <- MedianCleaning(timeseries = Data, breaks = Break_points, frequency = frequency, median_level = median_level)
  
  #Mean Cleaning
  Break_points <- MeanCleaning(timeseries = Data, breaks = Break_points, mean_level = mean_level)
  
  #Level check
  Break_points <- LevelCheck(timeseries = Data, breaks = Break_points, level_length = level_length)
  
  #Anomaly Detection
  Anom_output <- AnomalyDetection(timeseries = Data, frequency = frequency, breaks = Break_points, conf_level = conf_level)
  Cleanseries <- Anom_output$DeAnomalized_series
  Anomalies <- Anom_output$Anomalies
  
  #Smoothing of the time series
  Decomposedtrend <- Smoothing(timeseries = Cleanseries, breaks = Break_points, frequency = frequency)
  
  #Detrending the series
  Detrended_Data <- ts(Data - Decomposedtrend, frequency = frequency) 
  
  # --- Step X: Seasonal decomposition (robust to annual data or short series) ---
  
  seasonal <- rep(0, length(Detrended_Data))
  trend <- rep(NA, length(Detrended_Data))
  remainder <- rep(NA, length(Detrended_Data))
  Deseasonalized <- Detrended_Data
  
  if (frequency > 1 && length(Detrended_Data) >= 2 * frequency) 
  {
    # Try STL decomposition safely
    decomposed <- tryCatch(
      stl(ts(Detrended_Data, frequency = frequency), s.window = "periodic"),
      error = function(e) NULL
    )
    
    if (!is.null(decomposed) && !is.null(decomposed$time.series)) 
    {
      seasonal  <- decomposed$time.series[, "seasonal"]
      trend     <- decomposed$time.series[, "trend"]
      remainder <- decomposed$time.series[, "remainder"]
      
      # Remove seasonal component
      Deseasonalized <- Detrended_Data - seasonal
    } else 
    {
      warning("STL decomposition failed — using raw detrended data as deseasonalized series.")
    }
    
  } else 
  {
    # No valid seasonality for annual or short data
    warning("Frequency ≤ 1 or series too short — skipping STL decomposition.")
  }
  
  # --- Step Y: Construct output ---
  newList <- list(
    anomalies            = Anomalies,
    trend_line           = Decomposedtrend,
    Deseasonalized_Series = Deseasonalized,
    breakpoints          = Break_points,
    trend                = trend,
    seasonality          = seasonal,
    remainder            = remainder
  )
  
  return(newList)
}

