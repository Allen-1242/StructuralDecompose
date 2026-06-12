# StructuralDecompose
# Novel method to decompose a level-shifted time series

# NOTE: confirm the @format and @source fields below match the shipped data object
#' Nile River Dataset
#'
#' Annual flow measurements of the river Nile, included as the example
#' series for the package.
#'
#' @format A data frame with 100 rows and 1 numeric column
#' @source Derived from the classic Nile river flow series
#' @name Nile_dataset
#' @docType data
#' @keywords datasets
NULL

#' Generation of breakpoints
#'
#' @param timeseries Given time series
#' @param frequency Time series frequency, defaults to 52 points
#' @param break_algorithm Breakpoint algorithm to be used. One of
#'   'strucchange', 'changepoint' or 'segmented'. Defaults to 'strucchange'
#' @param break_level Additional parameter for the breakpoint algorithm.
#'   For 'strucchange' this is the minimal segment size h
#'
#' @return A numeric vector of breakpoints, always including 0 and the
#'   series length
#' @examples
#' BreakPoints(timeseries = seq(100), frequency = 52, break_level = 0.05)
#' BreakPoints(timeseries = StructuralDecompose::Nile_dataset[,1],
#'             frequency = 52, break_level = 0.1)
#' @importFrom stats lm
#' @export
#'
BreakPoints <- function(timeseries, frequency = 52, break_algorithm = 'strucchange', break_level = 0.01)
{
  if (!break_algorithm %in% c('strucchange', 'changepoint', 'segmented'))
  {
    stop("Unknown break_algorithm. Use 'strucchange', 'changepoint' or 'segmented'.")
  }

  bp <- NA

  # Strucchange algorithm
  if (break_algorithm == 'strucchange')
  {
    bp <- tryCatch(
      strucchange::breakpoints(timeseries ~ 1, h = break_level)$breakpoints,
      error = function(e) NA
    )
  }

  # Changepoint
  if (break_algorithm == 'changepoint')
  {
    changepoints <- tryCatch(
      changepoint::cpt.mean(timeseries, method = "BinSeg"),
      error = function(e) NULL
    )

    if (is.null(changepoints))
    {
      bp <- NA
    } else
    {
      bp <- changepoint::cpts(changepoints)
    }

    if (length(bp) == 0 || anyNA(bp))
    {
      bp <- NA
      warning('No changepoints found; consider adjusting the minimum segment size')
    }
  }

  # Segmented
  if (break_algorithm == 'segmented')
  {
    index <- seq_along(timeseries)
    fit <- stats::lm(timeseries ~ index)

    psi <- tryCatch(
      segmented::segmented(fit, seg.Z = ~index)$psi,
      error = function(e) NULL
    )

    if (is.null(psi) || NROW(psi) == 0)
    {
      bp <- NA
      warning('No segmented breakpoints found; consider a different algorithm')
    } else
    {
      est <- if (is.matrix(psi)) psi[, "Est."] else psi["Est."]
      bp <- sort(unique(round(as.numeric(est))))
      bp <- bp[bp > 0 & bp < length(timeseries)]

      if (length(bp) == 0)
      {
        bp <- NA
      }
    }
  }

  # Writing the breakpoints
  if (any(is.na(bp)))
  {
    breaks <- c(0, length(timeseries))
  } else
  {
    breaks <- sort(unique(c(0, bp, length(timeseries))))
  }

  return(breaks)
}

#' Seasonal breakpoint pruning
#'
#' Removes breakpoints that recur at approximately the seasonal frequency,
#' since these are likely seasonal artefacts rather than structural level
#' shifts. The first and last entries of \code{breaks} are treated as series
#' boundaries and are always retained.
#'
#' @param breaks Numeric vector of candidate breakpoints, including the
#'   series boundaries (0 and the series length)
#' @param frequency Time series frequency
#' @param tolerance Relative tolerance around \code{frequency} when matching
#'   seasonal spacing
#' @param min_repetitions Minimum number of linked breakpoints required to
#'   classify a group as seasonal
#' @examples
#' SeasonalBreaks(breaks = c(0, 52, 104, 156, 200), frequency = 52)
#' SeasonalBreaks(breaks = c(0, 30, 75, 200), frequency = 52)
#' @return Numeric vector of breakpoints with seasonal repeats removed
#' @export
#'
SeasonalBreaks <- function(breaks, frequency, tolerance = 0.10, min_repetitions = 3)
{
  brks <- sort(unique(as.numeric(breaks)))
  n_brks <- length(brks)

  if (n_brks <= 2)
  {
    return(brks)
  }

  # Boundaries never participate in seasonal clustering and are always kept
  bounds <- c(brks[1], brks[n_brks])
  interior <- brks[-c(1, n_brks)]
  m <- length(interior)

  if (m < min_repetitions)
  {
    return(brks)
  }

  # Tolerance band (relative)
  freq_min <- frequency * (1 - tolerance)
  freq_max <- frequency * (1 + tolerance)

  # Adjacency matrix of near-frequency links between interior breaks
  diff_mat <- abs(outer(interior, interior, "-"))
  adj <- (diff_mat >= freq_min) & (diff_mat <= freq_max)
  diag(adj) <- FALSE

  visited <- rep(FALSE, m)
  groups <- list()

  for (i in seq_len(m))
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

  # Remove groups that have enough seasonal hits
  to_remove <- rep(FALSE, m)

  for (g in groups)
  {
    if (length(g) >= min_repetitions)
    {
      to_remove[g] <- TRUE
    }
  }

  # Return filtered set with boundaries re-attached
  return(sort(unique(c(bounds, interior[!to_remove]))))
}

#' Median level checks
#'
#' @param timeseries Given time series
#' @param median_level Maximum relative median distance between two adjacent
#'   levels for them to be merged
#' @param breaks Breaks identified
#' @param frequency Time series frequency, defaults to 52 points
#' @param max_iter Maximum number of merge iterations, used as an infinite
#'   loop guard
#' @examples
#' MedianCleaning(timeseries = StructuralDecompose::Nile_dataset[,1], breaks = c(1,4,5))
#'
#' MedianCleaning(timeseries = runif(n = 50, min = 1, max = 10), breaks = c(1,4,5))
#' @return The cleaned set of breakpoints
#' @importFrom stats median
#' @export
#'
MedianCleaning <- function(timeseries, median_level = 0.1, breaks, frequency = 52, max_iter = 10000)
{
  # Initialize breakpoints
  if (any(is.na(breaks)))
  {
    brks <- c(0, length(timeseries))
  } else
  {
    brks <- sort(unique(breaks))
  }

  # Seasonal repetition pruning
  if (frequency != 1)
  {
    brks <- SeasonalBreaks(breaks = brks, frequency = frequency, tolerance = 0.10, min_repetitions = 3)
  }

  # Recursive median-based merging
  iter <- 0L
  repeat
  {
    iter <- iter + 1L
    if (iter > max_iter) stop("Max iterations reached; possible infinite loop")
    if (length(brks) <= 2) break

    merged <- FALSE

    # Scan triplets left to right and drop at most one boundary per pass
    for (i in seq_len(length(brks) - 2L))
    {
      L  <- brks[i]
      Lr <- brks[i + 1L]
      R  <- brks[i + 2L]

      # Guard invalid or empty segments
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

      # 0 is almost the same, 1 is very different
      dist <- abs(mid1 - mid2) / max(abs(mid1), abs(mid2), 1e-7)

      # Higher median_level means more merging, i.e. more breaks dropped
      if (dist < median_level)
      {
        # LEFT MERGE: drop the start of the right segment
        brks <- brks[-(i + 1L)]
        merged <- TRUE
        break
      }
    }

    if (!merged)
    {
      break
    }
  }

  return(brks)
}

#' Mean level checks
#'
#' @param timeseries Given time series
#' @param mean_level Minimum relative mean distance between the points on
#'   either side of a breakpoint for it to be retained
#' @param breaks Breakpoints returned
#' @param window Number of points on each side of a breakpoint used to
#'   compute the local means
#' @examples
#' MeanCleaning(timeseries = StructuralDecompose::Nile_dataset[,1], breaks = c(1,4,5))
#'
#' MeanCleaning(timeseries = runif(n = 50, min = 1, max = 10), breaks = c(1,4,5))
#' @return The cleaned set of breakpoints
#' @export
#'
MeanCleaning <- function(timeseries, mean_level = 0.25, breaks, window = 4)
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
#' @param breaks Breakpoints returned
#' @param level_length Minimum number of points required for a level
#' @param max_iter Maximum number of merge iterations, used as an infinite
#'   loop guard
#' @examples
#' LevelCheck(timeseries = StructuralDecompose::Nile_dataset[,1], breaks = c(1,4,5))
#'
#' LevelCheck(timeseries = runif(n = 50, min = 1, max = 10), breaks = c(1,4,5))
#' @return The cleaned set of breakpoints
#' @export
#'
LevelCheck <- function(timeseries, breaks, level_length = 10, max_iter = 100000)
{
  n <- length(timeseries)
  brks <- sort(unique(c(0, breaks, n)))  # ensure boundaries

  iter <- 0L
  repeat
  {
    iter <- iter + 1L
    if (iter > max_iter) stop("Max iterations reached; possible infinite loop")

    seg_lengths <- diff(brks)
    k <- which(seg_lengths < level_length)[1]  # leftmost short segment
    if (is.na(k)) break  # all segments long enough

    # Merge direction: right by default; left if it is the last segment
    drop_idx <- if (k == length(seg_lengths)) k else (k + 1L)

    # Safety: only drop internal boundaries (not 0 or n)
    if (drop_idx <= 1L || drop_idx >= length(brks)) break

    brks <- brks[-drop_idx]
  }

  return(brks)
}

#' Automatic Anomaly detection
#'
#' @param timeseries Given time series
#' @param breaks Breakpoints identified
#' @param frequency Time series frequency, defaults to 52 points
#' @param conf_level Multiplier applied to the MAD when flagging anomalies
#' @param window_len Window length for the moving median baseline
#' @examples
#' AnomalyDetection(timeseries = StructuralDecompose::Nile_dataset[,1], breaks = c(4, 50, 80))
#'
#' AnomalyDetection(timeseries = runif(n = 50, min = 1, max = 10),  breaks = c(4, 20, 30))
#' @return A list with the de-anomalized series and the anomaly positions
#' @importFrom stats mad median
#' @export
#'
AnomalyDetection <- function(timeseries, breaks, frequency = 52, conf_level = 1.5, window_len = 14)
{
  n <- length(timeseries)
  breaks <- sort(unique(c(0, breaks, n)))

  y_med <- numeric(n)

  # Build window-median baseline
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

      y_med[w_start:w_end] <- median(timeseries[w_start:w_end], na.rm = TRUE)
    }
  }

  # Compute residuals
  resid <- timeseries - y_med
  Anomalies <- integer(0)

  # MAD-based anomaly clipping
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

  # Recombine de-anomalized series
  de_anom <- resid + y_med

  return(list(
    DeAnomalized_series = de_anom,
    Anomalies = sort(unique(Anomalies))
  ))
}

#' Smoothing of the time series
#'
#' @param timeseries Given time series
#' @param breaks Breakpoints identified by the previous algorithms
#' @param frequency Time series frequency, defaults to 52 points
#' @param smoothing_algorithm Smoothing algorithm to use. One of 'lowess',
#'   'SMA' or 'loess'. Defaults to 'lowess'
#' @param span Smoother span passed to the lowess and loess algorithms
#' @examples
#' Smoothing(timeseries = StructuralDecompose::Nile_dataset[,1], breaks = c(4, 50, 80))
#'
#' Smoothing(timeseries = runif(n = 50, min = 1, max = 10), breaks = c(4, 20, 30))
#' @return The smoothed trend line
#' @importFrom stats lowess loess predict
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
      ma_n <- max(2L, floor(length(seg) / 5))
      smooth_seg <- TTR::SMA(seg, n = ma_n)
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
#' @param frequency Frequency of the time series
#' @param break_algorithm Breakpoints algorithm used. Defaults to strucchange
#' @param smoothening_algorithm Smoothing algorithm used. Defaults to lowess
#' @param break_level Break level for the breakpoints algorithm
#' @param median_level Average median distance between two levels
#' @param mean_level Average mean distance between a group of points near breakpoints
#' @param level_length Minimum number of points required to determine a level
#' @param conf_level Confidence level for Anomaly detection, best to keep this a static value
#' @param window_len Length of the moving window for Anomaly Detection
#' @param plot TRUE or FALSE indicating if a summary plot of the decomposition
#'   should be drawn
#' @examples
#' StructuralDecompose(Data = StructuralDecompose::Nile_dataset[,1])
#'
#' StructuralDecompose(Data = runif(n = 50, min = 1, max = 10))
#' @return A list containing the anomalies, the trend line, the deseasonalized
#'   series, the final breakpoints and the STL trend, seasonality and
#'   remainder components
#' @importFrom stats ts stl
#' @export
#'
StructuralDecompose <- function(Data, frequency = 12, break_algorithm = 'strucchange', smoothening_algorithm = 'lowess', break_level = 0.01, median_level = 0.1, mean_level = 0.25, level_length = 12, conf_level = 2, window_len = 12, plot = FALSE)
{

  # Initial sanity checks
  if (!is.numeric(Data))
  {
    stop('Data needs to be numeric')
  }

  if(!is.numeric(frequency)  || !is.numeric(break_level) || !is.numeric(mean_level) || !is.numeric(median_level) || !is.numeric(level_length) || !is.numeric(conf_level) || !is.numeric(window_len))
  {
    stop('Value needs to be numeric')
  }

  if(!is.logical(plot))
  {
    stop('Value needs to be boolean')
  }

  # Calling the main break-point algorithm
  Break_points <- BreakPoints(timeseries = Data, frequency = frequency, break_algorithm = break_algorithm, break_level = break_level)

  # Median Cleaning
  Break_points <- MedianCleaning(timeseries = Data, breaks = Break_points, frequency = frequency, median_level = median_level)

  # Mean Cleaning
  Break_points <- MeanCleaning(timeseries = Data, breaks = Break_points, mean_level = mean_level)

  # Level check
  Break_points <- LevelCheck(timeseries = Data, breaks = Break_points, level_length = level_length)

  # Anomaly Detection
  Anom_output <- AnomalyDetection(timeseries = Data, frequency = frequency, breaks = Break_points, conf_level = conf_level, window_len = window_len)
  Cleanseries <- Anom_output$DeAnomalized_series
  Anomalies <- Anom_output$Anomalies

  # Smoothing of the time series
  Decomposedtrend <- Smoothing(timeseries = Cleanseries, breaks = Break_points, frequency = frequency, smoothing_algorithm = smoothening_algorithm)

  # Detrending the series
  Detrended_Data <- ts(Data - Decomposedtrend, frequency = frequency)

  # Seasonal decomposition (robust to annual data or short series)
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
      warning("STL decomposition failed; using raw detrended data as deseasonalized series.")
    }

  } else
  {
    # No valid seasonality for annual or short data
    warning("Frequency <= 1 or series too short; skipping STL decomposition.")
  }

  # Optional summary plot
  if (plot)
  {
    graphics::plot(as.numeric(Data), type = "l", col = "grey40",
                   xlab = "Index", ylab = "Value", main = "StructuralDecompose")
    graphics::lines(Decomposedtrend, col = "blue", lwd = 2)

    internal_breaks <- setdiff(Break_points, c(0, length(Data)))
    if (length(internal_breaks) > 0)
    {
      graphics::abline(v = internal_breaks, col = "red", lty = 2)
    }

    if (length(Anomalies) > 0)
    {
      graphics::points(Anomalies, as.numeric(Data)[Anomalies], col = "darkorange", pch = 19)
    }

    graphics::legend("topright",
                     legend = c("Series", "Trend", "Breakpoints", "Anomalies"),
                     col = c("grey40", "blue", "red", "darkorange"),
                     lty = c(1, 1, 2, NA), pch = c(NA, NA, NA, 19),
                     bty = "n", cex = 0.8)
  }

  # Construct output
  newList <- list(
    anomalies             = Anomalies,
    trend_line            = Decomposedtrend,
    Deseasonalized_Series = Deseasonalized,
    breakpoints           = Break_points,
    trend                 = trend,
    seasonality           = seasonal,
    remainder             = remainder
  )

  return(newList)
}
