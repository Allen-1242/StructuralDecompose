"""
StructuralDecompose - Python port
Decompose a level-shifted time series into structural breaks,
trend, seasonality, and anomalies.
"""

import warnings
import numpy as np
from scipy.signal import savgol_filter
from scipy.ndimage import uniform_filter1d
from statsmodels.nonparametric.smoothers_lowess import lowess as sm_lowess

try:
    import ruptures as rpt
    HAS_RUPTURES = True
except ImportError:
    HAS_RUPTURES = False


# ---------------------------------------------------------------------------
# Breakpoint detection
# ---------------------------------------------------------------------------

def break_points(timeseries, frequency=12, algorithm="binseg",
                 min_size=None, n_bkps=None, pen=None):
    """Detect structural breakpoints using ruptures.

    Parameters
    ----------
    timeseries : array-like
        1-D numeric time series.
    frequency : int
        Series frequency (used as default min_size).
    algorithm : str
        One of 'binseg', 'pelt', 'dynp', 'kernel'.
    min_size : int or None
        Minimum segment length. Defaults to ``max(2, frequency // 4)``.
    n_bkps : int or None
        Expected number of breakpoints (used by binseg / dynp).
        If None and algorithm needs it, it is estimated internally.
    pen : float or None
        Penalty for PELT. If None a sensible default is chosen.

    Returns
    -------
    breaks : np.ndarray
        Sorted array of breakpoints including 0 and len(timeseries).
    """
    if not HAS_RUPTURES:
        raise ImportError("ruptures is required: pip install ruptures")

    y = np.asarray(timeseries, dtype=float)
    n = len(y)

    if min_size is None:
        min_size = max(2, frequency // 4)

    algo = algorithm.lower()

    try:
        if algo == "pelt":
            if pen is None:
                pen = np.log(n) * np.var(y)
            model = rpt.Pelt(model="l2", min_size=min_size).fit(y)
            bkps = model.predict(pen=pen)
        elif algo == "binseg":
            model = rpt.Binseg(model="l2", min_size=min_size).fit(y)
            if n_bkps is not None:
                bkps = model.predict(n_bkps=n_bkps)
            else:
                if pen is None:
                    pen = np.log(n) * np.var(y)
                bkps = model.predict(pen=pen)
        elif algo == "dynp":
            if n_bkps is None:
                n_bkps = max(1, n // (4 * frequency))
            model = rpt.Dynp(model="l2", min_size=min_size).fit(y)
            bkps = model.predict(n_bkps=n_bkps)
        elif algo == "kernel":
            model = rpt.KernelCPD(min_size=min_size).fit(y)
            if n_bkps is not None:
                bkps = model.predict(n_bkps=n_bkps)
            else:
                if pen is None:
                    pen = np.log(n) * np.var(y)
                bkps = model.predict(pen=pen)
        else:
            raise ValueError(
                f"Unknown algorithm '{algorithm}'. "
                "Use 'binseg', 'pelt', 'dynp', or 'kernel'."
            )
    except Exception:
        bkps = [n]

    # ruptures returns breakpoints with the last being n (end of series)
    # convert to our convention: sorted unique with 0 and n as boundaries
    bp = [b for b in bkps if 0 < b < n]
    breaks = sorted(set([0] + bp + [n]))
    return np.array(breaks, dtype=int)


# ---------------------------------------------------------------------------
# Seasonal breakpoint pruning
# ---------------------------------------------------------------------------

def seasonal_breaks(breaks, frequency, tolerance=0.10, min_repetitions=3):
    """Remove breakpoints that recur at approximately the seasonal frequency.

    Parameters
    ----------
    breaks : array-like
        Candidate breakpoints including boundaries (0 and n).
    frequency : int
        Seasonal frequency.
    tolerance : float
        Relative tolerance around *frequency* for matching.
    min_repetitions : int
        Minimum group size to classify as seasonal.

    Returns
    -------
    np.ndarray
        Filtered breakpoints.
    """
    brks = np.unique(np.asarray(breaks, dtype=int))
    if len(brks) <= 2:
        return brks

    bounds = np.array([brks[0], brks[-1]])
    interior = brks[1:-1]
    m = len(interior)

    if m < min_repetitions:
        return brks

    freq_min = frequency * (1 - tolerance)
    freq_max = frequency * (1 + tolerance)

    # adjacency: two interior breaks are linked if their distance is ~frequency
    diff_mat = np.abs(interior[:, None] - interior[None, :]).astype(float)
    adj = (diff_mat >= freq_min) & (diff_mat <= freq_max)
    np.fill_diagonal(adj, False)

    # connected components via BFS
    visited = np.zeros(m, dtype=bool)
    to_remove = np.zeros(m, dtype=bool)

    for i in range(m):
        if visited[i]:
            continue
        stack = [i]
        component = []
        while stack:
            node = stack.pop()
            if visited[node]:
                continue
            visited[node] = True
            component.append(node)
            neighbors = np.where(adj[node])[0]
            for nb in neighbors:
                if not visited[nb]:
                    stack.append(nb)
        if len(component) >= min_repetitions:
            to_remove[component] = True

    kept = interior[~to_remove]
    return np.sort(np.unique(np.concatenate([bounds, kept])))


# ---------------------------------------------------------------------------
# Median cleaning
# ---------------------------------------------------------------------------

def median_cleaning(timeseries, breaks, frequency=12, median_level=0.1,
                    max_iter=10000):
    """Merge adjacent segments whose medians are close.

    Parameters
    ----------
    timeseries : np.ndarray
    breaks : np.ndarray
    frequency : int
    median_level : float
        Maximum relative median distance for merging.
    max_iter : int

    Returns
    -------
    np.ndarray
        Cleaned breakpoints.
    """
    y = np.asarray(timeseries, dtype=float)
    brks = np.sort(np.unique(breaks))

    # seasonal pruning first
    if frequency > 1:
        brks = seasonal_breaks(brks, frequency)

    for iteration in range(max_iter):
        if len(brks) <= 2:
            break

        merged = False
        for i in range(len(brks) - 2):
            l, lr, r = int(brks[i]), int(brks[i + 1]), int(brks[i + 2])
            if lr <= l or r <= lr:
                continue

            seg1 = y[l:lr]
            seg2 = y[lr:r]

            if len(seg1) == 0 or len(seg2) == 0:
                continue

            mid1 = np.nanmedian(seg1)
            mid2 = np.nanmedian(seg2)
            denom = max(abs(mid1), abs(mid2), 1e-7)
            dist = abs(mid1 - mid2) / denom

            if dist < median_level:
                brks = np.delete(brks, i + 1)
                merged = True
                break

        if not merged:
            break
    else:
        raise RuntimeError("Max iterations reached in median_cleaning")

    return brks


# ---------------------------------------------------------------------------
# Mean cleaning
# ---------------------------------------------------------------------------

def mean_cleaning(timeseries, breaks, mean_level=0.25, window=4):
    """Drop breakpoints where the local means on each side are similar.

    Parameters
    ----------
    timeseries : np.ndarray
    breaks : np.ndarray
    mean_level : float
        Minimum relative mean distance to retain a break.
    window : int
        Points on each side used to compute local means.

    Returns
    -------
    np.ndarray
        Cleaned breakpoints.
    """
    y = np.asarray(timeseries, dtype=float)
    n = len(y)

    if len(breaks) <= 2:
        return breaks

    new_breaks = [0]

    for b in breaks[1:-1]:  # skip boundaries
        b = int(b)
        start_before = max(0, b - window)
        end_before = max(0, b)
        start_after = min(n, b)
        end_after = min(n, b + window)

        before = y[start_before:end_before]
        after = y[start_after:end_after]

        if len(before) == 0 or len(after) == 0:
            continue

        m1 = np.nanmean(before)
        m2 = np.nanmean(after)
        denom = max(abs(m1), abs(m2), 1e-7)
        dist = abs(m1 - m2) / denom

        if dist > mean_level:
            new_breaks.append(b)

    new_breaks.append(n)
    return np.sort(np.unique(new_breaks))


# ---------------------------------------------------------------------------
# Level length check
# ---------------------------------------------------------------------------

def level_check(timeseries, breaks, level_length=10, max_iter=100000):
    """Merge segments shorter than *level_length*.

    Parameters
    ----------
    timeseries : np.ndarray
    breaks : np.ndarray
    level_length : int
    max_iter : int

    Returns
    -------
    np.ndarray
        Cleaned breakpoints.
    """
    n = len(timeseries)
    brks = np.sort(np.unique(np.concatenate([[0], breaks, [n]])))

    for _ in range(max_iter):
        seg_lengths = np.diff(brks)
        short = np.where(seg_lengths < level_length)[0]

        if len(short) == 0:
            break

        k = short[0]

        # merge right by default; left if last segment
        if k == len(seg_lengths) - 1:
            drop_idx = k
        else:
            drop_idx = k + 1

        if drop_idx <= 0 or drop_idx >= len(brks) - 1:
            break

        brks = np.delete(brks, drop_idx)
    else:
        raise RuntimeError("Max iterations reached in level_check")

    return brks


# ---------------------------------------------------------------------------
# Anomaly detection
# ---------------------------------------------------------------------------

def anomaly_detection(timeseries, breaks, conf_level=1.5, window_len=14):
    """MAD-based anomaly detection within each segment.

    Parameters
    ----------
    timeseries : np.ndarray
    breaks : np.ndarray
    conf_level : float
        Multiplier on the MAD for anomaly thresholds.
    window_len : int
        Window length for the moving-median baseline.

    Returns
    -------
    de_anomalized : np.ndarray
        Series with anomalies clipped.
    anomalies : np.ndarray
        Indices of detected anomalies.
    """
    y = np.asarray(timeseries, dtype=float).copy()
    n = len(y)
    brks = np.sort(np.unique(np.concatenate([[0], breaks, [n]])))

    y_med = np.zeros(n)

    # build piecewise window-median baseline
    for i in range(len(brks) - 1):
        seg_start = int(brks[i])
        seg_end = int(brks[i + 1])
        seg_len = seg_end - seg_start
        num_windows = int(np.ceil(seg_len / window_len))

        for w in range(num_windows):
            w_start = seg_start + w * window_len
            w_end = min(seg_start + (w + 1) * window_len, seg_end)
            y_med[w_start:w_end] = np.nanmedian(y[w_start:w_end])

    resid = y - y_med
    anomalies = []

    for i in range(len(brks) - 1):
        seg_start = int(brks[i])
        seg_end = int(brks[i + 1])
        seg_resid = resid[seg_start:seg_end]

        seg_median = np.nanmedian(seg_resid)
        seg_mad = np.nanmedian(np.abs(seg_resid - seg_median)) * 1.4826

        upper = seg_median + conf_level * seg_mad
        lower = seg_median - conf_level * seg_mad

        out_idx = np.where((seg_resid > upper) | (seg_resid < lower))[0]

        if len(out_idx) > 0:
            resid[seg_start + out_idx] = np.clip(
                seg_resid[out_idx], lower, upper
            )
            anomalies.extend((seg_start + out_idx).tolist())

    de_anomalized = resid + y_med
    return de_anomalized, np.sort(np.unique(anomalies))


# ---------------------------------------------------------------------------
# Smoothing
# ---------------------------------------------------------------------------

def smoothing(timeseries, breaks, algorithm="lowess", span=0.3):
    """Piecewise smoothing within each segment.

    Parameters
    ----------
    timeseries : np.ndarray
    breaks : np.ndarray
    algorithm : str
        One of 'lowess', 'sma', 'savgol'.
    span : float
        Smoother span for lowess.

    Returns
    -------
    np.ndarray
        Smoothed trend line.
    """
    y = np.asarray(timeseries, dtype=float)
    n = len(y)
    brks = np.sort(np.unique(np.concatenate([[0], breaks, [n]])))
    trend = np.empty(n)

    for i in range(len(brks) - 1):
        seg_start = int(brks[i])
        seg_end = int(brks[i + 1])
        seg = y[seg_start:seg_end]
        seg_len = len(seg)

        if seg_len < 3:
            trend[seg_start:seg_end] = np.nanmean(seg)
            continue

        algo = algorithm.lower()
        if algo == "lowess":
            x = np.arange(seg_len, dtype=float)
            smoothed = sm_lowess(seg, x, frac=span, return_sorted=False)
        elif algo == "sma":
            win = max(2, seg_len // 5)
            smoothed = uniform_filter1d(seg, size=win, mode="nearest")
        elif algo == "savgol":
            win = min(seg_len, max(3, seg_len // 5))
            if win % 2 == 0:
                win += 1
            smoothed = savgol_filter(seg, window_length=win, polyorder=2)
        else:
            raise ValueError(
                f"Unknown smoothing algorithm '{algorithm}'. "
                "Use 'lowess', 'sma', or 'savgol'."
            )

        trend[seg_start:seg_end] = smoothed

    return trend


# ---------------------------------------------------------------------------
# Main decomposition
# ---------------------------------------------------------------------------

def structural_decompose(data, frequency=12, break_algorithm="binseg",
                         smoothing_algorithm="lowess",
                         break_level=None, min_size=None, n_bkps=None,
                         pen=None,
                         median_level=0.1, mean_level=0.25,
                         level_length=12, conf_level=2.0,
                         window_len=12, span=0.3, periods=None):
    """Full structural decomposition pipeline.

    Parameters
    ----------
    data : array-like
        1-D numeric time series.
    frequency : int
        Seasonal frequency.
    break_algorithm : str
        Algorithm for ruptures ('binseg', 'pelt', 'dynp', 'kernel').
    smoothing_algorithm : str
        Smoothing method ('lowess', 'sma', 'savgol').
    break_level : float or None
        Deprecated, kept for R API compat. Use *pen* or *n_bkps*.
    min_size : int or None
        Minimum segment size for ruptures.
    n_bkps : int or None
        Number of breakpoints (for binseg/dynp).
    pen : float or None
        Penalty (for pelt/binseg).
    median_level : float
        Threshold for median-based merging.
    mean_level : float
        Threshold for mean-based pruning.
    level_length : int
        Minimum segment length.
    conf_level : float
        MAD multiplier for anomaly detection.
    window_len : int
        Window for anomaly baseline.
    span : float
        Lowess smoother span.
    periods : list of int or None
        Seasonal periods for MSTL. Defaults to [frequency].
        Pass multiple values for multi-seasonal series (e.g. [24, 168]).

    Returns
    -------
    dict with keys:
        breakpoints, anomalies, trend_line, de_anomalized,
        detrended, seasonality, remainder
    """
    y = np.asarray(data, dtype=float)

    if y.ndim != 1:
        raise ValueError("data must be 1-D")

    n = len(y)

    # 1. breakpoint detection
    bps = break_points(y, frequency=frequency, algorithm=break_algorithm,
                       min_size=min_size, n_bkps=n_bkps, pen=pen)

    # 2. median cleaning
    bps = median_cleaning(y, bps, frequency=frequency,
                          median_level=median_level)

    # 3. mean cleaning
    bps = mean_cleaning(y, bps, mean_level=mean_level)

    # 4. level length check
    bps = level_check(y, bps, level_length=level_length)

    # 5. anomaly detection
    clean_series, anomalies = anomaly_detection(
        y, bps, conf_level=conf_level, window_len=window_len
    )

    # 6. smoothing
    trend_line = smoothing(clean_series, bps,
                           algorithm=smoothing_algorithm, span=span)

    # 7. detrend
    detrended = y - trend_line

    # 8. seasonal decomposition via MSTL (if possible)
    seasonality = np.zeros(n)
    remainder = np.full(n, np.nan)

    if frequency > 1 and n >= 2 * frequency:
        try:
            from statsmodels.tsa.seasonal import MSTL
            import pandas as pd
            mstl_periods = periods if periods is not None else [frequency]
            min_len = 2 * max(mstl_periods)
            if n < min_len:
                warnings.warn("Series too short for MSTL; skipping.")
            else:
                mstl_result = MSTL(pd.Series(detrended), periods=mstl_periods).fit()
                seasonal_out = mstl_result.seasonal
                if hasattr(seasonal_out, 'values'):
                    seasonal_out = seasonal_out.values
                if seasonal_out.ndim == 2:
                    seasonality = seasonal_out.sum(axis=1)
                else:
                    seasonality = seasonal_out
                resid_out = mstl_result.resid
                remainder = resid_out.values if hasattr(resid_out, 'values') else resid_out
        except Exception:
            warnings.warn("MSTL decomposition failed; skipping.")
    else:
        warnings.warn(
            "Frequency <= 1 or series too short; skipping MSTL."
        )

    deseasonalized = detrended - seasonality

    return {
        "breakpoints": bps,
        "anomalies": anomalies,
        "trend_line": trend_line,
        "de_anomalized": clean_series,
        "detrended": detrended,
        "deseasonalized": deseasonalized,
        "seasonality": seasonality,
        "remainder": remainder,
    }


# ---------------------------------------------------------------------------
# statsforecast-compatible wrapper
# ---------------------------------------------------------------------------

class StructuralAutoARIMA:
    """Decompose-first forecasting model compatible with statsforecast.

    Detects structural breaks, truncates to the last regime,
    then delegates to AutoARIMA for the forecast.
    """

    def __init__(self, frequency=12, break_algorithm="binseg",
                 median_level=0.1, mean_level=0.25, level_length=12,
                 season_length=None, alias=None):
        self.frequency = frequency
        self.break_algorithm = break_algorithm
        self.median_level = median_level
        self.mean_level = mean_level
        self.level_length = level_length
        self.season_length = season_length or frequency
        self.alias = alias or "StructuralAutoARIMA"
        self.uses_exog = False

    def new(self):
        return StructuralAutoARIMA(
            frequency=self.frequency,
            break_algorithm=self.break_algorithm,
            median_level=self.median_level,
            mean_level=self.mean_level,
            level_length=self.level_length,
            season_length=self.season_length,
            alias=self.alias,
        )

    def _get_last_regime(self, y):
        """Run decomposition and return the last regime's data + level offset."""
        bps = break_points(y, frequency=self.frequency,
                           algorithm=self.break_algorithm)
        bps = median_cleaning(y, bps, frequency=self.frequency,
                              median_level=self.median_level)
        bps = mean_cleaning(y, bps, mean_level=self.mean_level)
        bps = level_check(y, bps, level_length=self.level_length)

        # last regime
        last_start = int(bps[-2])
        last_segment = y[last_start:]

        # level offset = median of last segment
        level = np.nanmedian(last_segment)

        return last_segment, level, bps

    def forecast(self, h, y, X=None, X_future=None, fitted=False, **kwargs):
        from statsforecast.models import AutoARIMA

        last_seg, level, bps = self._get_last_regime(y)

        inner = AutoARIMA(season_length=self.season_length)
        result = inner.forecast(h=h, y=last_seg, X=X_future,
                                fitted=False, **kwargs)

        out = {"mean": result["mean"]}

        # pass through prediction intervals if present
        for key in result:
            if key.startswith("lo-") or key.startswith("hi-"):
                out[key] = result[key]

        if fitted:
            # re-fit on last segment to get fitted values, pad with NaN
            inner_fitted = inner.forecast(h=h, y=last_seg, fitted=True)
            full_fitted = np.full(len(y), np.nan)
            if "fitted" in inner_fitted:
                seg_fitted = inner_fitted["fitted"]
                full_fitted[-len(seg_fitted):] = seg_fitted
            out["fitted"] = full_fitted

        return out

    def __repr__(self):
        return self.alias


# ---------------------------------------------------------------------------
# Quick demo / smoke test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    np.random.seed(42)

    # synthetic series with a level shift at t=80
    n = 200
    t = np.arange(n, dtype=float)
    y = np.where(t < 80, 10.0, 20.0) + np.random.normal(0, 1, n)

    result = structural_decompose(y, frequency=12)
    print("Breakpoints:", result["breakpoints"])
    print("Anomalies:", result["anomalies"])
    print("Trend (first 10):", result["trend_line"][:10])

    # test the wrapper
    try:
        wrapper = StructuralAutoARIMA(frequency=12)
        seg, lvl, bps = wrapper._get_last_regime(y)
        print(f"\nLast regime: {len(seg)} points, level={lvl:.2f}")
        print(f"Breaks found: {bps}")
    except Exception as e:
        print(f"Wrapper test skipped: {e}")