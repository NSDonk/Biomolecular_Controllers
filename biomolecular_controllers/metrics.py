"""
Performance metrics calculator for controller systems.

Computes standard control theory metrics:
- Overshoot
- Settling time
- Rise time
- Steady-state error
- Tracking error (MAE, RMSE)
"""

from typing import Dict, Optional, Tuple
import numpy as np


class MetricsCalculator:
    """
    Calculate performance metrics from time series data.
    
    All methods return dictionaries with numeric values only.
    Visualization is handled separately by VisualizationPipeline.
    
    Examples
    --------
    >>> calc = MetricsCalculator()
    >>> metrics = calc.overshoot(time, y, ref=1.0, pert_start=20, pert_end=100)
    >>> print(f"Overshoot: {metrics['magnitude']:.2%}")
    >>> print(f"Time to peak: {metrics['time_to_peak']:.1f}")
    """
    
    @staticmethod
    def find_settling_window(
        time: np.ndarray,
        y: np.ndarray,
        pert_start: float,
        pert_end: float,
        window_frac: float = 0.2,    # size of rolling window as fraction of post-pert time
        std_threshold: float = 0.02, # std < X * signal_range = settled
    ) -> Tuple[float, bool]:
        """
        Find the time at which the signal has settled by looking for
        a region where rolling std drops below threshold.
        
        Returns (settle_time, is_oscillating)
        """
        mask = (time >= pert_start) & (time <= pert_end)
        time_post = time[mask]
        y_post    = y[mask]
        
        if len(y_post) == 0:
            return pert_start, False
        
        n_window  = max(2, int(window_frac * len(y_post)))
        y_range   = float(np.max(y_post) - np.min(y_post))
        threshold = std_threshold * y_range if y_range > 1e-9 else 1e-9
        
        # rolling std
        rolling_std = np.array([
            np.std(y_post[max(0, i - n_window):i + 1])
            for i in range(len(y_post))
        ])
        
        # find first point where rolling std stays below threshold
        settled_mask = rolling_std < threshold
        # require it to stay settled (no large jumps after)
        for i in range(len(settled_mask)):
            if settled_mask[i] and np.all(settled_mask[i:]):
                return float(time_post[i]), False
        
        # never settled = oscillating
        return float(time_post[len(time_post)//2]), True  # use second half as best guess

    @staticmethod
    def overshoot(
        time: np.ndarray,
        y: np.ndarray,
        ref: float,
        pert_start: float,
        pert_end: float,
        window_frac: float = 0.2,
        std_threshold: float = 0.02,
    ) -> Dict[str, float]:

        settle_time, is_oscillating = MetricsCalculator.find_settling_window(
            time, y, pert_start, window_frac, std_threshold
        )
        
        # ss from post-settling window
        ss_mask  = (time >= settle_time) & (time <= pert_end)
        ss_value = float(np.mean(y[ss_mask])) if ss_mask.any() else ref
        
        # peak from full post-perturbation window
        mask       = (time >= pert_start) & (time <= pert_end)
        y_window   = y[mask]
        time_window = time[mask]
        
        if len(y_window) == 0:
            return {'magnitude': 0.0, 'peak_value': ss_value,
                    'time_to_peak': 0.0, 'peak_time': pert_start,
                    'ss_value': ss_value, 'oscillating': is_oscillating}
        
        peak_idx   = np.argmax(y_window)
        peak_value = float(y_window[peak_idx])
        peak_time  = float(time_window[peak_idx])
        
        raw       = (peak_value - ss_value) / abs(ss_value) * 100 if ss_value != 0 else 0.0
        magnitude = max(0.0, float(raw))
        
        return {
            'magnitude':    magnitude,
            'peak_value':   peak_value,
            'ss_value':     ss_value,
            'time_to_peak': peak_time - pert_start,
            'peak_time':    peak_time,
            'oscillating':  is_oscillating,
        }
    
    @staticmethod
    def rise_time(
        time: np.ndarray,
        y: np.ndarray,
        ref: float,
        pert_start: float,
        pert_end: float,
        from_pct: float = 0.1,
        to_pct: float = 0.9,
        window_frac: float = 0.2,
        std_threshold: float = 0.02,
    ) -> Dict[str, float]:

        if pert_start is None:
            pert_start = time[0]
        
        settle_time, is_oscillating = MetricsCalculator.find_settling_window(
            time, y, pert_start, window_frac, std_threshold
        )
        
        # ss from post-settling window
        ss_mask  = (time >= settle_time) & (time <= pert_end)
        ss_value = float(np.mean(y[ss_mask])) if ss_mask.any() else ref
        
        # rise measured from pert_start
        mask        = (time >= settle_time) & (time <= pert_end)
        time_window = time[mask]
        y_window    = y[mask]
        
        if len(y_window) == 0:
            return {'rise_time': np.inf, 'from_time': np.inf,
                    'to_time': np.inf, 'valid': False}
        
        initial_value = float(y_window[0])
        total_change  = ss_value - initial_value
        
        if abs(total_change) < 1e-9:
            return {'rise_time': np.inf, 'from_time': np.inf,
                    'to_time': np.inf, 'valid': False}
        
        from_threshold = initial_value + from_pct * total_change
        to_threshold   = initial_value + to_pct   * total_change
        
        from_idx = np.where(y_window >= from_threshold)[0]
        to_idx   = np.where(y_window >= to_threshold)[0]
        
        if len(from_idx) == 0 or len(to_idx) == 0:
            return {
                'rise_time': np.inf,
                'from_time': np.inf if len(from_idx) == 0
                            else float(time_window[from_idx[0]]),
                'to_time':   np.inf if len(to_idx) == 0
                            else float(time_window[to_idx[0]]),
                'ss_value':  ss_value,
                'valid':     False,
            }
        
        return {
            'rise_time': float(time_window[to_idx[0]] - time_window[from_idx[0]]),
            'from_time': float(time_window[from_idx[0]]),
            'to_time':   float(time_window[to_idx[0]]),
            'ss_value':  ss_value,
            'oscillating': is_oscillating,
            'valid':     True,
        }
    
    @staticmethod
    def steady_state(
        time: np.ndarray,
        y: np.ndarray,
        pert_start: float,
        pert_end: float,
        window_frac: float = 0.2,
        std_threshold: float = 0.02,
    ) -> Dict[str, float]:
        """
        Compute steady state robustly for both deterministic and stochastic sims.
        - Finds settling time via rolling std
        - If oscillating: ss_value = mean of second half (stochastic-friendly),
        variance captures oscillation amplitude
        - If settled: ss_value = mean of post-settling window
        """
        settle_time, is_oscillating = MetricsCalculator.find_settling_window(
            time, y, pert_start, pert_end, window_frac, std_threshold
        )
        
        mask     = (time >= settle_time) & (time <= pert_end)
        y_window = y[mask]
        
        if len(y_window) == 0:
            return {'ss_value': np.nan, 'ss_variance': np.nan,
                    'ss_std': np.nan, 'oscillating': is_oscillating,
                    'settle_time': np.nan}
        
        y_mean = float(np.mean(y_window))
        y_std  = float(np.std(y_window))
        y_var  = float(np.var(y_window))
        
        if is_oscillating:
            # peak - mean captures oscillation amplitude as effective ss error
            y_peak   = float(np.max(np.abs(y_window - y_mean)))
            ss_value = y_mean + y_peak   # worst case deviation
        else:
            ss_value = y_mean
        
        return {
            'ss_value':    ss_value,
            'ss_variance': y_var,
            'ss_std':      y_std,
            'oscillating': is_oscillating,
            'settle_time': settle_time,
        }
    
    @staticmethod
    def tracking_error(
        time: np.ndarray,
        y: np.ndarray,
        ref: np.ndarray,
    ) -> Dict[str, float]:
        """
        Calculate tracking error metrics.
        
        Parameters
        ----------
        time : np.ndarray
            Time points
        y : np.ndarray
            Output signal
        ref : np.ndarray or float
            Reference signal (array or constant)
            
        Returns
        -------
        dict
            {
                'mae': float (mean absolute error),
                'rmse': float (root mean squared error),
                'max_error': float (maximum absolute error),
                'mean_error': float (mean signed error, shows bias)
            }
        """
        # Handle constant reference
        if np.isscalar(ref):
            ref = np.full_like(y, ref)
        
        error = y - ref
        abs_error = np.abs(error)
        
        return {
            'mae': float(np.mean(abs_error)),
            'rmse': float(np.sqrt(np.mean(error**2))),
            'max_error': float(np.max(abs_error)),
            'mean_error': float(np.mean(error)),
        }
    
    @staticmethod
    def settling_time(
        time: np.ndarray,
        y: np.ndarray,
        ref: float,
        pert_start: float,
        pert_end: float,
        tolerance: float = 0.02,
        window_frac: float = 0.2,
        std_threshold: float = 0.02,
    ) -> Dict[str, float]:
        """
        Settling time: first time the signal enters and stays within
        tolerance band around ss_value.
        tolerance=0.02 means ±2% of ss_value.
        """
        if pert_start is None:
            pert_start = time[0]
        
        settle_time, is_oscillating = MetricsCalculator.find_settling_window(
            time, y, pert_start, pert_end, window_frac, std_threshold
        )
        
        ss_mask  = (time >= settle_time) & (time <= pert_end)
        ss_value = float(np.mean(y[ss_mask])) if ss_mask.any() else ref
        
        mask        = (time >= settle_time) & (time <= pert_end)
        time_window = time[mask]
        y_window    = y[mask]
        
        if len(y_window) == 0:
            return {'settling_time': np.inf, 'settled': False,
                    'ss_value': ss_value, 'oscillating': is_oscillating}
        
        band_low  = ss_value * (1 - tolerance)
        band_high = ss_value * (1 + tolerance)
        
        within_band = (y_window >= band_low) & (y_window <= band_high)
        
        # find first index where it enters and stays within band
        for i in range(len(within_band)):
            if within_band[i] and np.all(within_band[i:]):
                return {
                    'settling_time': float(time_window[i] - pert_start),
                    'settled':       True,
                    'ss_value':      ss_value,
                    'oscillating':   is_oscillating,
                }
        
        return {
            'settling_time': np.inf,
            'settled':       False,
            'ss_value':      ss_value,
            'oscillating':   is_oscillating,
        }
