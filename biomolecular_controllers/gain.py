from dataclasses import dataclass
from typing import Dict, List, Tuple, Literal, Callable, Optional
import numpy as np

Params = Dict[str, float]
Axis = Literal["alpha_1", "alpha_2", "theta_1", "theta_2", "beta", "gamma"]  # extend as needed

@dataclass
class ClosedLoopGain:
    params: Params
    gain: float
    metrics: Dict[str, float]          # e.g. {"settling_time":..., "overshoot":...}
    stable_margin: Optional[float] = None  # e.g. max Re(eigs)

def pairwise_sweep_cloud(
    model_key: str,
    base_params: Params,
    x_axis: Axis, x_values: np.ndarray,
    y_axis: Axis, y_values: np.ndarray,
    gain_compute: Callable[[Params], float],
    simulate_and_metrics: Callable[[Params], Dict[str, float]],
    stable_margin_fn: Optional[Callable[[Params], float]] = None,
) -> List[ClosedLoopGain]:
    cloud: List[ClosedLoopGain] = []
    for y in y_values:
        for x in x_values:
            p = dict(base_params)
            p[x_axis] = float(x)
            p[y_axis] = float(y)
            g = float(gain_compute(p))
            m = simulate_and_metrics(p)
            sm = float(stable_margin_fn(p)) if stable_margin_fn else None
            cloud.append(ClosedLoopGain(params=p, gain=g, metrics=m, stable_margin=sm))
    return cloud

def bin_by_gain(
    cloud: List[ClosedLoopGain],
    metric: str,
    gain_bins: np.ndarray,
    require_stable: bool = True,
) -> Dict[str, np.ndarray]:
    gains = np.array([pt.gain for pt in cloud], dtype=float)
    vals  = np.array([pt.metrics[metric] for pt in cloud], dtype=float)

    if require_stable:
        mask = np.array([(pt.stable_margin is None) or (pt.stable_margin < 0) for pt in cloud], dtype=bool)
        gains, vals = gains[mask], vals[mask]

    idx = np.digitize(gains, gain_bins) - 1  # bin index
    nb = len(gain_bins) - 1

    x_mid = 0.5 * (gain_bins[:-1] + gain_bins[1:])
    med = np.full(nb, np.nan)
    q25 = np.full(nb, np.nan)
    q75 = np.full(nb, np.nan)
    best = np.full(nb, np.nan)  # “best” = min by default for time/error metrics

    for b in range(nb):
        inb = vals[idx == b]
        if inb.size == 0:
            continue
        med[b] = np.median(inb)
        q25[b] = np.quantile(inb, 0.25)
        q75[b] = np.quantile(inb, 0.75)
        best[b] = np.min(inb)

    return {"gain_mid": x_mid, "median": med, "q25": q25, "q75": q75, "best": best}

# gain.py (append)

def gain_pc(p: Params) -> float:
    """P+LPF: closed loop gain = alpha_1*theta_1*theta_2"""
    return p["alpha_1"] * p["theta_1"] * p["theta_2"]

def gain_pd1(p: Params) -> float:
    """PD1 (crossed): KP = (alpha_1 - beta*alpha_2)*theta_1*theta_2"""
    return (p["alpha_1"] - p["beta"] * p["alpha_2"]) * p["theta_1"] * p["theta_2"]

def gain_pd2(p: Params) -> float:
    """PD2 (uncrossed): KP = (alpha_1 + beta*alpha_2)*theta_1*theta_2"""
    return (p["alpha_1"] + p["beta"] * p["alpha_2"]) * p["theta_1"] * p["theta_2"]

def gain_pi1(p: Params) -> float:
    """PI1: proportional numerator = (alpha_1 - alpha_2*beta)*theta_1*theta_2
    DC gain = KP / (beta*k - 1), zero SS error when beta*k = 1"""
    KP = (p["alpha_1"] - (p["alpha_2"] * p["beta"])) * p["theta_1"] * p["theta_2"]
    bk = p["beta"] * p["k"]
    return KP / (bk - 1) if abs(bk - 1) > 1e-9 else np.inf

def gain_pi2(p: Params) -> float:
    """PI2: proportional numerator = (alpha_1 + alpha_2*beta)*theta_1*theta_2
    DC gain = KP / (beta*k - 1), zero SS error when beta*k = 1"""
    KP = (p["alpha_1"] + p["alpha_2"] * p["beta"]) * p["theta_1"] * p["theta_2"]
    bk = p["beta"] * p["k"]
    return KP / (bk - 1) if abs(bk - 1) > 1e-9 else np.inf

def gain_pid1(p: Params) -> float:
    """PID1: KI = (alpha_1 - alpha_2*beta*rho)*theta_1*theta_2
    DC gain = KI / (beta*k - 1), zero SS error when beta*k = 1"""
    KI = (p["alpha_1"] - p["alpha_2"] * p["beta"] * p["rho"]) * p["theta_1"] * p["theta_2"]
    bk = p["beta"] * p["k"]
    return KI / (bk - 1) if abs(bk - 1) > 1e-9 else np.inf

def gain_pid2(p: Params) -> float:
    """PID2: KI = (alpha_1 + alpha_2*beta*rho)*theta_1*theta_2
    DC gain = KI / (beta*k - 1), zero SS error when beta*k = 1"""
    KI = (p["alpha_1"] + p["alpha_2"] * p["beta"] * p["rho"]) * p["theta_1"] * p["theta_2"]
    bk = p["beta"] * p["k"]
    return KI / (bk - 1) if abs(bk - 1) > 1e-9 else np.inf

GAIN_FNS: Dict[str, Callable[[Params], float]] = {
    "PC":    gain_pc,
    "PI_1":  gain_pi1,
    "PI_2":  gain_pi2,
    "PD_1":  gain_pd1,
    "PD_2":  gain_pd2,
    "PID_1": gain_pid1,
    "PID_2": gain_pid2,
}

def get_gain_fn(model_key: str) -> Callable[[Params], float]:
    try:
        return GAIN_FNS[model_key]
    except KeyError as e:
        raise KeyError(f"No gain function registered for model_key={model_key}") from e


# ── Analytical stability boundaries ──────────────────────────────────
#
# Each returns a callable ``f(alpha_1_array) -> theta_1_array`` that
# traces the stability boundary in (α₁, θ₁) parameter space.
#
# Derivation summary (Routh–Hurwitz):
#   PC / PD  — 3rd-order characteristic polynomial → G_crit = 8
#   PI / PID — 4th/5th-order with integrator pole  → G_crit = 6
#
# In all cases the boundary is the hyperbola  θ₁ = G_crit / (G_eff · θ₂)
# where G_eff depends on the controller variant.

def pc_stability_boundary(
    theta_2: float,
    G_crit: float = 8.0,
) -> Callable[[np.ndarray], np.ndarray]:
    """
    Stability boundary for the PC controller in (α₁, θ₁) space.

    G_eff = α₁  →  θ₁_crit = G_crit / (α₁ · θ₂)
    """
    def _boundary(alpha_1: np.ndarray) -> np.ndarray:
        return G_crit / (np.asarray(alpha_1, dtype=float) * theta_2)
    return _boundary


def pd1_stability_boundary(
    theta_2: float,
    alpha_2: float,
    beta: float,
    G_crit: float = 8.0,
) -> Callable[[np.ndarray], np.ndarray]:
    """
    Stability boundary for the PD_1 (crossed) controller in (α₁, θ₁) space.

    G_eff = α₁ − β·α₂  →  θ₁_crit = G_crit / ((α₁ − β·α₂) · θ₂)
    """
    def _boundary(alpha_1: np.ndarray) -> np.ndarray:
        alpha_1 = np.asarray(alpha_1, dtype=float)
        G_eff = alpha_1 - beta * alpha_2
        with np.errstate(divide='ignore', invalid='ignore'):
            result = np.where(G_eff > 0, G_crit / (G_eff * theta_2), np.inf)
        return result
    return _boundary


def pd2_stability_boundary(
    theta_2: float,
    alpha_2: float,
    beta: float,
    G_crit: float = 8.0,
) -> Callable[[np.ndarray], np.ndarray]:
    """
    Stability boundary for the PD_2 (uncrossed) controller in (α₁, θ₁) space.

    G_eff = α₁ + β·α₂  →  θ₁_crit = G_crit / ((α₁ + β·α₂) · θ₂)
    """
    def _boundary(alpha_1: np.ndarray) -> np.ndarray:
        alpha_1 = np.asarray(alpha_1, dtype=float)
        G_eff = alpha_1 + beta * alpha_2
        with np.errstate(divide='ignore', invalid='ignore'):
            result = np.where(G_eff > 0, G_crit / (G_eff * theta_2), np.inf)
        return result
    return _boundary


def pi1_stability_boundary(
    theta_2: float,
    alpha_2: float,
    beta: float,
    k: float,
    G_crit: float = 8.0,   # same upper bound as P+LPF from e1 condition
) -> Tuple[Callable, Callable]:
    """
    Stability boundaries for PI_1 in (alpha_1, theta_1) space.
    
    Upper (e1): alpha_1 * theta_1 * theta_2 < 8
        → theta_1 < 8 / (alpha_1 * theta_2)
    
    Lower (a0): (alpha_1 - alpha_2*beta) * theta_1 * theta_2 > beta*k - 1
        → theta_1 > (beta*k - 1) / ((alpha_1 - alpha_2*beta) * theta_2)
        only binding when beta*k > 1
    """
    def upper(alpha_1: np.ndarray) -> np.ndarray:
        alpha_1 = np.asarray(alpha_1, dtype=float)
        return np.where(alpha_1 > 0, G_crit / (alpha_1 * theta_2), np.inf)

    def lower(alpha_1: np.ndarray) -> np.ndarray:
        alpha_1 = np.asarray(alpha_1, dtype=float)
        F = alpha_1 - alpha_2 * beta
        G = beta * k - 1
        with np.errstate(divide='ignore', invalid='ignore'):
            result = np.where(F > 0, G / (F * theta_2), -np.inf)
        return np.where(G > 0, result, np.zeros_like(alpha_1))

    return upper, lower


def pid_stability_boundary(
    theta_2: float,
    G_crit: float = 6.0,
) -> Callable[[np.ndarray], np.ndarray]:
    """
    Stability boundary for PID_1 / PID_2 controllers in (α₁, θ₁) space.

    Integrator pole reduces phase margin → G_crit = 6.
    G_eff = α₁  →  θ₁_crit = G_crit / (α₁ · θ₂)
    """
    def _boundary(alpha_1: np.ndarray) -> np.ndarray:
        return G_crit / (np.asarray(alpha_1, dtype=float) * theta_2)
    return _boundary