"""
Stability analysis for biomolecular controllers.

Performs eigenvalue analysis at steady states, parameter sweeps,
and bifurcation detection (saddle-node, Hopf).
"""

from typing import Dict, List, Optional, TypedDict
import numpy as np
import numpy.typing as npt
from dataclasses import dataclass

from .model_library import Models, DEFAULT_PARAMS
from .gain import get_gain_fn


@dataclass
class StabilityResult:
    """Results from stability analysis at a single parameter point."""
    eigenvalues: npt.NDArray[np.complex128]  # Complex eigenvalues
    max_real: float  # Maximum real part
    stable: bool  # True if all Re(λ) < 0
    has_complex: bool  # True if any complex eigenvalues exist
    gain: float = 0.0  # Closed-loop gain computed from gain.py
    bifurcation_type: Optional[str] = None  # 'saddle-node', 'hopf', or None


@dataclass
class BifurcationPoint:
    """Detected bifurcation point."""
    param_value: float  # Parameter value where bifurcation occurs
    param_name: str  # Which parameter was varied
    bifurcation_type: str  # 'saddle-node', 'hopf', or 'unknown'
    eigenvalues: np.ndarray  # Eigenvalues at bifurcation
    max_real: float  # max(Re(λ)) at bifurcation

@dataclass
class ParameterSweepResult(TypedDict):
    param_values: npt.NDArray[np.float64]
    max_real: npt.NDArray[np.float64]
    stable: npt.NDArray[np.bool_]
    has_complex: npt.NDArray[np.bool_]
    eigenvalues: list[npt.NDArray[np.complex128]]


class StabilityAnalyzer:
    """
    Perform linear stability analysis on controller models.
    
    Uses Tellurium's Jacobian computation at steady state to determine
    stability and detect bifurcations.
    
    Examples
    --------
    >>> analyzer = StabilityAnalyzer()
    >>> 
    >>> # Single point analysis
    >>> result = analyzer.analyze_stability("PC", params={"g": 0.1})
    >>> print(f"Stable: {result.stable}, max(Re(λ)): {result.max_real}")
    >>> 
    >>> # Parameter sweep
    >>> sweep = analyzer.parameter_sweep(
    ...     "PC", 
    ...     param_name="g",
    ...     param_range=(0.01, 1.0),
    ...     n_points=50
    ... )
    """
    
    def __init__(self):
        self.factory = Models()
    
    def analyze_stability(
        self,
        model: str,
        params: Optional[Dict[str, float]] = None,
        steady_state_time: float = 1000.0,
        tolerance: float = 1e-6,
    ) -> StabilityResult:
        """
        Analyze stability at steady state.
        
        Parameters
        ----------
        model : str
            Model name
        params : dict, optional
            Parameter values
        steady_state_time : float
            How long to simulate to reach steady state
        tolerance : float
            Threshold for determining complex vs real eigenvalues
            
        Returns
        -------
        StabilityResult
            Eigenvalues and stability classification
        """
        # Create model
        rr = self.factory.create_roadrunner(model, params=params)

        # Simulate to steady state
        rr.simulate(0, steady_state_time, 10)

        # Get full Jacobian at current state
        # Tellurium returns Jacobian as 2D array
        jacobian = rr.getFullJacobian()

        # Compute eigenvalues
        eigenvalues = np.linalg.eigvals(jacobian)

        # Extract real parts
        real_parts = eigenvalues.real
        max_real = float(np.max(real_parts))

        # Check stability (all Re(λ) < 0)
        stable = max_real < -tolerance

        # Check for complex eigenvalues (indicating oscillations)
        has_complex = bool(np.any(np.abs(eigenvalues.imag) > tolerance))

        # Compute closed-loop gain: merge defaults with any overrides so the
        # gain function always has all required keys.
        full_params = {**DEFAULT_PARAMS.get(model, {}), **(params or {})}
        gain = get_gain_fn(model)(full_params)

        return StabilityResult(
            eigenvalues=eigenvalues,
            max_real=max_real,
            stable=stable,
            has_complex=has_complex,
            gain=gain,
        )
    
    def parameter_sweep(
        self,
        model: str,
        param_name: str,
        param_range: tuple[float, float],
        n_points: int = 50,
        log_space: bool = True,
        fixed_params: Optional[Dict[str, float]] = None,
    ) -> ParameterSweepResult:
        """
        Sweep one parameter and analyze stability at each point.
        
        Parameters
        ----------
        model : str
            Model name
        param_name : str
            Parameter to vary
        param_range : tuple of (min, max)
            Range to sweep
        n_points : int
            Number of points in sweep
        log_space : bool
            Use logarithmic spacing (recommended for rate constants)
        fixed_params : dict, optional
            Other parameters to fix at specific values
            
        Returns
        -------
        dict
            {
                'param_values': np.ndarray,
                'max_real': np.ndarray (max real part at each point),
                'stable': np.ndarray (bool array),
                'has_complex': np.ndarray (bool array),
                'eigenvalues': list of np.ndarray
            }
        """
        # Generate parameter values
        if log_space:
            param_values = np.logspace(
                np.log10(param_range[0]),
                np.log10(param_range[1]),
                n_points
            )
        else:
            param_values = np.linspace(param_range[0], param_range[1], n_points)
        
        # Storage
        max_reals = np.zeros(n_points)
        stable = np.zeros(n_points, dtype=bool)
        has_complex = np.zeros(n_points, dtype=bool)
        all_eigenvalues = []
        
        # Sweep
        for i, param_val in enumerate(param_values):
            # Build parameter dict
            params = fixed_params.copy() if fixed_params else {}
            params[param_name] = param_val
            
            # Analyze
            try:
                result = self.analyze_stability(model, params=params)
                max_reals[i] = result.max_real
                stable[i] = result.stable
                has_complex[i] = result.has_complex
                all_eigenvalues.append(result.eigenvalues)
            except Exception as e:
                # If simulation fails (e.g., blows up), mark as unstable
                max_reals[i] = np.inf
                stable[i] = False
                has_complex[i] = False
                all_eigenvalues.append(np.array([]))
        
        return {
            'param_values': param_values,
            'max_real': max_reals,
            'stable': stable,
            'has_complex': has_complex,
            'eigenvalues': all_eigenvalues,
        }
    
    def detect_bifurcations(
        self,
        sweep_results: Dict[str, np.ndarray],
        param_name: str,
    ) -> List[BifurcationPoint]:
        """
        Detect bifurcation points from parameter sweep results.
        
        Identifies where stability changes (max(Re(λ)) crosses zero) and
        classifies bifurcation type.
        
        Parameters
        ----------
        sweep_results : dict
            Output from parameter_sweep()
        param_name : str
            Name of parameter that was swept
            
        Returns
        -------
        list of BifurcationPoint
            Detected bifurcations with classification
        """
        param_values = sweep_results['param_values']
        max_reals = sweep_results['max_real']
        has_complex_arr = sweep_results['has_complex']
        all_eigenvalues = sweep_results['eigenvalues']
        
        bifurcations = []
        
        # Look for sign changes in max_real
        for i in range(len(max_reals) - 1):
            # Skip if either point is invalid
            if np.isinf(max_reals[i]) or np.isinf(max_reals[i+1]):
                continue
            
            # Check for zero crossing
            if max_reals[i] * max_reals[i+1] < 0:
                # Bifurcation occurred between i and i+1
                
                # Interpolate to find approximate crossing point
                # Linear interpolation: λ(p) ≈ λ₁ + (λ₂-λ₁)*(p-p₁)/(p₂-p₁)
                # Solve for p where λ(p) = 0
                p1, p2 = param_values[i], param_values[i+1]
                r1, r2 = max_reals[i], max_reals[i+1]
                
                # Bifurcation parameter value (linear interp)
                bif_param = p1 - r1 * (p2 - p1) / (r2 - r1)
                
                # Determine bifurcation type
                # If eigenvalues at crossing have imaginary parts → Hopf
                # If eigenvalues are real → Saddle-node
                
                # Check eigenvalues at both sides
                complex_before = has_complex_arr[i]
                complex_after = has_complex_arr[i+1]
                
                if complex_before or complex_after:
                    # Complex eigenvalues involved → likely Hopf
                    # (stable spiral → unstable spiral → limit cycle)
                    bif_type = "hopf"
                else:
                    # Real eigenvalues → saddle-node
                    # (node appears/disappears)
                    bif_type = "saddle-node"
                
                # Use eigenvalues from the unstable side (i+1 if going unstable)
                if max_reals[i] < 0 and max_reals[i+1] > 0:
                    eigs = all_eigenvalues[i+1]
                else:
                    eigs = all_eigenvalues[i]
                
                bifurcations.append(BifurcationPoint(
                    param_value=bif_param,
                    param_name=param_name,
                    bifurcation_type=bif_type,
                    eigenvalues=eigs,
                    max_real=0.0,  # At bifurcation, max(Re(λ)) ≈ 0
                ))
        
        return bifurcations
    
    def two_parameter_sweep(
        self,
        model: str,
        param1_name: str,
        param1_range: tuple[float, float],
        param2_name: str,
        param2_range: tuple[float, float],
        n_points: tuple[int, int] = (20, 20),
        log_space: tuple[bool, bool] = (True, True),
        fixed_params: Optional[Dict[str, float]] = None,
    ) -> Dict[str, np.ndarray]:
        """
        2D parameter sweep for stability analysis.
        
        Parameters
        ----------
        model : str
            Model name
        param1_name, param2_name : str
            Parameters to vary
        param1_range, param2_range : tuple of (min, max)
            Ranges for each parameter
        n_points : tuple of (n1, n2)
            Number of points for each parameter
        log_space : tuple of (bool, bool)
            Use log spacing for each parameter
        fixed_params : dict, optional
            Other parameters to fix
            
        Returns
        -------
        dict
            {
                'param1_values': np.ndarray (1D),
                'param2_values': np.ndarray (1D),
                'max_real_grid': np.ndarray (2D, shape (n2, n1)),
                'stable_grid': np.ndarray (2D bool),
            }
        """
        # Generate parameter grids
        if log_space[0]:
            param1_vals = np.logspace(
                np.log10(param1_range[0]),
                np.log10(param1_range[1]),
                n_points[0]
            )
        else:
            param1_vals = np.linspace(param1_range[0], param1_range[1], n_points[0])
        
        if log_space[1]:
            param2_vals = np.logspace(
                np.log10(param2_range[0]),
                np.log10(param2_range[1]),
                n_points[1]
            )
        else:
            param2_vals = np.linspace(param2_range[0], param2_range[1], n_points[1])
        
        # Storage (2D grids)
        max_real_grid = np.zeros((n_points[1], n_points[0]))
        stable_grid = np.zeros((n_points[1], n_points[0]), dtype=bool)
        
        # Sweep both parameters
        for i, p2_val in enumerate(param2_vals):
            for j, p1_val in enumerate(param1_vals):
                params = fixed_params.copy() if fixed_params else {}
                params[param1_name] = p1_val
                params[param2_name] = p2_val
                
                try:
                    result = self.analyze_stability(model, params=params)
                    max_real_grid[i, j] = result.max_real
                    stable_grid[i, j] = result.stable
                except Exception:
                    max_real_grid[i, j] = np.inf
                    stable_grid[i, j] = False
        
        return {
            'param1_values': param1_vals,
            'param2_values': param2_vals,
            'max_real_grid': max_real_grid,
            'stable_grid': stable_grid,
        }
    
    def get_steady_state(
        self,
        model: str,
        params: Optional[Dict[str, float]] = None,
        steady_state_time: float = 1000.0,
    ) -> Dict[str, float]:
        """
        Get steady-state values for all species.
        
        Parameters
        ----------
        model : str
            Model name
        params : dict, optional
            Parameter values
        steady_state_time : float
            How long to simulate
            
        Returns
        -------
        dict
            Steady-state values for each species
        """
        rr = self.factory.create_roadrunner(model, params=params)
        
        # Simulate to steady state
        result = rr.simulate(0, steady_state_time, 10)
        
        # Extract final values
        state_names = self.factory.get_state_names(model)
        steady_state = {}
        
        for i, name in enumerate(state_names):
            steady_state[name] = float(result[-1, i + 1])  # +1 for time column
        
        return steady_state
    
    def bifurcation_diagram(
        self,
        model: str,
        param_name: str,
        param_range: tuple[float, float],
        n_points: int = 100,
        log_space: bool = True,
        state_variable: str = "y",
        fixed_params: Optional[Dict[str, float]] = None,
    ) -> Dict[str, np.ndarray]:
        """
        Generate bifurcation diagram (steady state vs parameter).
        
        Shows how steady-state values change as parameter is varied,
        useful for visualizing bifurcations.
        
        Parameters
        ----------
        model : str
            Model name
        param_name : str
            Parameter to vary
        param_range : tuple of (min, max)
            Parameter range
        n_points : int
            Number of points
        log_space : bool
            Use log spacing
        state_variable : str
            Which state to track (default: output "y")
        fixed_params : dict, optional
            Fixed parameters
            
        Returns
        -------
        dict
            {
                'param_values': np.ndarray,
                'steady_states': np.ndarray,
                'stable': np.ndarray (bool),
            }
        """
        # Generate parameter values
        if log_space:
            param_values = np.logspace(
                np.log10(param_range[0]),
                np.log10(param_range[1]),
                n_points
            )
        else:
            param_values = np.linspace(param_range[0], param_range[1], n_points)
        
        steady_states = np.zeros(n_points)
        stable = np.zeros(n_points, dtype=bool)
        
        for i, param_val in enumerate(param_values):
            params = fixed_params.copy() if fixed_params else {}
            params[param_name] = param_val
            
            try:
                # Get steady state
                ss = self.get_steady_state(model, params=params)
                steady_states[i] = ss[state_variable]
                
                # Check stability
                stab_result = self.analyze_stability(model, params=params)
                stable[i] = stab_result.stable
            except Exception:
                steady_states[i] = np.nan
                stable[i] = False
        
        return {
            'param_values': param_values,
            'steady_states': steady_states,
            'stable': stable,
        }
