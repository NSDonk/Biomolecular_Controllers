"""
Sensitivity analysis for biomolecular controllers.

Performs Sobol sensitivity analysis and parameter sweeps to quantify
how parameter variations affect controller performance metrics.
"""

from typing import Dict, List, Optional, Callable, Any
import numpy as np
import sympy as sp
from SALib.sample import saltelli
from SALib.analyze import sobol
from scipy.optimize import differential_evolution, NonlinearConstraint

from .simulation import SimulationRunner
from .metrics import MetricsCalculator


class SensitivityAnalyzer:
    """
    Perform variance-based sensitivity analysis using Sobol indices.
    
    Quantifies how much each parameter contributes to variance in
    performance metrics (overshoot, settling time, etc.).
    
    Examples
    --------
    >>> analyzer = SensitivityAnalyzer()
    >>> 
    >>> # Define parameter ranges for sensitivity
    >>> param_ranges = {
    ...     'g': (0.01, 1.0),
    ...     'delta': (0.05, 0.5),
    ...     'theta_1': (0.1, 2.0),
    ... }
    >>> 
    >>> # Analyze overshoot sensitivity
    >>> results = analyzer.sobol_analysis(
    ...     "PC",
    ...     param_ranges,
    ...     metric_function=lambda result: calc.overshoot(...)['magnitude'],
    ...     n_samples=1024
    ... )
    >>> 
    >>> print(f"Most important parameter: {results['ranking'][0]}")
    """
    
    def __init__(self):
        self.runner = SimulationRunner()
        self.calc = MetricsCalculator()
    
    def sobol_analysis(
        self,
        model: str,
        param_ranges: Dict[str, tuple[float, float]],
        metric_function: Callable[[Dict[str, float]], float],
        n_samples: int = 1024,
        fixed_params: Optional[Dict[str, float]] = None,
        calc_second_order: bool = True,
        seed: Optional[int] = None,
    ) -> Dict[str, Any]:
        """
        Perform Sobol sensitivity analysis.
        
        Uses Saltelli sampling to efficiently estimate first-order (S1) and
        total-order (ST) Sobol indices for each parameter.
        
        Parameters
        ----------
        model : str
            Model name
        param_ranges : dict
            Dictionary mapping param_name -> (min, max)
            Example: {'g': (0.01, 1.0), 'delta': (0.1, 0.5)}
        metric_function : callable
            Function that takes simulation result dict and returns a float metric
            Example: lambda r: MetricsCalculator().overshoot(r['time'], r['states']['y'], ...)['magnitude']
        n_samples : int
            Base sample size (actual: n_samples * (2*n_params + 2))
            Larger = more accurate but slower
        fixed_params : dict, optional
            Parameters to hold constant
        calc_second_order : bool
            Calculate second-order indices (interactions)
        seed : int, optional
            Random seed for reproducibility
            
        Returns
        -------
        dict
            {
                'S1': dict mapping param -> first-order index,
                'ST': dict mapping param -> total-order index,
                'S2': dict (if calc_second_order=True),
                'ranking': list of params sorted by importance (ST),
                'samples_used': int,
                'problem': SALib problem dict,
            }
        """
        # Define SALib problem
        param_names = list(param_ranges.keys())
        n_params = len(param_names)
        
        problem = {
            'num_vars': n_params,
            'names': param_names,
            'bounds': [param_ranges[name] for name in param_names]
        }
        
        # Before calling saltelli.sample, set numpy seed
        if seed is not None:
            np.random.seed(seed)

        param_values = saltelli.sample(
            problem,
            n_samples,
            calc_second_order=calc_second_order
        )
        
        n_total_samples = param_values.shape[0]
        
        # Evaluate metric for each parameter set
        Y = np.zeros(n_total_samples)
        
        for i, params_array in enumerate(param_values):
            # Build parameter dictionary
            params = fixed_params.copy() if fixed_params else {}
            for j, name in enumerate(param_names):
                params[name] = params_array[j]
            
            try:
                # Run simulation (you'll define how based on metric needs)
                # For now, assume metric_function handles simulation internally
                # This requires metric_function to be a closure that captures model/runner
                
                # Actually, let's make this work by having user pass pre-configured metric
                # OR we provide helper methods below
                
                # For flexibility, metric_function should accept (model, params, runner, calc)
                Y[i] = metric_function(params)
                
            except Exception as e:
                # If simulation fails, use NaN (will be filtered)
                Y[i] = np.nan
        
        # Remove NaN values (failed simulations)
        valid_mask = ~np.isnan(Y)
        if not np.all(valid_mask):
            print(f"Warning: {np.sum(~valid_mask)} simulations failed and were excluded")
            # SALib can't handle NaN, so we need to handle this carefully
            # For now, replace with mean (not ideal but prevents crash)
            Y[~valid_mask] = np.nanmean(Y)
        
        # Analyze Sobol indices
        Si = sobol.analyze(
            problem,
            Y,
            calc_second_order=calc_second_order,
            print_to_console=False
        )
        
        # Format results
        S1 = {name: float(Si['S1'][i]) for i, name in enumerate(param_names)}
        ST = {name: float(Si['ST'][i]) for i, name in enumerate(param_names)}
        
        # Rank parameters by total-order index (most important first)
        ranking = sorted(param_names, key=lambda p: ST[p], reverse=True)
        
        results = {
            'S1': S1,  # First-order: individual parameter effect
            'ST': ST,  # Total-order: including interactions
            'ranking': ranking,
            'samples_used': n_total_samples,
            'problem': problem,
        }
        
        if calc_second_order:
            # Second-order indices (pairwise interactions)
            S2_dict = {}
            for i, p1 in enumerate(param_names):
                for j, p2 in enumerate(param_names):
                    if i < j:  # Upper triangle only
                        key = f"{p1}*{p2}"
                        S2_dict[key] = float(Si['S2'][i, j])
            results['S2'] = S2_dict
        
        return results
    
    @staticmethod
    def get_feasible_range(conditions, sweep_param, fixed_vals):
        lower, upper = 0.0, np.inf

        for name, cond in conditions.items():
            expr = cond.subs({p: v for p, v in fixed_vals.items() 
                            if p != sweep_param})

            if sweep_param not in expr.free_symbols:
                val = float(expr)
                if val <= 0:
                    print(f"  WARNING: {name} always violated with these fixed values")
                continue

            boundary_vals = sp.solve(sp.Eq(expr, 0), sweep_param)
            
            for b in boundary_vals:
                try:
                    b_float = float(b)
                except (TypeError, ValueError):
                    continue
                    
                if b_float <= 0:
                    continue

                test_below = float(expr.subs(sweep_param, b_float * 0.99))
                test_above = float(expr.subs(sweep_param, b_float * 1.01))
                
                if test_below > 0 and test_above < 0:
                    upper = min(upper, b_float)
                    print(f"  {name}: upper bound at {b_float:.4f}")
                elif test_below < 0 and test_above > 0:
                    lower = max(lower, b_float)
                    print(f"  {name}: lower bound at {b_float:.4f}")

        return lower, upper

    @staticmethod
    def build_sweep(lower, upper, n_inside=8, sigma_frac=0.15):
        """
        sigma_frac of feasible width as one unit outside the boundary.
        If lower == 0.0 (no lower bound found), outside_low is empty and
        the sweep starts from just above 0.
        """
        has_lower = lower > 0.0
        width  = upper - lower
        sigma  = sigma_frac * width

        inside       = np.linspace(lower + 0.05*width,
                                    upper - 0.05*width, n_inside)
        outside_low  = np.array([lower - sigma]) if has_lower else np.array([])
        outside_high = np.array([upper + sigma])
        combined = np.sort(np.concatenate([outside_low, inside, outside_high]))

        return inside, outside_low, outside_high, [float(x) for x in combined]

    @staticmethod
    def make_constrained_sweep(free_param_vals, free_param_sym,
                                gain_constraints, all_params_sym,
                                fixed_params):
        known = set(fixed_params.keys()) | {free_param_sym}
        dependent = [p for p in all_params_sym if p not in known]
        
        assert len(dependent) == len(gain_constraints), \
            f"Need exactly as many constraints as dependent params, " \
            f"got {len(gain_constraints)} constraints for {dependent}"
        
        equations = [
            expr - target 
            for expr, target in gain_constraints.values()
        ]
        
        # solve symbolically
        raw_solutions = sp.solve(equations, dependent)
        
        # normalize to dict regardless of what sp.solve returns
        if isinstance(raw_solutions, list):
            if len(raw_solutions) == 0:
                raise ValueError("No symbolic solution found")
            # list of tuples case: [(val1, val2), ...]
            if isinstance(raw_solutions[0], (list, tuple)):
                raw_solutions = raw_solutions[0]
            # map to dependent symbols
            solutions = dict(zip(dependent, 
                                raw_solutions if isinstance(raw_solutions, (list, tuple)) 
                                else [raw_solutions]))
        elif isinstance(raw_solutions, dict):
            solutions = raw_solutions
        else:
            # single value for single dependent param
            solutions = {dependent[0]: raw_solutions}
        
        print("Symbolic solutions for dependent params:")
        for sym, sol in solutions.items():
            print(f"  {sym} = ", end="")
            sp.pprint(sol)
        
        results = []
        for val in free_param_vals:
            point = dict(fixed_params)
            point[free_param_sym] = val
            
            solved = {}
            for sym, sol in solutions.items():
                solved_val = float(sol.subs(point))
                if solved_val < 0:
                    print(f"  WARNING: {sym}={solved_val:.4f} < 0 at "
                        f"{free_param_sym}={val:.4f}, skipping")
                    break
                solved[sym] = solved_val
            else:
                point.update(solved)
                results.append({str(sym_key): float(sym_val) 
                    for sym_key, sym_val in point.items()})
        
        return results

    @staticmethod
    def build_optimizer(controller_type, fixed_params, rh_conditions, 
                    ss_error_expr, param_syms):
        """
        controller_type: string label
        fixed_params:    {sym: value} for truly fixed params (theta1, theta2)
        rh_conditions:   {name: sympy_expr > 0} stability constraints
        ss_error_expr:   sympy expression for e_ss in fast sequestration limit
        param_syms:      list of syms to optimize over
        """
        # convert symbolic expressions to numerical functions
        ss_error_fn = sp.lambdify(param_syms, 
                                ss_error_expr.subs(fixed_params), 
                                'numpy')
        
        constraint_fns = {
            name: sp.lambdify(param_syms,
                            cond.subs(fixed_params),
                            'numpy')
            for name, cond in rh_conditions.items()
        }
        
        # scipy constraints: each RH condition > 0
        scipy_constraints = [
            NonlinearConstraint(
                fun=lambda x, fn=fn: fn(*x),
                lb=0.0,
                ub=np.inf
            )
            for fn in constraint_fns.values()
        ]
            
        # bounds: all params positive
        bounds = [(1e-6, 10.0)] * len(param_syms)
        
        def multi_objective(x, alpha_sse=1.0, alpha_stability=0.5, alpha_gain=0.1):
            """
            Weighted sum of:
            1. steady state error (minimize)
            2. proximity to stability boundary (maximize margin)
            3. total gain magnitude (penalize very large gains)
            """
            sse    = ss_error_fn(*x)
            margin = min(fn(*x) for fn in constraint_fns.values())
            gain   = np.sum(np.array(x)**2)
            
            return (alpha_sse * sse 
                    - alpha_stability * margin  # negative: want large margin
                    + alpha_gain * gain)
        
        return multi_objective, scipy_constraints, bounds

    
    @staticmethod
    def create_metric_function(
        metric_name: str,
        t_span: tuple[float, float] = (0, 200),
        perturbations: Optional[List[Dict]] = None,
        metric_kwargs: Optional[Dict] = None,
    ) -> Callable:
        """
        Helper to create metric functions for Sobol analysis.
        
        Parameters
        ----------
        metric_name : str
            Name of metric from MetricsCalculator
            ('overshoot', 'settling_time', 'tracking_error', etc.)
        t_span : tuple
            Simulation time span
        perturbations : list of dict, optional
            Perturbations to apply
        metric_kwargs : dict, optional
            Additional kwargs for metric calculation
            
        Returns
        -------
        callable
            Function that takes (model, params, runner, calc) and returns metric value
        """
        if perturbations is None:
            perturbations = [
                {'time': t_span[1] * 0.25, 'type': 'parameter', 'param': 'ref', 'value': 1.5}
            ]
        
        if metric_kwargs is None:
            metric_kwargs = {}
        
        def metric_func(model, params, runner, calc):
            # Run simulation
            result = runner.run_with_perturbations(
                model,
                t_span=t_span,
                points=int((t_span[1] - t_span[0]) * 5),  # 5 points per time unit
                perturbations=perturbations,
                params=params,
            )
            
            # Calculate metric
            time = result['time']
            y = result['states']['y']
            
            if metric_name == 'overshoot':
                pert_time = perturbations[0]['time']
                ref = perturbations[0].get('value', 1.0)
                metric_result = calc.overshoot(
                    time, y, ref=ref,
                    pert_start=pert_time,
                    pert_end=t_span[1],
                    **metric_kwargs
                )
                return abs(metric_result['magnitude'])
            
            elif metric_name == 'settling_time':
                pert_time = perturbations[0]['time']
                ref = perturbations[0].get('value', 1.0)
                metric_result = calc.settling_time(
                    time, y, ref=ref,
                    pert_start=pert_time,
                    **metric_kwargs
                )
                return metric_result['settling_time'] if metric_result['settled'] else t_span[1]
            
            elif metric_name == 'rise_time':
                pert_time = perturbations[0]['time']
                ref = perturbations[0].get('value', 1.0)
                metric_result = calc.rise_time(
                    time, y, ref=ref,
                    pert_start=pert_time,
                    **metric_kwargs
                )
                return metric_result['rise_time'] if metric_result['valid'] else t_span[1]
            
            elif metric_name == 'steady_state_error':
                ss = calc.steady_state(time, y, window=(-1, None))
                ref = perturbations[0].get('value', 1.0) if perturbations else 1.0
                return abs(ss['ss_value'] - ref)
            
            elif metric_name == 'iae':
                ref = perturbations[0].get('value', 1.0) if perturbations else 1.0
                pert_time = perturbations[0]['time'] if perturbations else 0
                metric_result = calc.integral_metrics(
                    time, y, ref=ref,
                    pert_start=pert_time,
                    pert_end=t_span[1]
                )
                return metric_result['iae']
            
            else:
                raise ValueError(f"Unknown metric: {metric_name}")
        
        return metric_func
    
    def parameter_sweep_with_metric(
        self,
        model: str,
        param_name: str,
        param_range: tuple[float, float],
        metric_function: Callable,
        n_points: int = 50,
        log_space: bool = True,
        fixed_params: Optional[Dict[str, float]] = None,
    ) -> Dict[str, np.ndarray]:
        """
        Sweep one parameter and evaluate a metric at each point.
        
        Similar to stability parameter sweep but for arbitrary metrics.
        
        Parameters
        ----------
        model : str
            Model name
        param_name : str
            Parameter to vary
        param_range : tuple of (min, max)
            Parameter range
        metric_function : callable
            Function(model, params, runner, calc) -> float
        n_points : int
            Number of points
        log_space : bool
            Use log spacing
        fixed_params : dict, optional
            Fixed parameters
            
        Returns
        -------
        dict
            {
                'param_values': np.ndarray,
                'metric_values': np.ndarray,
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
        
        metric_values = np.zeros(n_points)
        
        for i, param_val in enumerate(param_values):
            params = fixed_params.copy() if fixed_params else {}
            params[param_name] = param_val
            
            try:
                metric_values[i] = metric_function(model, params, self.runner, self.calc)
            except Exception as e:
                metric_values[i] = np.nan
        
        return {
            'param_values': param_values,
            'metric_values': metric_values,
        }
