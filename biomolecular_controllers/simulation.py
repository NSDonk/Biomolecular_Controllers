"""
Simulation runner for biomolecular controller models.

Handles deterministic and stochastic simulations with perturbation support.
"""

from typing import Dict, List, Optional, Tuple, Union, TypedDict
from itertools import groupby
from collections import defaultdict
from numpy.typing import NDArray
import numpy as np
from tellurium import roadrunner
import tellurium as te

from .model_library import Models

class SimulationResult(TypedDict):
    time: NDArray[np.float64]
    states: dict[str, NDArray[np.float64]]

class PerturbedSimulationResult(TypedDict):
    time: NDArray[np.float64]
    states: dict[str, NDArray[np.float64]]
    perturbation_times: list[float]

class SimulationRunner:
    """
    Execute simulations with consistent interface.
    
    Supports:
    - Deterministic ODE simulation
    - Stochastic Gillespie simulation
    - Segmented perturbation handling
    
    Examples
    --------
    >>> runner = SimulationRunner()
    >>> result = runner.run_deterministic("PC", t_span=(0, 100), points=1000)
    >>> print(result['time'].shape)  # (1000,)
    >>> print(result['states']['y'])  # Output trajectory
    
    >>> # With perturbation
    >>> result = runner.run_with_perturbations(
    ...     "PI_1",
    ...     t_span=(0, 100),
    ...     perturbations=[
    ...         {'time': 20, 'type': 'step', 'species': 'r', 'value': 2.0}
    ...     ]
    ... )
    """
    
    def __init__(self):
        self.model_constructor = Models()
    
    def prepare_runner(
        self,
        model: str,
        params: Optional[Dict[str, float]] = None,
        ic: Optional[Dict[str, float]] = None,
        rr: Optional[roadrunner.ExtendedRoadRunner] = None,
    ) -> roadrunner.ExtendedRoadRunner:
        """
        Return a configured RoadRunner instance.

        If rr is None, create a fresh runner using model defaults plus any
        parameter/IC overrides. If rr is provided, reset it to baseline and
        then apply overrides in-place.
        """
        if rr is None:
            return self.model_constructor.create_roadrunner(model, params=params, ic=ic)

        rr.reset()

        if params is not None:
            for k, v in params.items():
                rr[k] = v

        if ic is not None:
            for k, v in ic.items():
                rr[k] = v

        return rr
    
    def run_deterministic(
        self,
        model: str,
        t_span: Tuple[float, float] = (0, 100),
        points: int = 1000,
        params: Optional[Dict[str, float]] = None,
        ic: Optional[Dict[str, float]] = None,
        rr: Optional[roadrunner.ExtendedRoadRunner] = None,
    ) -> SimulationResult:
        """
        Run deterministic ODE simulation.
        
        Parameters
        ----------
        model : str
            Model name ("PC", "PD_1", etc.)
        t_span : tuple of float
            (start_time, end_time)
        points : int
            Number of time points
        params : dict, optional
            Parameter overrides
        ic : dict, optional
            Initial condition overrides
            
        Returns
        -------
        dict
            {
                'time': np.ndarray of shape (points,),
                'states': dict mapping state_name -> np.ndarray of shape (points,)
            }
        """
        # Create simulator
        rr = self.prepare_runner(model, params=params, ic=ic, rr=rr)
        rr.integrator.setValue('maximum_num_steps', 100000)
        rr.integrator.setValue('relative_tolerance', 1e-8)
        rr.integrator.setValue('absolute_tolerance', 1e-10)
        rr.integrator.setValue('stiff', True)
        
        # Run simulation
        result = rr.simulate(t_span[0], t_span[1], points)
        
        # Extract results
        time = result[:, 0]
        state_names = self.model_constructor.get_state_names(model)
        
        states: dict[str, NDArray[np.float64]] = {}
        for i, name in enumerate(state_names):
            # Column 0 is time, states start at column 1
            states[name] = result[:, i + 1]
        
        return {
            'time': time,
            'states': states,
        }
    
    def run_stochastic(
        self,
        model: str,
        t_span: Tuple[float, float] = (0, 100),
        points: int = 1000,
        n_trajectories: int = 50,
        params: Optional[Dict[str, float]] = None,
        ic: Optional[Dict[str, float]] = None,
        seed: Optional[int] = None,
        rr: Optional[roadrunner.ExtendedRoadRunner] = None,
    ) -> List[SimulationResult]:
        """
        Run stochastic Gillespie simulations.
        
        Parameters
        ----------
        model : str
            Model name
        t_span : tuple of float
            (start_time, end_time)
        points : int
            Number of output time points
        n_trajectories : int
            Number of stochastic realizations
        params : dict, optional
            Parameter overrides
        ic : dict, optional
            Initial condition overrides
        seed : int, optional
            Random seed for reproducibility
            
        Returns
        -------
        list of dict
            List of n_trajectories results, each with format:
            {
                'time': np.ndarray,
                'states': dict[state_name -> np.ndarray]
            }
        """
        results = []

        for i in range(n_trajectories):

            # Create fresh simulator for each trajectory
            rr = self.prepare_runner(model, params=params, ic=ic, rr=rr)
            
            # Set seed if provided
            traj_seed = None if seed is None else seed + i
            
            # Run Gillespie
            result = rr.gillespie(t_span[0], t_span[1], points, seed=traj_seed)
            
            # Extract results
            time = result[:, 0]
            state_names = self.model_constructor.get_state_names(model)
            
            states = {}
            for j, name in enumerate(state_names):
                states[name] = result[:, j + 1]
            
            results.append({
                'time': time,
                'states': states,
            })
        
        return results
    
    def run_with_perturbations(
        self,
        model: str,
        t_span: Tuple[float, float] = (0, 100),
        points: int = 1000,
        perturbations: Optional[List[Dict]] = None,
        params: Optional[Dict[str, float]] = None,
        ic: Optional[Dict[str, float]] = None,
        method: str = 'deterministic',
        rr: Optional[roadrunner.ExtendedRoadRunner] = None,
    ) -> PerturbedSimulationResult:
        """
        Run simulation with segmented perturbations.
        
        Perturbations are applied by stopping the simulation, modifying the state
        or parameters, and continuing from the perturbed state.
        
        Parameters
        ----------
        model : str
            Model name
        t_span : tuple of float
            (start_time, end_time)
        points : int
            Total number of output time points
        perturbations : list of dict, optional
            Each perturbation is a dict with:
            - 'time': float - when to apply
            - 'type': str - 'step' (change species), 'knockout' (scale species), 
                           or 'parameter' (change parameter)
            - For 'step': 'species': str, 'value': float
            - For 'knockout': 'species': str, 'factor': float (e.g., 0.5 for 50% reduction)
            - For 'parameter': 'param': str, 'value': float
            
            Example: [
                {'time': 20, 'type': 'step', 'species': 'r', 'value': 2.0},
                {'time': 50, 'type': 'knockout', 'species': 'u_1', 'factor': 0.5},
                {'time': 80, 'type': 'parameter', 'param': 'gamma', 'value': 0.5}
            ]
        params : dict, optional
            Parameter overrides
        ic : dict, optional
            Initial condition overrides
        method : str
            'deterministic' or 'stochastic'
            
        Returns
        -------
        dict
            {
                'time': np.ndarray,
                'states': dict[state_name -> np.ndarray],
                'perturbation_times': list of float
            }
        """
        if perturbations is None:
            perturbations = []
        
        # Sort perturbations by time
        perturbations = sorted(perturbations, key=lambda p: p['time'])
        
        # Validate perturbation times
        for p in perturbations:
            if not (t_span[0] < p['time'] < t_span[1]):
                raise ValueError(
                    f"Perturbation time {p['time']} outside simulation range {t_span}"
                )
        
        # Deduplicate segment boundaries (handles same-time perturbations)
        unique_times = [t for t, _ in groupby(p['time'] for p in perturbations)]
        segment_times = [t_span[0]] + unique_times + [t_span[1]]
        
        # Group perturbations by time for application
        perts_by_time: dict[float, list] = defaultdict(list)
        for p in perturbations:
            perts_by_time[p['time']].append(p)
        
        # Allocate arrays for results
        time_full = np.linspace(t_span[0], t_span[1], points)
        state_names = self.model_constructor.get_state_names(model)
        states_full: dict[str, NDArray[np.float64]] = {
            name: np.zeros(points, dtype=float) for name in state_names
        }
        
        # Create initial simulator
        rr = self.prepare_runner(model, params=params, ic=ic, rr=rr)
        rr.integrator.setValue('maximum_num_steps', 100000)
        rr.integrator.setValue('relative_tolerance', 1e-8)
        rr.integrator.setValue('absolute_tolerance', 1e-10)
        rr.integrator.setValue('stiff', True)
        rr.integrator.setValue('initial_time_step', 1e-3)
        rr.integrator.setValue('minimum_time_step', 1e-8)
        
        # Track current parameter values (for parameter perturbations)
        current_params = dict(params) if params is not None else {}
        
        # Simulate each segment
        for seg_idx in range(len(segment_times) - 1):
            t_start = segment_times[seg_idx]
            t_end = segment_times[seg_idx + 1]
            
            # Find time points in this segment
            if t_end == segment_times[-1]:
                seg_mask = (time_full >= t_start) & (time_full <= t_end)   # final segment
            else:
                seg_mask = (time_full >= t_start) & (time_full < t_end)    # avoid duplicate boundary
            
            if len(segment_times) == 0:
                continue
            
            # Run segment
            if method == 'deterministic':
                seg_result = rr.simulate(t_start, t_end, int(np.sum(seg_mask)))
            elif method == 'stochastic':
                seg_result = rr.gillespie(t_start, t_end, int(np.sum(seg_mask)))
            else:
                raise ValueError(f"Unknown method '{method}'")
            
            # Store results
            for i, name in enumerate(state_names):
                states_full[name][seg_mask] = seg_result[:, i + 1]
            
            # Apply perturbation at end of segment (if not last segment)
            if seg_idx < len(unique_times):
                t_pert = unique_times[seg_idx]
                for pert in perts_by_time[t_pert]:
                    if pert['type'] == 'step':
                        setattr(rr, pert['species'], pert['value'])
                    elif pert['type'] == 'knockout':
                        current_value = getattr(rr, pert['species'])
                        setattr(rr, pert['species'], current_value * pert['factor'])
                    elif pert['type'] == 'parameter':
                        setattr(rr, pert['param'], pert['value'])
                        current_params[pert['param']] = pert['value']
                    else:
                        raise ValueError(f"Unknown perturbation type '{pert['type']}'")
            # After applying perturbations, before next segment     
            rr.integrator.setValue('maximum_num_steps', 100000)
            rr.integrator.setValue('stiff', True)
            rr.integrator.setValue('initial_time_step', 1e-3)
            rr.integrator.setValue('minimum_time_step', 1e-8)
        
        # Reset integrator after perturbation 
        rr.integrator.setValue('maximum_num_steps', 100000)
        rr.integrator.setValue('stiff', True)
        
        return {
            'time': time_full,
            'states': states_full,
            'perturbation_times': [p['time'] for p in perturbations],
        }
