"""
Tellurium factory and supporting metadata for biomolecular controller models.

Contains:
  Models       -- factory class for creating configured RoadRunner instances
  GAIN_LABELS  -- x-axis label strings for metric plots
  STATE_VARIABLES -- state variable names per model (for validation)
"""

from typing import Dict, Optional, cast
from tellurium import roadrunner
import tellurium as te

from .registry import MODELS, DEFAULT_PARAMS, DEFAULT_INITIAL_CONDITIONS, STATE_VARIABLES

# Internal alias kept for clarity within this module
ANTIMONY_MODELS = MODELS

# X-axis labels for metric plots, showing the gain/SS-error formula per model.
# Update non-PC entries once the SS-error expressions are confirmed from the paper.
GAIN_LABELS = {
    "PC":    "Gain (G = α₁·θ₁·θ₂)",
    "PD_1":  "Gain (PD_1 — update with SS-error formula)",
    "PD_2":  "Gain (PD_2 — update with SS-error formula)",
    "PI_1":  "Gain (PI_1 — update with SS-error formula)",
    "PI_2":  "Gain (PI_2 — update with SS-error formula)",
    "PID_1": "Gain (PID_1 — update with SS-error formula)",
    "PID_2": "Gain (PID_2 — update with SS-error formula)",
    "HPT_dimensionless":  "HPT axis dimensionless",
    "HPT_full": "HPT axis full model",
}


class Models:
    """
    Factory for creating configured Tellurium RoadRunner instances.

    Examples
    --------
    >>> model = Models()
    >>> rr = model.create_roadrunner("PC")
    >>> result = rr.simulate(0, 100, 1000)

    >>> # Custom parameters
    >>> rr = model.create_roadrunner("PD_1", params={"gamma": 0.5})

    >>> # Override reference signal
    >>> rr = model.create_roadrunner("PI_1", params={"ref": 2.0})
    """

    def __init__(self):
        self.available_models = list(ANTIMONY_MODELS.keys())

    def create_roadrunner(
        self,
        model: str,
        params: Optional[Dict[str, float]] = None,
        ic: Optional[Dict[str, float]] = None,
    ) -> roadrunner.ExtendedRoadRunner:
        """
        Create a configured RoadRunner instance.

        Parameters
        ----------
        model : str
            Model name ("PC", "PD_1", "PD_2", "PI_1", "PI_2", "PID_1", "PID_2")
        params : dict, optional
            Parameter overrides. If not provided, uses DEFAULT_PARAMS[model]
        ic : dict, optional
            Initial condition overrides. If not provided, uses DEFAULT_INITIAL_CONDITIONS[model]

        Returns
        -------
        roadrunner.RoadRunner
            Configured simulator instance

        Raises
        ------
        ValueError
            If model name is invalid or unknown parameters/states provided
        """
        if model not in ANTIMONY_MODELS:
            raise ValueError(
                f"Unknown model '{model}'. Available: {self.available_models}"
            )

        # Merge parameters with defaults
        full_params = dict(DEFAULT_PARAMS[model])
        if params is not None:
            unknown_params = set(params.keys()) - set(full_params.keys())
            if unknown_params:
                raise ValueError(
                    f"Unknown parameters for model '{model}': {unknown_params}"
                )
            full_params.update(params)

        # Merge initial conditions with defaults
        full_ic = dict(DEFAULT_INITIAL_CONDITIONS[model])
        if ic is not None:
            unknown_states = set(ic.keys()) - set(full_ic.keys())
            if unknown_states:
                raise ValueError(
                    f"Unknown states for model '{model}': {unknown_states}"
                )
            full_ic.update(ic)

        # Format Antimony string and create RoadRunner instance
        antimony_str = ANTIMONY_MODELS[model].format(**full_params, **full_ic)
        rr = cast(roadrunner.ExtendedRoadRunner, te.loada(antimony_str))

        return rr

    def get_state_names(self, model: str) -> list:
        """Get list of state variable names for a model."""
        if model not in STATE_VARIABLES:
            raise ValueError(f"Unknown model '{model}'")
        return STATE_VARIABLES[model]

    def get_default_params(self, model: str) -> Dict[str, float]:
        """Get default parameters for a model."""
        if model not in DEFAULT_PARAMS:
            raise ValueError(f"Unknown model '{model}'")
        return dict(DEFAULT_PARAMS[model])

    def get_default_ic(self, model: str) -> Dict[str, float]:
        """Get default initial conditions for a model."""
        if model not in DEFAULT_INITIAL_CONDITIONS:
            raise ValueError(f"Unknown model '{model}'")
        return dict(DEFAULT_INITIAL_CONDITIONS[model])
