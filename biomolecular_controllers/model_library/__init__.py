"""model_library — modular registry of Antimony models, defaults, and factory.

Exports the canonical package-wide registries:
  MODELS                     dict[str, str]
  DEFAULT_PARAMS             dict[str, dict[str, float]]
  DEFAULT_INITIAL_CONDITIONS dict[str, dict[str, float]]

And the factory / metadata:
  Models          -- RoadRunner factory class
  GAIN_LABELS     -- x-axis label strings for metric plots
  STATE_VARIABLES -- state variable names per model

Family modules:
  hpt_models.py / hpt_defaults.py          -- HPT-axis models
  controller_models.py / controller_defaults.py -- P/PD/PI/PID controller models
  models.py                                -- Models factory class and metadata
"""

from .registry import MODELS, DEFAULT_PARAMS, DEFAULT_INITIAL_CONDITIONS
from .models import Models, GAIN_LABELS, STATE_VARIABLES

__all__ = [
    "MODELS",
    "DEFAULT_PARAMS",
    "DEFAULT_INITIAL_CONDITIONS",
    "Models",
    "GAIN_LABELS",
    "STATE_VARIABLES",
]
