"""Compatibility shim — all definitions have moved to model_library.

Import from biomolecular_controllers.model_library instead.
"""

from .model_library import (
    MODELS,
    DEFAULT_PARAMS,
    DEFAULT_INITIAL_CONDITIONS,
    Models,
    GAIN_LABELS,
    STATE_VARIABLES,
)

# Legacy alias
ANTIMONY_MODELS = MODELS

__all__ = [
    "MODELS",
    "ANTIMONY_MODELS",
    "DEFAULT_PARAMS",
    "DEFAULT_INITIAL_CONDITIONS",
    "Models",
    "GAIN_LABELS",
    "STATE_VARIABLES",
]
