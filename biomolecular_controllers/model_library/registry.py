"""Central registry for all Antimony models and defaults.

Merges all family-level dictionaries into the canonical package-wide registries:
  MODELS                  -- dict[str, str]          Antimony model strings
  DEFAULT_PARAMS          -- dict[str, dict[str, float]]
  DEFAULT_INITIAL_CONDITIONS -- dict[str, dict[str, float]]
  STATE_VARIABLES         -- dict[str, list[str]]

Also validates at import time that DEFAULT_PARAMS and DEFAULT_INITIAL_CONDITIONS
contain only keys that exist in MODELS, and that no model name appears in more
than one family module.
"""

from .hpt_models import MODELS_HPT
from .hpt_defaults import DEFAULT_PARAMS_HPT, DEFAULT_INITIAL_CONDITIONS_HPT, STATE_VARIABLES_HPT

from .controller_models import MODELS_CONTROLLER
from .controller_defaults import (
    DEFAULT_PARAMS_CONTROLLER,
    DEFAULT_INITIAL_CONDITIONS_CONTROLLER,
    STATE_VARIABLES_CONTROLLER,
)


def _merge_unique(*named_dicts):
    """Merge dictionaries, raising ValueError on duplicate keys.

    Parameters
    ----------
    *named_dicts : tuple[str, dict]
        Each element is a (name, dict) pair used for error reporting.

    Returns
    -------
    dict
        Merged dictionary with all key-value pairs.

    Raises
    ------
    ValueError
        If any key appears in more than one source dictionary.
    """
    merged = {}
    seen = set()
    for name, d in named_dicts:
        overlap = seen.intersection(d.keys())
        if overlap:
            raise ValueError(
                f"Duplicate model keys while merging {name}: {sorted(overlap)}"
            )
        merged.update(d)
        seen.update(d.keys())
    return merged


MODELS = _merge_unique(
    ("MODELS_HPT", MODELS_HPT),
    ("MODELS_CONTROLLER", MODELS_CONTROLLER),
)

DEFAULT_PARAMS = _merge_unique(
    ("DEFAULT_PARAMS_HPT", DEFAULT_PARAMS_HPT),
    ("DEFAULT_PARAMS_CONTROLLER", DEFAULT_PARAMS_CONTROLLER),
)

DEFAULT_INITIAL_CONDITIONS = _merge_unique(
    ("DEFAULT_INITIAL_CONDITIONS_HPT", DEFAULT_INITIAL_CONDITIONS_HPT),
    ("DEFAULT_INITIAL_CONDITIONS_CONTROLLER", DEFAULT_INITIAL_CONDITIONS_CONTROLLER),
)

STATE_VARIABLES = _merge_unique(
    ("STATE_VARIABLES_HPT", STATE_VARIABLES_HPT),
    ("STATE_VARIABLES_CONTROLLER", STATE_VARIABLES_CONTROLLER),
)

# Validate consistency: every key in the defaults dicts must exist in MODELS
_missing_params = sorted(set(DEFAULT_PARAMS) - set(MODELS))
_missing_ics = sorted(set(DEFAULT_INITIAL_CONDITIONS) - set(MODELS))

if _missing_params:
    raise ValueError(
        f"DEFAULT_PARAMS contains keys missing from MODELS: {_missing_params}"
    )
if _missing_ics:
    raise ValueError(
        f"DEFAULT_INITIAL_CONDITIONS contains keys missing from MODELS: {_missing_ics}"
    )

_missing_sv = sorted(set(STATE_VARIABLES) - set(MODELS))
if _missing_sv:
    raise ValueError(
        f"STATE_VARIABLES contains keys missing from MODELS: {_missing_sv}"
    )
