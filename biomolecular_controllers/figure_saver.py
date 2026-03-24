"""
figure_saver.py — save Bokeh figures with consistent filenames.

Naming convention:  {output_dir}/{controller}_{metric}.{ext}
  e.g.  figures/PI_1_steady_state_error.html
        figures/PC_stability_diagram.png

Usage (drop this at the bottom of any notebook cell after show()):
    from figure_saver import save_fig

    save_fig(layout, 'PI_1', 'steady_state_error')           # HTML (default)
    save_fig(layout, 'PC',   'root_locus', fmt='png')         # PNG  (needs driver)
    save_fig(layout, 'PD_1', 'overshoot',  output_dir='../my_figs')  # custom dir
"""

import os
import re
from pathlib import Path
from typing import Union
from bokeh.io import show as bokeh_show
from bokeh.resources import CDN
from bokeh.io import save as bokeh_save, output_notebook
from bokeh.io import export_png

PathLike = Union[str, Path]

DEFAULT_DIR = Path(
    os.getenv(
        "BIOCONTROLLERS_OUTPUT_DIR",
        Path(__file__).resolve().parents[1] / "figures",
    )
)

# ── helpers ──────────────────────────────────────────────────────────────────
def resolve_output_dir(output_dir: PathLike | None = None) -> Path:
    out = Path(output_dir) if output_dir is not None else DEFAULT_DIR
    out.mkdir(parents=True, exist_ok=True)
    return out

def _slugify(text: str) -> str:
    """'Steady-State Error (α₁)' → 'steady_state_error_1'"""
    text = text.lower()
    # replace non-ascii with nothing, spaces/hyphens/dots with underscore
    text = text.encode('ascii', 'ignore').decode()
    text = re.sub(r'[\s\-\.]+', '_', text)
    text = re.sub(r'[^a-z0-9_]', '', text)
    return text.strip('_')


def _make_path(output_dir: PathLike | None, controller: str, metric: str, ext: str) -> Path:
    out = resolve_output_dir(output_dir)
    fname = f"{_slugify(controller)}_{_slugify(metric)}.{ext}"
    return out / fname



# ── main API ─────────────────────────────────────────────────────────────────

def save_fig(
    fig,
    controller: str,
    metric: str,
    output_dir: PathLike | None = None,
    fmt: str = 'html',
    show: bool = False,
) -> Path:
    """
    Save a Bokeh figure/layout to disk.

    Parameters
    ----------
    fig         : Bokeh figure or layout returned by plotter.*
    controller  : e.g. 'PI_1', 'PC', 'PD_1', 'PID_2'
    metric      : e.g. 'steady_state_error', 'root_locus', 'stability_diagram'
    output_dir  : directory to write into (created if absent), default 'figures'
    fmt         : 'html' (default, always works) or 'png' (needs selenium/playwright)
    show        : if True, display the figure inline in the notebook after saving

    Returns
    -------
    Path to the saved file.
    """
    
    fmt = fmt.lower()

    if fmt == 'html':
        path = _save_html(fig, controller, metric, output_dir)
    elif fmt == 'png':
        path = _save_png(fig, controller, metric, output_dir)
    else:
        raise ValueError(f"fmt must be 'html' or 'png', got '{fmt}'")

    if show:
        bokeh_show(fig)

    return path


def save_all(
    figures: dict,
    output_dir: str = 'figures',
    fmt: str = 'html',
    show: bool = False,
) -> list:
    """
    Save multiple figures at once.

    Parameters
    ----------
    figures : dict mapping (controller, metric) tuples to Bokeh figures, e.g.
                {
                    ('PI_1', 'steady_state_error'): layout1,
                    ('PI_1', 'overshoot'):           layout2,
                    ('PC',   'root_locus'):           layout3,
                }
    output_dir, fmt : same as save_fig

    Returns
    -------
    List of saved Path objects.
    """
    saved = []
    for (controller, metric), fig in figures.items():
        path = save_fig(fig, controller, metric, output_dir=output_dir, fmt=fmt, show=show)
        saved.append(path)
    return saved


# ── format-specific helpers ───────────────────────────────────────────────────

def _save_html(fig, controller, metric, output_dir):

    path = _make_path(output_dir, controller, metric, 'html')
    bokeh_save(
        fig,
        filename=str(path),
        resources=CDN,
        title=f'{controller} — {metric}',
    )
    output_notebook(hide_banner=True)  # restore inline display after bokeh_save resets it
    print(f'  Saved HTML → {path}')
    return path


def _save_png(fig, controller, metric, output_dir):

    path = _make_path(output_dir, controller, metric, 'png')
    try:
        export_png(fig, filename=str(path))
        print(f'  Saved PNG  → {path}')
        return path
    except Exception as e:
        msg = str(e)
        # give a helpful hint for the most common failure (missing driver)
        if 'selenium' in msg.lower() or 'webdriver' in msg.lower() or 'playwright' in msg.lower() or 'driver' in msg.lower():
            print(
                f'\n  PNG export failed — no headless browser driver found.\n'
                f'  Install one of:\n'
                f'      pip install selenium     # then install chromedriver or geckodriver\n'
                f'      pip install playwright   # then: playwright install chromium\n'
                f'  Falling back to HTML.\n'
            )
        else:
            print(f'\n  PNG export failed ({msg})\n  Falling back to HTML.\n')
        return _save_html(fig, controller, metric, output_dir)