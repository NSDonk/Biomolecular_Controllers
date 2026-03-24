"""
figure_gifs.py — export Bokeh sweep layouts as animated GIFs.

Designed to work with VisualizationPipeline.plot_metric_interactive() layouts.
Pass traj_renderers / ref_renderers / metric_source in directly alongside the layout. 


Typical usage
-------------
    from biomolecular_controllers.figure_gifs import save_gif

    save_gif(
        layout           = layout,
        traj_renderers   = traj_renderers,
        ref_renderers    = ref_renderers,
        metric_source    = metric_source,
        controller       = 'PC',
        metric           = 'steady_state_error',
    )
"""

from __future__ import annotations

import re
import os
import shutil
from pathlib import Path
from typing import List, Sequence, Union

import imageio.v3 as iio
from bokeh.io import export_png
from bokeh.models import ColumnDataSource, GlyphRenderer


# ── default output directory ───────────────────────
PathLike = Union[str, Path]
DEFAULT_DIR = Path(os.getenv("BIOCONTROLLERS_OUTPUT_DIR", Path(__file__).resolve().parents[1] / "figures"))

# ── helpers ───────────────────────────────────────────────────────────────────
def resolve_output_dir(output_dir: PathLike | None = None) -> Path:
    out = Path(output_dir) if output_dir is not None else DEFAULT_DIR
    out.mkdir(parents=True, exist_ok=True)
    return out

def _slugify(text: str) -> str:
    text = text.lower().encode('ascii', 'ignore').decode()
    text = re.sub(r'[\s\-\.]+', '_', text)
    text = re.sub(r'[^a-z0-9_]', '', text)
    return text.strip('_')


def _make_gif_path(output_dir: PathLike | None, controller: str, metric: str) -> Path:
    out = resolve_output_dir(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    return out / f"{_slugify(controller)}_{_slugify(metric)}.gif"


# ── main API ──────────────────────────────────────────────────────────────────

def save_gif(
    layout,
    traj_renderers: Sequence[GlyphRenderer],
    ref_renderers: Sequence[GlyphRenderer],
    metric_source: ColumnDataSource,
    controller: str,
    metric: str,
    output_dir: PathLike | None = None,
    fps: int = 4,
    cleanup_frames: bool = True,
) -> Path:
    """
    Export a sweep layout as an animated GIF.

    Each frame shows one trajectory (hiding all others) and highlights the
    corresponding point in the metric panel.

    Parameters
    ----------
    layout
        Bokeh layout returned by plot_metric_interactive().
    traj_renderers
        List of trajectory line renderers (one per sweep point), in the same
        order as metric_source rows.
    ref_renderers
        List of reference line renderers (one per sweep point).
    metric_source
        ColumnDataSource from the metric panel (has 'x', 'y', 'color' columns).
    controller
        e.g. 'PC', 'PI_1' — used in the output filename.
    metric
        e.g. 'steady_state_error' — used in the output filename.
    output_dir
        Directory to write the GIF (and temporary frames) into.
    fps
        Frames per second.
    cleanup_frames
        Delete temporary PNG frames after stitching.

    Returns
    -------
    Path to the saved GIF.
    """
    n = len(traj_renderers)
    if n == 0:
        raise ValueError("traj_renderers is empty.")
    if len(ref_renderers) != n:
        raise ValueError("traj_renderers and ref_renderers must have the same length.")

    output_dir = resolve_output_dir(output_dir)
    frames_dir = Path(output_dir) / f"_gif_frames_{_slugify(controller)}_{_slugify(metric)}"
    frames_dir.mkdir(parents=True, exist_ok=True)

    # Preserve original visibility so we can restore after export
    orig_traj_vis = [r.visible for r in traj_renderers]
    orig_ref_vis  = [r.visible for r in ref_renderers]

    frame_files: List[Path] = []

    try:
        for i in range(n):
            # Show only the i-th trajectory
            for j in range(n):
                traj_renderers[j].visible = (j == i)
                ref_renderers[j].visible  = (j == i)

            frame_path = frames_dir / f"frame_{i:04d}.png"
            export_png(layout, filename=str(frame_path))
            frame_files.append(frame_path)

    finally:
        # Restore original visibility even if export fails mid-way
        for j in range(n):
            traj_renderers[j].visible = orig_traj_vis[j]
            ref_renderers[j].visible  = orig_ref_vis[j]

    # Stitch frames into GIF
    gif_path = _make_gif_path(output_dir, controller, metric)
    images = [iio.imread(str(f)) for f in frame_files]
    iio.imwrite(str(gif_path), images, extension='.gif', duration=int(1000/fps), loop=0)
    print(f'  Saved GIF  → {gif_path}')

    if cleanup_frames:
        shutil.rmtree(frames_dir, ignore_errors=True)

    return gif_path


def save_all_gifs(
    gif_specs: list,
    output_dir: PathLike | None = None,
    fps: int = 4,
) -> list:
    """
    Save multiple GIFs at once.

    Parameters
    ----------
    gif_specs : list of dicts, each with keys:
        layout, traj_renderers, ref_renderers, metric_source, controller, metric

    Returns
    -------
    List of saved Path objects.
    """
    return [
        save_gif(output_dir=output_dir, fps=fps, **spec)
        for spec in gif_specs
    ]