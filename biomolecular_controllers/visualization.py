"""
Visualization pipeline for biomolecular controller analysis.

Uses Bokeh for interactive plots with optional SVG export.
"""

from typing import Dict, List, Optional, Tuple, Any, cast, Literal, Union
import numpy as np
import bokeh.plotting
import bokeh.layouts
import bokeh.models
import bokeh.palettes
from bokeh.models import (
    FactorRange, Range1d, DataRange1d, Span, Label, Arrow, VeeHead,
    HoverTool, CrosshairTool, LinearColorMapper, ColorBar,
    ColumnDataSource, CustomJS, Button, TapTool, Div, LayoutDOM,
    Legend, LegendItem, GlyphRenderer, Band, Whisker
)
from bokeh.transform import factor_cmap
from bokeh.palettes import Viridis256
from bokeh.io import export_svgs


class VisualizationPipeline:
    """
    Create publication-ready Bokeh figures for controller analysis.
    
    Examples
    --------
    >>> viz = VisualizationPipeline()
    >>> p = viz.plot_time_series(time, states, title="PC Controller Response")
    >>> bokeh.plotting.show(p)
    
    >>> # With perturbation and metrics
    >>> p = viz.plot_overshoot(time, y, overshoot_metric, pert_time=20)
    >>> bokeh.plotting.save(p, "overshoot.html")
    """
    
    def __init__(
        self,
        width: int = 800,
        height: int = 400,
        palette: str = "Category10",
    ):
        """
        Initialize visualization pipeline.
        
        Parameters
        ----------
        width : int
            Default plot width in pixels
        height : int
            Default plot height in pixels
        palette : str
            Bokeh palette name for colors
        """
        self.width = width
        self.height = height
        self.palette = getattr(bokeh.palettes, palette)
        
    def plot_time_series(
        self,
        time: np.ndarray,
        states: Dict[str, np.ndarray],
        title: str = "Time Series",
        highlight_regions: Optional[List[Tuple[float, float, str]]] = None,
        perturbation_times: Optional[List[float]] = None,
        width: Optional[int] = None,
        height: Optional[int] = None,
    ) -> bokeh.plotting.figure:
        """
        Plot multiple state variables over time.
        
        Parameters
        ----------
        time : np.ndarray
            Time points
        states : dict
            Dictionary mapping state name to array
        title : str
            Plot title
        highlight_regions : list of (start, end, label), optional
            Regions to highlight with shading
        perturbation_times : list of float, optional
            Times to mark with vertical lines
        width, height : int, optional
            Override default dimensions
            
        Returns
        -------
        bokeh.plotting.figure
        """
        w = width or self.width
        h = height or self.height
        
        p = bokeh.plotting.figure(
            width=w,
            height=h,
            title=title,
            x_axis_label="Time",
            y_axis_label="Concentration",
            tools="pan,wheel_zoom,box_zoom,reset,save",
        )
        
        # Add hover tool
        hover = HoverTool(tooltips=[("Time", "@x{0.00}"), ("Value", "@y{0.000}")])
        p.add_tools(hover)
        p.add_tools(CrosshairTool())
        
        # Plot each state
        colors = self.palette[max(10, len(states))]
        for i, (name, values) in enumerate(states.items()):
            color = colors[i % len(colors)]
            p.line(
                time,
                values,
                legend_label=name,
                color=color,
                line_width=2,
                alpha=0.8,
            )
        
        # Add perturbation markers
        if perturbation_times:
            for t in perturbation_times:
                vline = Span(
                    location=t,
                    dimension='height',
                    line_color='red',
                    line_dash='dashed',
                    line_width=2,
                    line_alpha=0.5,
                )
                p.add_layout(vline)

                label = Label(
                    x=t,
                    y=0.0,
                    text=f" Pert @ t={t}",
                    text_font_size="10pt",
                    text_color="red",
                    text_baseline="bottom",
                    angle=90,
                    angle_units="deg",
                )
                p.add_layout(label)
        
        # Add highlighted regions
        if highlight_regions:
            from bokeh.models import BoxAnnotation
            for start, end, label_text in highlight_regions:
                box = BoxAnnotation(
                    left=start,
                    right=end,
                    fill_alpha=0.1,
                    fill_color='gray',
                )
                p.add_layout(box)
        
        p.legend.click_policy = "hide"
        p.legend.location = "top_right"
        
        return p
    
    def plot_overshoot(
        self,
        time: np.ndarray,
        y: np.ndarray,
        overshoot_metric: Dict[str, float],
        ref: float,
        pert_start: float,
        title: str = "Overshoot Analysis",
    ) -> bokeh.plotting.figure:
        """
        Plot time series with overshoot annotation.
        
        Parameters
        ----------
        time : np.ndarray
            Time points
        y : np.ndarray
            Output signal
        overshoot_metric : dict
            Output from MetricsCalculator.overshoot()
        ref : float
            Reference value
        pert_start : float
            Perturbation start time
        title : str
            Plot title
            
        Returns
        -------
        bokeh.plotting.figure
        """
        p = bokeh.plotting.figure(
            width=self.width,
            height=self.height,
            title=title,
            x_axis_label="Time",
            y_axis_label="Output",
            tools="pan,wheel_zoom,box_zoom,reset,save",
        )
        
        # Plot output
        p.line(time, y, legend_label="Output y", color="navy", line_width=2)
        
        # Plot reference
        p.line(time, np.full_like(time, ref), 
               legend_label=f"Reference (r={ref})",
               color="red", line_dash="dashed", line_width=2)
        
        # Mark perturbation start
        vline = Span(location=pert_start, dimension='height',
                     line_color='orange', line_dash='dotted', line_width=2)
        p.add_layout(vline)
        
        # Mark peak (overshoot point)
        peak_time = overshoot_metric['peak_time']
        peak_value = overshoot_metric['peak_value']
        
        p_any = cast(Any, p)
        p_any.scatter(
            [peak_time],
            [peak_value],
            size=10,
            color="red",
            marker="circle",
            legend_label=f"Peak (OS={overshoot_metric['magnitude']:.1%})",
        )
        
        # Arrow from ref to peak
        arrow = Arrow(
            end=VeeHead(size=10, fill_color="red"),
            x_start=peak_time,
            y_start=ref,
            x_end=peak_time,
            y_end=peak_value,
        )

        arrow_any = cast(Any, arrow)
        arrow_any.line_color = "red"
        arrow_any.line_width = 2
        
        # Add text annotation
        overshoot_pct = overshoot_metric['magnitude'] * 100
        time_to_peak = overshoot_metric['time_to_peak']
        
        label = Label(
            x=peak_time + (time[-1] - time[0]) * 0.02,
            y=(ref + peak_value) / 2,
            text=f"Overshoot: {overshoot_pct:.1f}%\nTime: {time_to_peak:.1f}s",
            text_font_size="11pt",
            background_fill_color="white",
            background_fill_alpha=0.8,
        )
        p.add_layout(label)
        
        p.legend.location = "bottom_right"
        
        return p
    
    def plot_settling_time(
        self,
        time: np.ndarray,
        y: np.ndarray,
        settling_metric: Dict[str, float],
        ref: float,
        tolerance: float = 0.02,
        pert_start: Optional[float] = None,
        title: str = "Settling Time Analysis",
    ) -> bokeh.plotting.figure:
        """
        Plot time series with settling time annotation.
        
        Parameters
        ----------
        time : np.ndarray
            Time points
        y : np.ndarray
            Output signal
        settling_metric : dict
            Output from MetricsCalculator.settling_time()
        ref : float
            Reference value
        tolerance : float
            Tolerance band (fractional)
        pert_start : float, optional
            Perturbation start time
        title : str
            Plot title
            
        Returns
        -------
        bokeh.plotting.figure
        """
        p = bokeh.plotting.figure(
            width=self.width,
            height=self.height,
            title=title,
            x_axis_label="Time",
            y_axis_label="Output",
            tools="pan,wheel_zoom,box_zoom,reset,save",
        )
        
        # Plot output
        p.line(time, y, legend_label="Output y", color="navy", line_width=2)
        
        # Plot reference and tolerance band
        threshold = tolerance * abs(ref)
        upper = ref + threshold
        lower = ref - threshold
        
        p.line(time, np.full_like(time, ref), 
               legend_label=f"Reference (r={ref})",
               color="red", line_dash="dashed", line_width=2)
        
        # Tolerance band
        from bokeh.models import Band
        p.line(time, np.full_like(time, upper), 
               color="green", line_dash="dotted", line_width=1, alpha=0.5)
        p.line(time, np.full_like(time, lower), 
               color="green", line_dash="dotted", line_width=1, alpha=0.5)
        
        # Mark settling time
        if settling_metric['settled']:
            ts_abs = settling_metric['settling_time_abs']
            vline = Span(location=ts_abs, dimension='height',
                         line_color='green', line_width=2)
            p.add_layout(vline)
            
            ts_rel = settling_metric['settling_time']
            label = Label(
                x=ts_abs,
                y=ref,
                text=f" Settling time: {ts_rel:.1f}s\n (±{tolerance*100:.0f}% band)",
                text_font_size="11pt",
                background_fill_color="white",
                background_fill_alpha=0.8,
            )
            p.add_layout(label)
        else:
            # Did not settle
            label = Label(
                x=time[len(time)//2],
                y=ref,
                text="Did not settle within time window",
                text_font_size="12pt",
                text_color="red",
                background_fill_color="white",
                background_fill_alpha=0.8,
            )
            p.add_layout(label)
        
        # Mark perturbation start if provided
        if pert_start is not None:
            vline = Span(location=pert_start, dimension='height',
                         line_color='orange', line_dash='dotted', line_width=2)
            p.add_layout(vline)
        
        p.legend.location = "bottom_right"
        
        return p
    
    def plot_model_comparison(
        self,
        results: Dict[str, Dict],
        metric: str = 'overshoot',
        title: Optional[str] = None,
    ) -> bokeh.plotting.figure:
        """
        Compare a metric across multiple models.
        
        Parameters
        ----------
        results : dict
            Dictionary mapping model_name -> metrics_dict
            Example: {"PC": {"overshoot": 0.2, "settling_time": 45}, ...}
        metric : str
            Which metric to plot (key in metrics_dict)
        title : str, optional
            Plot title
            
        Returns
        -------
        bokeh.plotting.figure
        """
        models = list(results.keys())
        values = [results[m][metric] for m in models]
        
        if title is None:
            title = f"{metric.replace('_', ' ').title()} Comparison"
        
        p = bokeh.plotting.figure(
            x_range=FactorRange(*models),
            width=self.width,
            height=self.height,
            title=title,
            x_axis_label="Model",
            y_axis_label=metric.replace('_', ' ').title(),
            tools="save",
        )
        
        # Bar chart
        colors = self.palette[max(10, len(models))]
        p.vbar(
            x=models,
            top=values,
            width=0.7,
            color=colors[:len(models)],
            alpha=0.8,
        )
        
        # Add value labels on bars
        for i, (model, value) in enumerate(zip(models, values)):
            label = Label(
                x=i,
                y=value,
                text=f"{value:.2f}",
                text_align="center",
                text_baseline="bottom",
                text_font_size="10pt",
            )
            p.add_layout(label)
        
        p.xgrid.grid_line_color = None
        p.y_range = DataRange1d(range_padding=0.1)
        
        return p
    
    @staticmethod
    def save_svg(plot: bokeh.plotting.figure, filename: str):
        """
        Export plot to SVG (requires selenium and geckodriver/chromedriver).
        
        Parameters
        ----------
        plot : bokeh.plotting.figure
            Plot to export
        filename : str
            Output filename (should end in .svg)
        """
        plot.output_backend = "svg"
        export_svgs(plot, filename=filename)
    
    @staticmethod
    def create_dashboard(
        plots: List[bokeh.plotting.figure],
        ncols: int = 2,
    ) -> LayoutDOM:
        """
        Arrange multiple plots in a grid layout.
        
        Parameters
        ----------
        plots : list of figure
            Plots to arrange
        ncols : int
            Number of columns in grid
            
        Returns
        -------
        bokeh.layouts.LayoutDOM
            Grid layout that can be shown or saved
        """
        rows = []
        for i in range(0, len(plots), ncols):
            row_plots = plots[i:i+ncols]
            rows.append(bokeh.layouts.row(*row_plots))
        
        return bokeh.layouts.column(*rows)
    
    def plot_bifurcation_diagram(
        self,
        param_values: np.ndarray,
        steady_states: np.ndarray,
        stable: np.ndarray,
        param_name: str,
        state_name: str = "y",
        title: Optional[str] = None,
    ) -> bokeh.plotting.figure:
        """
        Plot bifurcation diagram (steady state vs parameter).
        
        Parameters
        ----------
        param_values : np.ndarray
            Parameter values
        steady_states : np.ndarray
            Steady-state values of tracked variable
        stable : np.ndarray
            Boolean stability array
        param_name : str
            Parameter being varied
        state_name : str
            Name of state variable
        title : str, optional
            Plot title
            
        Returns
        -------
        bokeh.plotting.figure
        """
        if title is None:
            title = f"Bifurcation Diagram: {state_name} vs {param_name}"
        
        p = bokeh.plotting.figure(
            width=self.width,
            height=self.height,
            title=title,
            x_axis_label=param_name,
            y_axis_label=f"Steady State {state_name}",
            tools="pan,wheel_zoom,box_zoom,reset,save",
            x_axis_type="log",
        )
        
        # Separate stable and unstable branches
        stable_params = param_values[stable]
        stable_ss = steady_states[stable]
        unstable_params = param_values[~stable]
        unstable_ss = steady_states[~stable]
        
        # Plot stable branch (solid)
        if len(stable_params) > 0:
            p.line(stable_params, stable_ss, line_width=2, color="blue",
                   legend_label="Stable")
            p.scatter(stable_params, stable_ss, size=4, color="blue", marker="circle")
        
        # Plot unstable branch (dashed)
        if len(unstable_params) > 0:
            p.line(unstable_params, unstable_ss, line_width=2, color="red",
                   line_dash="dashed", legend_label="Unstable")
            p.scatter(unstable_params, unstable_ss, size=4, color="red", marker="circle")
        
        p.legend.location = "top_right"
        
        return p
    
    def plot_sobol_indices(
        self,
        sobol_results: Dict[str, Dict[str, float]],
        title: str = "Sobol Sensitivity Indices",
    ) -> bokeh.plotting.figure:
        """
        Plot Sobol sensitivity indices as bar chart.
        
        Parameters
        ----------
        sobol_results : dict
            Output from SensitivityAnalyzer.sobol_analysis()
            Must contain 'S1' and 'ST' dicts
        title : str
            Plot title
            
        Returns
        -------
        bokeh.plotting.figure
        """
        S1 = sobol_results['S1']
        ST = sobol_results['ST']
        params = sobol_results['ranking']  # Already sorted by importance
        
        # Construct hierarchical categorical coordinates for grouped bars.
        # Each bar is indexed by (parameter, Sobol_index_type).
        # Example: ("alpha_1","S1"), ("beta","S1"), ("alpha_1","ST"), ("beta","ST").
        x = [(p, "S1") for p in params] + [(p, "ST") for p in params]
        source_s1 = ColumnDataSource(
            data=dict(
                x=[(p, 'S1') for p in params],
                counts=[S1[p] for p in params],
                )
            )

        # Corresponding bar heights aligned with the x tuples above.
        # First block contains first-order Sobol indices, second block total indices.
        source_st = ColumnDataSource(
            data=dict(
                x=[(p, 'ST') for p in params],
                counts=[ST[p] for p in params],
                )
            )

        p = bokeh.plotting.figure(
            x_range=FactorRange(*x),
            width=self.width,
            height=self.height,
            title=title,
            x_axis_label="Parameter",
            y_axis_label="Sobol Index",
            tools="save",
        )
        # Color scheme
        colors = ["#3182bd", "#e6550d"]  # Blue for S1, orange for ST
        
        # Add bars
        p.vbar(
            x='x',
            top='counts',
            width=0.9,
            source=source_s1,
            color="#3182bd",
            legend_label="First-order (S1)",
        )

        p.vbar(
            x='x',
            top='counts',
            width=0.9,
            source=source_st,
            color="#e6550d",
            legend_label="Total (ST)",
        )

        # Formatting
        p.y_range = Range1d(0, 1) # sobol indices are theoretically in [0,1]
        p.xaxis.major_label_orientation = 1
        p.xgrid.grid_line_color = None
        p.legend.location = "top_right"

        return p

    def plot_stability_diagram(
        self,
        x_vals: np.ndarray,
        y_vals: np.ndarray,
        stable_condition: Any,
        boundary_fns: Optional[List[Dict[str, Any]]] = None,
        x_name: str = "x",
        y_name: str = "y",
        title: str = "Stability Diagram",
        width: Optional[int] = None,
        height: Optional[int] = None,
        show_sample_points: bool = False,
    ) -> bokeh.plotting.figure:
        """
        Plot a stability diagram with shaded regions and analytical boundary curves.

        Parameters
        ----------
        x_vals : np.ndarray
            1-D array of x-axis parameter values (must be strictly increasing).
        y_vals : np.ndarray
            1-D array of y-axis parameter values (must be strictly increasing).
        stable_condition : callable
            Function ``(X_mesh, Y_mesh) -> bool_array`` where X_mesh and Y_mesh
            are 2-D arrays from ``np.meshgrid(x_vals, y_vals)`` (shape
            ``(len(y_vals), len(x_vals))``).  Returns a boolean array of the
            same shape where ``True`` means stable.
        boundary_fns : list of dict, optional
            Analytical boundary curves to overlay.  Each dict may contain:
            - ``'fn'``: callable ``(x_array) -> y_array``  *(required)*
            - ``'label'``: str legend label (default ``'Boundary'``)
            - ``'color'``: line color (default ``'black'``)
            - ``'dash'``: Bokeh dash style (default ``'dashed'``)
            - ``'line_width'``: float (default ``2.5``)
        x_name, y_name : str
            Axis labels.
        title : str
            Plot title.
        width, height : int, optional
            Override instance defaults.
        show_sample_points : bool
            Overlay the raw grid points.

        Returns
        -------
        bokeh.plotting.figure
        """
        # image_rgba uint32 uses little-endian byte order: bytes are [R, G, B, A]
        # so uint32 = R + G*256 + B*65536 + A*16777216 = 0xAABBGGRR
        _GREEN = np.uint32(0xA58ECF8E)  # #8ECF8E green,  alpha≈0.65
        _GRAY  = np.uint32(0xA5C8C8C8)  # #C8C8C8 gray,   alpha≈0.65

        x_vals = np.asarray(x_vals, dtype=float)
        y_vals = np.asarray(y_vals, dtype=float)

        # Build meshgrid and classify
        X, Y = np.meshgrid(x_vals, y_vals)
        stable_mask = np.asarray(stable_condition(X, Y), dtype=bool)

        # RGBA image (shape: n_y × n_x, uint32)
        rgba = np.where(stable_mask, _GREEN, _GRAY).astype(np.uint32)

        # Pixel edges so the image covers [x0-dx/2, x_end+dx/2]
        dx = float(x_vals[1] - x_vals[0]) if len(x_vals) > 1 else 1.0
        dy = float(y_vals[1] - y_vals[0]) if len(y_vals) > 1 else 1.0
        x0 = float(x_vals[0])  - dx / 2
        y0 = float(y_vals[0])  - dy / 2
        dw = float(x_vals[-1]) - float(x_vals[0]) + dx
        dh = float(y_vals[-1]) - float(y_vals[0]) + dy

        w = width  or self.width
        h = height or self.height

        p = bokeh.plotting.figure(
            width=w, height=h,
            title=title,
            x_axis_label=x_name,
            y_axis_label=y_name,
            x_range=Range1d(float(x_vals[0]), float(x_vals[-1])),
            y_range=Range1d(float(y_vals[0]), float(y_vals[-1])),
            tools="pan,wheel_zoom,box_zoom,reset,save",
        )

        p.image_rgba(image=[rgba], x=x0, y=y0, dw=dw, dh=dh) # type: ignore[arg-type]

        # Proxy legend squares (image_rgba does not participate in the legend)
        p.scatter([], [], marker="square", size=14,
                  fill_color="#8ECF8E", fill_alpha=0.65,
                  line_color=None, legend_label="Stable")
        p.scatter([], [], marker="square", size=14,
                  fill_color="#C8C8C8", fill_alpha=0.65,
                  line_color=None, legend_label="Unstable")

        # Analytical boundary curves
        if boundary_fns:
            for bfn in boundary_fns:
                fn         = bfn["fn"]
                label      = bfn.get("label",      "Boundary")
                color      = bfn.get("color",      "black")
                dash       = bfn.get("dash",       "dashed")
                line_width = bfn.get("line_width",  2.5)

                y_curve = fn(x_vals)
                # Mask points outside the y_vals range
                in_range = (y_curve >= float(y_vals[0])) & (y_curve <= float(y_vals[-1]))
                xs = x_vals[in_range]
                ys = y_curve[in_range]
                if len(xs) >= 2:
                    p.line(xs, ys,
                           line_color=color, line_dash=dash,
                           line_width=line_width, alpha=0.95,
                           legend_label=label)

        if show_sample_points:
            Xc, Yc = np.meshgrid(x_vals, y_vals)
            p.scatter(Xc.ravel(), Yc.ravel(), size=4, color="#444444", alpha=0.30)

        p.add_tools(HoverTool(tooltips=[(x_name, "$x{0.000}"), (y_name, "$y{0.000}")]))
        p.add_tools(CrosshairTool())
        p.legend.location = "top_right"
        p.legend.click_policy = "hide"

        return p

    def plot_root_locus(
        self,
        eigenvalue_paths: Dict[str, Dict],
        title: str = "Root Locus (Dominant Pair)",
        width: Optional[int] = None,
        height: Optional[int] = None,
    ) -> bokeh.plotting.figure:
        """
        Bokeh root locus: eigenvalue paths in the complex plane.

        Points are colored by closed-loop gain (shared Viridis gradient);
        multiple paths are distinguished by line dash style and color.
        Open circle = low-gain start, × = high-gain end.

        Parameters
        ----------
        eigenvalue_paths : dict
            Maps label -> {'real': array, 'imag': array, 'gain': array}.
            'gain' is optional; if present on any path, a shared Viridis
            ColorBar is added and all scatter points are colored by gain.
        title : str
            Plot title.
        width, height : int, optional
            Override instance defaults.
        """
        w = width or self.width
        h = height or self.height

        p = bokeh.plotting.figure(
            width=w, height=h, title=title,
            x_axis_label='Real Part (Re(λ))',
            y_axis_label='Imaginary Part (Im(λ))',
            tools='pan,wheel_zoom,box_zoom,reset,save',
        )

        # Stability boundary Re=0 and real axis
        p.add_layout(Span(location=0, dimension='height',
                          line_color='green', line_dash='dashed',
                          line_width=2.5, line_alpha=0.7))
        p.add_layout(Span(location=0, dimension='width',
                          line_color='black', line_width=0.5, line_alpha=0.3))

        # Build shared gain color mapper if any path has gain data
        paths_with_gain = [ep for ep in eigenvalue_paths.values() if 'gain' in ep]
        if paths_with_gain:
            all_gains = np.concatenate([np.asarray(ep['gain']) for ep in paths_with_gain])
            valid = all_gains[np.isfinite(all_gains)]
            gain_mapper = LinearColorMapper(
                palette=Viridis256,
                low=float(valid.min()),
                high=float(valid.max()),
            )
            color_bar = ColorBar(
                color_mapper=gain_mapper,
                label_standoff=8,
                title='Gain',
                location=(0, 0),
            )
            p.add_layout(color_bar, 'right')
        else:
            gain_mapper = None

        dash_cycle = ['solid', 'dashed', 'dotted', 'dashdot']
        palette = self.palette[max(10, len(eigenvalue_paths))]
        hover_renderers = []

        for idx, (label, eig_path) in enumerate(eigenvalue_paths.items()):
            real_parts = np.asarray(eig_path['real'])
            imag_parts = np.asarray(eig_path['imag'])
            gains = np.asarray(eig_path['gain']) if 'gain' in eig_path else None
            color = palette[idx % len(palette)]
            dash = dash_cycle[idx % len(dash_cycle)]

            # Path line: color + dash distinguish paths
            p.line(real_parts, imag_parts,
                   line_width=2.0, color=color, line_dash=dash,
                   alpha=0.6, legend_label=label)

            # Scatter points: colored by gain if mapper available, else path color
            if gain_mapper is not None and gains is not None:
                src = ColumnDataSource({
                    'x': real_parts, 'y': imag_parts,
                    'gain': gains, 'label': [label] * len(real_parts),
                })
                r = p.scatter('x', 'y', source=src, size=8, alpha=0.9, marker='circle',
                              color={'field': 'gain', 'transform': gain_mapper})
            else:
                src = ColumnDataSource({
                    'x': real_parts, 'y': imag_parts,
                    'label': [label] * len(real_parts),
                })
                r = p.scatter('x', 'y', source=src, size=8,
                              color=color, alpha=0.9, marker='circle')
            hover_renderers.append(r)

            # Start marker: open circle (low gain)
            p.scatter([real_parts[0]], [imag_parts[0]],
                      size=14, color=color, fill_color='white',
                      line_width=2.5, marker='circle')
            # End marker: × (high gain)
            p.scatter([real_parts[-1]], [imag_parts[-1]],
                      size=16, color=color, line_width=2.5, marker='x')

        tooltips = [('Path', '@label'), ('Re(λ)', '@x{0.0000}'), ('Im(λ)', '@y{0.0000}')]
        if gain_mapper is not None:
            tooltips.append(('Gain', '@gain{0.0000}'))
        p.add_tools(HoverTool(renderers=hover_renderers, tooltips=tooltips))
        p.add_tools(CrosshairTool())

        p.legend.location = 'top_left'
        p.legend.click_policy = 'hide'
        return p

    def plot_combined_figure(
        self,
        fig_trajectories: LayoutDOM,
        fig_phase: LayoutDOM,
        fig_overshoot: LayoutDOM,
        fig_sse: LayoutDOM,
        fig_timing: LayoutDOM,
        controller_name: str,
        circuit_diagram_path: Optional[str] = None,
        panel_width: int = 500,
        panel_height: int = 350,
    ) -> LayoutDOM:
        """
        Assemble 6 pre-built Bokeh figures into a 3×2 publication grid.

        Layout::

            [Circuit     ] [fig_trajectories] 
            [fig_phase  ] [fig_sse]
            [fig_overshoot]  [fig_timing ]

        Parameters
        ----------
        fig_trajectories, fig_phase, fig_overshoot, fig_sse, fig_timing
            Pre-built Bokeh figures/layouts (e.g. from plot_metric_interactive).
        controller_name : str
            Used in the header, e.g. 'PI_1'.
        circuit_diagram_path : str, optional
            Path to circuit image (PNG). Placeholder drawn if None.
        panel_width, panel_height : int
            Dimensions for the circuit placeholder panel.
        """
        if circuit_diagram_path is not None:
            circuit_panel: LayoutDOM = Div(
                text=(
                    f'<div style="width:{panel_width}px;height:{panel_height}px;'
                    f'display:flex;align-items:center;justify-content:center;">'
                    f'<img src="{circuit_diagram_path}" '
                    f'style="max-width:100%;max-height:100%;object-fit:contain;"/>'
                    f'</div>'
                ),
                width=panel_width, height=panel_height,
            )
        else:
            circuit_panel = Div(
                text=(
                    f'<div style="width:{panel_width}px;height:{panel_height}px;'
                    f'border:2px solid #888;display:flex;align-items:center;'
                    f'justify-content:center;font-size:14px;color:#888;'
                    f'font-style:italic;text-align:center;">'
                    f'{controller_name}<br>circuit diagram</div>'
                ),
                width=panel_width, height=panel_height,
            )

        header = Div(
            text=(
                f"<h3 style='margin:6px 0'>"
                f"{controller_name} — controller performance summary</h3>"
            )
        )
        grid = bokeh.layouts.gridplot(
            [[circuit_panel,   fig_trajectories, ],
            [fig_phase,  fig_sse,],         
             [fig_overshoot,  fig_timing,]],
            merge_tools=False,
        )
        return bokeh.layouts.column(header, grid)

    def plot_metric_interactive(
        self,
        param_vals: np.ndarray,
        trajectories: Dict[float, Dict],
        metric_name: str,
        gain_vals: np.ndarray,
        metric_vals: np.ndarray,
        metric_stds: Optional[np.ndarray] = None,
        param_name: str = 'α₁',                  
        output_units: str = 'μM',         
        time_units: str = 's',            
        title: str = "Metric vs Gain",
        x_label: str = "Gain",
        x_scale: str = "log",
        y_scale: str = "log",
        return_sources: bool = False,
    ) -> Union[LayoutDOM, Tuple[LayoutDOM, List, List, ColumnDataSource]]:
        """
        Interactive two-panel Bokeh layout: trajectory viewer (left) + metric vs gain (right).

        Clicking a point in the metric panel updates the trajectory panel to show
        only that trajectory.  A Reset button restores the default multi-trajectory view.

        Parameters
        ----------
        param_vals : np.ndarray
            Parameter values swept (one per simulation).
        trajectories : dict
            Maps param -> {'time': array, 'y': array, 'ref': array (optional)}.
        metric_name : str
            Metric label (e.g. 'ss_error', 'settling_time').
        gain_vals : np.ndarray
            Gain value for each param (x-axis of metric panel).
        metric_vals : np.ndarray
            Metric value for each param (y-axis of metric panel).
        title : str
            Overall layout title.
        x_label : str
            X-axis label for the metric panel.
        x_scale, y_scale : str
            'linear' or 'log' for metric panel axes.
        """
        panel_w = self.width
        panel_h = self.height

        n = len(param_vals)
        hex_colors = (
            [Viridis256[128]] if n == 1
            else [Viridis256[int(round(i * 255 / (n - 1)))] for i in range(n)]
        )

        # Keep only values for which a trajectory was recorded
        valid = [(i, float(a)) for i, a in enumerate(param_vals) if float(a) in trajectories]
        v_idx     = [i for i, _ in valid]
        print(v_idx)
        v_values  = [a for _, a in valid]
        v_gains   = [float(gain_vals[i]) for i in v_idx]
        v_metrics = [float(metric_vals[i]) for i in v_idx]
        v_colors  = [hex_colors[i] for i in v_idx]

        # All trajectory data as plain Python lists (needed for CustomJS serialisation)
        v_times = [trajectories[a]['time'].tolist() for a in v_values]
        v_ys    = [trajectories[a]['y'].tolist()    for a in v_values]
        v_refs  = [
            trajectories[a]['ref'].tolist() if 'ref' in trajectories[a]
            else [0.0] * len(trajectories[a]['time'])
            for a in v_values
        ]
        v_stds = [float(metric_stds[i]) if metric_stds is not None else 0.0 for i in v_idx]


        # ── Metric panel source ───────────────────────────────────────────────
        metric_source = ColumnDataSource(data={
            'x':     v_gains,
            'y':     v_metrics,
            'label': [f'{param_name}={a:.2f}' for a in v_values],
            'color': list(v_colors),
            'std':    v_stds,
            'y_upper': [m + s for m, s in zip(v_metrics, v_stds)],
            'y_lower': [max(m - s, 0.0) for m, s in zip(v_metrics, v_stds)],
        })

        # ── LEFT: Trajectory panel ────────────────────────────────────────────
        # One line renderer per param so Bokeh's native Legend can be used and
        # toggling renderer.visible automatically hides/shows the legend item.
        traj_fig = bokeh.plotting.figure(
            width=panel_w, height=panel_h,
            title='Trajectories',
            x_axis_label=f'Time ({time_units})',
            y_axis_label=f'Concentration y ({output_units})',
            tools='pan,wheel_zoom,box_zoom,reset,save',
        )

        traj_renderers: List = []
        ref_renderers:  List = []
        legend_items:   List = []

        for k, (param, color, ts, ys, rs) in enumerate(
                zip(v_values, v_colors, v_times, v_ys, v_refs)):
            src = ColumnDataSource({'t': ts, 'y': ys, 'r': rs})
            # fade out higher params (which are drawn earlier, behind)
            alpha = 1.0 - 0.45 * (k / max(len(v_values) - 1, 1))
            
            r_traj = traj_fig.line('t', 'y', source=src,
                                line_color=color, line_width=2.0, line_alpha=alpha)
            r_ref  = traj_fig.line('t', 'r', source=src,
                                   line_color='black', line_alpha=0.45,
                                   line_width=1.2, line_dash='dashed')
            traj_renderers.append(r_traj)
            ref_renderers.append(r_ref)
            legend_items.append(LegendItem(label=f'{param_name}={param:.2f}', 
                                        renderers=[cast(GlyphRenderer, r_traj)]))

        traj_legend = Legend(items=legend_items, location='top_left', click_policy='hide')
        traj_fig.add_layout(traj_legend)
        traj_fig.add_tools(
            CrosshairTool(),
            HoverTool(tooltips=[('Time', '$x{0.00}'), ('y', '$y{0.0000}')]),
        )

        # ── RIGHT: Metric panel ───────────────────────────────────────────────
        metric_fig = bokeh.plotting.figure(
            width=panel_w, height=panel_h,
            title=f"{metric_name.replace('_', ' ').title()} vs Gain",
            x_axis_label=x_label,
            y_axis_label=metric_name.replace('_', ' ').title(),
            x_axis_type='log' if x_scale == 'log' else 'linear',
            y_axis_type='log' if y_scale == 'log' else 'linear',
            tools='pan,wheel_zoom,box_zoom,reset,save',
        )
        metric_fig.line('x', 'y', source=metric_source,
                        line_color='navy', line_width=2, alpha=0.4)
        scatter = metric_fig.scatter(
            'x', 'y', source=metric_source,
            size=10, color='color', line_color='navy', line_width=1,
            alpha=0.85,
            nonselection_alpha=0.3, nonselection_line_color='navy',
            marker='circle',
        )
        if metric_stds is not None:
            # # shaded band for ±1 std
            # band = Band(
            #     base='x', lower='y_lower', upper='y_upper',
            #     source=metric_source,
            #     fill_alpha=0.2, fill_color='navy',
            #     line_color='navy', line_alpha=0.4, line_width=0.5,
            # )
            # metric_fig.add_layout(band)
            
            # whisker error bars
            whisker = Whisker(
                base='x', upper='y_upper', lower='y_lower',
                source=metric_source,
                level='overlay',
                line_color='navy', line_alpha=0.6,
            )
            metric_fig.add_layout(whisker)
    
        metric_fig.add_tools(
            CrosshairTool(),
            TapTool(renderers=[scatter]),
            HoverTool(
                renderers=[scatter],
                tooltips=[
                    ('α₁',  '@label'),
                    ('Gain', '@x{0.0000}'),
                    (metric_name.replace('_', ' ').title(), '@y{0.0000}'),
                ],
            ),
        )

        # ── CustomJS: tap → hide all renderers except selected ────────────────
        tap_cb = CustomJS(args=dict(
            metric_source=metric_source,
            traj_renderers=traj_renderers,
            ref_renderers=ref_renderers,
        ), code="""
            const sel = metric_source.selected.indices;
            if (sel.length === 0) return;
            const i = sel[0];
            for (let k = 0; k < traj_renderers.length; k++) {
                traj_renderers[k].visible = (k === i);
                ref_renderers[k].visible  = (k === i);
            }
        """)
        metric_source.selected.js_on_change('indices', tap_cb)

        # ── Reset button ──────────────────────────────────────────────────────
        reset_btn = Button(label='↺  Reset View', button_type='default', width=140)
        reset_cb = CustomJS(args=dict(
            metric_source=metric_source,
            traj_renderers=traj_renderers,
            ref_renderers=ref_renderers,
        ), code="""
            for (let k = 0; k < traj_renderers.length; k++) {
                traj_renderers[k].visible = true;
                ref_renderers[k].visible  = true;
            }
            metric_source.selected.indices = [];
            metric_source.change.emit();
        """)
        reset_btn.js_on_click(reset_cb)

        # ── Assemble layout ───────────────────────────────────────────────────
        header = Div(text=f"<h3 style='margin:4px 0'>{title}</h3>")
        layout = bokeh.layouts.column(
            header,
            bokeh.layouts.row(traj_fig, metric_fig),
            reset_btn,
        )
        if return_sources:
            return layout, traj_renderers, ref_renderers, metric_source
        return layout
