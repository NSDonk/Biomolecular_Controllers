"""
Corrected analysis plots for biomolecular controllers.

Implements:
1-4. Metric vs Gain (2-panel: trajectories + metric plot)
5-7. Stability Phase Diagrams (2D parameter sweeps)
8. Root Locus (eigenvalue paths in complex plane)
"""

from typing import Dict, List, Optional, Tuple, Callable
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib import colormaps


class MetricPlotter:
    """Create trajectory + metric plots for controller analysis."""
    
    @staticmethod
    def plot_metric_with_trajectories(
        vals: np.ndarray,
        trajectories: Dict[float, Dict],  # {alpha: {'time': ..., 'y': ..., 'metric': ...}}
        metric_name: str,
        gain_vals: np.ndarray,
        metric_vals: np.ndarray,
        title: str = "Metric vs Gain",
        figsize: Tuple[int, int] = (16, 5),
        x_scale: str = 'log',
        y_scale: str = 'log',
        x_label: Optional[str] = None,
    ) -> Figure:
        """
        Two-panel plot: left shows trajectories, right shows metric vs gain.
        
        Parameters
        ----------
        vals : np.ndarray
            Alpha values used in sweep
        trajectories : dict
            Maps alpha -> {'time': array, 'y': array, 'ref': array, 'metric': float}
        metric_name : str
            Name of metric ('ss_error', 'settling_time', 'overshoot', 'rise_time')
        gain_vals : np.ndarray
            Computed gain values (G = α₁·θ₁·θ₂)
        metric_vals : np.ndarray
            Computed metric values
        title : str
            Figure title
        figsize : tuple
            Figure size
            
        Returns
        -------
        matplotlib.figure.Figure
        """
        fig, (ax_traj, ax_metric) = plt.subplots(1, 2, figsize=figsize)
        
        # LEFT PANEL: Trajectories
        cmap = colormaps["tab10"]
        colors = cmap(np.linspace(0, 1, len(vals)))
        
        for i, alpha in enumerate(vals):
            if alpha in trajectories:
                traj = trajectories[alpha]
                ax_traj.plot(traj['time'], traj['y'], 
                           color=colors[i], linewidth=2, alpha=0.7,
                           label=f'α₁={alpha:.1f}')
                
                # Plot reference too
                if 'ref' in traj:
                    ax_traj.plot(traj['time'], traj['ref'], 
                               color=colors[i], linestyle='--', 
                               linewidth=1, alpha=0.4)
        
        ax_traj.set_xlabel('Time', fontsize=12, fontweight='bold')
        ax_traj.set_ylabel('Output (y)', fontsize=12, fontweight='bold')
        ax_traj.set_title('Trajectories at Different Gains', fontsize=13, fontweight='bold')
        ax_traj.legend(loc='best', fontsize=9)
        ax_traj.grid(True, alpha=0.3)
        
        # RIGHT PANEL: Metric vs Gain
        # When using log scale, shift any zero gain values by 1e-8 to avoid log(0)
        plot_gains = gain_vals.copy().astype(float)
        if x_scale == 'log':
            plot_gains[plot_gains == 0] += 1e-8

        ax_metric.plot(plot_gains, metric_vals, 'o-', linewidth=2.5,
                       markersize=8, color='navy', alpha=0.7)
        ax_metric.set_xscale(x_scale)
        ax_metric.set_yscale(y_scale)

        xlabel = x_label if x_label is not None else 'Gain'
        ax_metric.set_xlabel(xlabel, fontsize=12, fontweight='bold')
        ax_metric.set_ylabel(metric_name.replace('_', ' ').title(), fontsize=12, fontweight='bold')
        ax_metric.set_title(f'{metric_name.replace("_", " ").title()} vs Gain',
                           fontsize=13, fontweight='bold')
        grid_kwargs = {'alpha': 0.3, 'which': 'both'} if (x_scale == 'log' or y_scale == 'log') else {'alpha': 0.3}
        ax_metric.grid(True, **grid_kwargs)
        
        fig.suptitle(title, fontsize=15, fontweight='bold', y=1.00)
        fig.tight_layout()
        
        return fig
    
    @staticmethod
    def plot_stability_phase_diagram(
        x_vals: np.ndarray,
        y_vals: np.ndarray,
        stability_grid: np.ndarray,  # 2D array: shape (len(y_vals), len(x_vals))
        x_name: str = "α₁",
        y_name: str = "θ₁",
        title: str = "Stability Phase Diagram",
        figsize: Tuple[int, int] = (10, 8),
    ) -> Figure:
        """
        2D stability phase diagram with stable/unstable regions.
        
        Parameters
        ----------
        x_vals : np.ndarray
            X-axis parameter values
        y_vals : np.ndarray
            Y-axis parameter values
        stability_grid : np.ndarray
            2D array of λ_max values, shape (len(y_vals), len(x_vals))
        x_name, y_name : str
            Parameter names
        title : str
            Plot title
        figsize : tuple
            Figure size
            
        Returns
        -------
        matplotlib.figure.Figure
        """
        fig, ax = plt.subplots(figsize=figsize)
        
        X, Y = np.meshgrid(x_vals, y_vals)
        
        # Create custom colormap: red for unstable, green for stable
        # Values < 0 are stable (green), > 0 are unstable (red)
        vmin = stability_grid.min()
        vmax = stability_grid.max()
        vcenter = 0
        
        # Plot stable region (green)
        stable_mask = stability_grid < 0
        ax.contourf(X, Y, stability_grid, levels=[-1e6, 0], colors=['green'], alpha=0.3)
        
        # Plot unstable region (red)
        unstable_mask = stability_grid >= 0
        ax.contourf(X, Y, stability_grid, levels=[0, 1e6], colors=['red'], alpha=0.3)
        
        # Plot stability boundary (λ_max = 0)
        cs = ax.contour(X, Y, stability_grid, levels=[0], colors=['black'], linewidths=2.5)
        ax.clabel(cs, inline=True, fontsize=10, fmt='Boundary')
        
        # Labels and scaling
        ax.set_xlabel(x_name, fontsize=12, fontweight='bold')
        ax.set_ylabel(y_name, fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3, linestyle='--')
        
        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='green', alpha=0.3, label='Stable'),
            Patch(facecolor='red', alpha=0.3, label='Unstable')
        ]
        ax.legend(handles=legend_elements, loc='best', fontsize=11)
        
        fig.tight_layout()
        return fig
    
    @staticmethod
    def plot_root_locus(
        gain_vals: np.ndarray,
        eigenvalue_paths: Dict[str, Dict],  # {eig_label: {'real': [...], 'imag': [...]}}
        title: str = "Root Locus (Dominant Pair)",
        figsize: Tuple[int, int] = (10, 8),
    ) -> Figure:
        """
        Plot root locus showing eigenvalue paths in complex plane.
        
        Parameters
        ----------
        gain_vals : np.ndarray
            Gain values along path
        eigenvalue_paths : dict
            Maps eigenvalue identifier -> {'real': array, 'imag': array}
            (One entry for dominant eigenvalue pair)
        title : str
            Plot title
        figsize : tuple
            Figure size
            
        Returns
        -------
        matplotlib.figure.Figure
        """
        fig, ax = plt.subplots(figsize=figsize)
        
        # Plot stability boundary (Re = 0)
        y_lim = [-3, 3]  # Will adjust based on data
        ax.axvline(x=0, color='green', linestyle='--', linewidth=2.5, 
                  label='Stability Boundary (Re=0)', alpha=0.7)
        
        # Plot eigenvalue paths
        cmap = colormaps["tab10"]
        colors = cmap(np.linspace(0, 1, len(eigenvalue_paths)))
        
        # setting default y-axis limits
        y_min = -3.0
        y_max = 3.0
        
        for idx, (label, eig_path) in enumerate(eigenvalue_paths.items()):
            real_parts = eig_path['real']
            imag_parts = eig_path['imag']
            
            # Plot path
            ax.plot(real_parts, imag_parts, 'o-', linewidth=2.5, 
                   markersize=6, color=colors[idx], alpha=0.7, label=label)
            
            # Mark start (low gain) with 'o'
            ax.plot(real_parts[0], imag_parts[0], 'o', markersize=10, 
                   color=colors[idx], markerfacecolor='white', markeredgewidth=2)
            
            # Mark end (high gain) with 'x'
            ax.plot(real_parts[-1], imag_parts[-1], 'x', markersize=12, 
                   color=colors[idx], markeredgewidth=2.5)
            
            # Update y limits
            if len(imag_parts) > 0:
                y_min = min(y_lim[0], imag_parts.min() - 0.5)
                y_max = max(y_lim[1], imag_parts.max() + 0.5)

        
        ax.set_xlabel('Real Part (Re(λ))', fontsize=12, fontweight='bold')
        ax.set_ylabel('Imaginary Part (Im(λ))', fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best', fontsize=11)
        ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)
        ax.set_ylim(y_min, y_max)
        
        fig.tight_layout()
        return fig
