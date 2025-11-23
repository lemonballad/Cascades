"""
Visualization functions for 2D Raman spectroscopy data.

Provides plotting utilities for cascade and direct signals, spectral data,
and cascade-to-direct ratio maps.
"""

import numpy as np
from numpy.typing import NDArray
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from typing import Optional


def plot_2d_spectrum(
    data: NDArray[np.complex128],
    freq: NDArray[np.float64],
    ax: Optional[Axes] = None,
    title: str = "",
    xlabel: str = r"$\omega_1/2\pi c$ (cm$^{-1}$)",
    ylabel: str = r"$\omega_2/2\pi c$ (cm$^{-1}$)",
    n_contours: int = 50,
    cmap: str = "jet",
    normalize: bool = True
) -> Axes:
    """
    Plot a 2D contour spectrum.

    Parameters
    ----------
    data : NDArray[np.complex128]
        2D spectral data (complex, will plot absolute value).
    freq : NDArray[np.float64]
        Frequency axis in cm^-1.
    ax : Axes, optional
        Matplotlib axes to plot on. Creates new figure if None.
    title : str, optional
        Plot title.
    xlabel : str, optional
        X-axis label.
    ylabel : str, optional
        Y-axis label.
    n_contours : int, optional
        Number of contour levels (default 50).
    cmap : str, optional
        Colormap name (default "jet").
    normalize : bool, optional
        Normalize data to maximum (default True).

    Returns
    -------
    Axes
        The matplotlib axes with the plot.
    """
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(6, 6))

    plot_data = np.abs(data)
    if normalize:
        plot_data = plot_data / np.max(plot_data)

    cs = ax.contour(freq, freq, plot_data, n_contours, cmap=cmap)
    ax.set_xlabel(xlabel, fontsize=10, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=10, fontweight='bold')
    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.set_aspect('equal')
    plt.colorbar(cs, ax=ax)

    return ax


def plot_ratio_map(
    ratio: NDArray[np.float64],
    detunings: NDArray[np.float64],
    displacements: NDArray[np.float64],
    ax: Optional[Axes] = None,
    title: str = "Cascade/Direct Ratio",
    xlabel: str = r"$(\omega_L - \omega_{eg})/2\pi c$ (cm$^{-1}$)",
    ylabel: str = "Mode Displacement $d$",
    n_contours: int = 50,
    cmap: str = "jet"
) -> Axes:
    """
    Plot cascade-to-direct ratio as a filled contour map.

    Parameters
    ----------
    ratio : NDArray[np.float64]
        2D array of cascade/direct ratios.
    detunings : NDArray[np.float64]
        Laser detuning values (cm^-1).
    displacements : NDArray[np.float64]
        Mode displacement values.
    ax : Axes, optional
        Matplotlib axes to plot on. Creates new figure if None.
    title : str, optional
        Plot title.
    xlabel : str, optional
        X-axis label.
    ylabel : str, optional
        Y-axis label.
    n_contours : int, optional
        Number of contour levels (default 50).
    cmap : str, optional
        Colormap name (default "jet").

    Returns
    -------
    Axes
        The matplotlib axes with the plot.
    """
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(6, 6))

    cs = ax.contourf(
        detunings, displacements, ratio,
        n_contours, cmap=cmap
    )
    ax.set_xlabel(xlabel, fontsize=10, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=10, fontweight='bold')
    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.set_aspect('equal')
    plt.colorbar(cs, ax=ax)

    return ax


def plot_2drr_comparison(
    cascade: NDArray[np.complex128],
    direct: NDArray[np.complex128],
    ratio: NDArray[np.float64],
    freq: NDArray[np.float64],
    detunings: NDArray[np.float64],
    displacements: NDArray[np.float64],
    figsize: tuple[float, float] = (15, 5)
) -> Figure:
    """
    Create a 3-panel comparison plot of cascade, direct, and ratio data.

    Parameters
    ----------
    cascade : NDArray[np.complex128]
        2D cascade signal data.
    direct : NDArray[np.complex128]
        2D direct signal data.
    ratio : NDArray[np.float64]
        2D cascade/direct ratio map.
    freq : NDArray[np.float64]
        Frequency axis in cm^-1.
    detunings : NDArray[np.float64]
        Laser detuning values (cm^-1).
    displacements : NDArray[np.float64]
        Mode displacement values.
    figsize : tuple, optional
        Figure size (width, height) in inches.

    Returns
    -------
    Figure
        The matplotlib figure with three subplots.
    """
    fig, axes = plt.subplots(1, 3, figsize=figsize)

    # Normalize to common scale
    max_signal = max(np.max(np.abs(cascade)), np.max(np.abs(direct)))

    plot_2d_spectrum(
        cascade / max_signal, freq, ax=axes[0],
        title="Cascade", normalize=False
    )
    plot_2d_spectrum(
        direct / max_signal, freq, ax=axes[1],
        title="Direct", normalize=False
    )
    plot_ratio_map(
        ratio, detunings, displacements, ax=axes[2]
    )

    fig.tight_layout()
    return fig


def save_figure(
    fig: Figure,
    filename: str,
    dpi: int = 300,
    format: str = "tif"
) -> None:
    """
    Save figure to file.

    Parameters
    ----------
    fig : Figure
        Matplotlib figure to save.
    filename : str
        Output filename (without extension).
    dpi : int, optional
        Resolution in dots per inch (default 300).
    format : str, optional
        Output format: tif, png, pdf, etc. (default "tif").
    """
    fig.savefig(f"{filename}.{format}", dpi=dpi, bbox_inches='tight')
