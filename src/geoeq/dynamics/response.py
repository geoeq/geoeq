"""
Response-spectrum plotting helpers.

A minimal utility: given period T and spectral acceleration Sa arrays
(typically from a code/standard such as ASCE 7 or IBC), draws a clean
log-log response spectrum plot.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np


def response_spectrum(T: Sequence[float], Sa: Sequence[float],
                       site_class: str = None, ax=None,
                       save_as: str = None, **kwargs):
    """Plot a response spectrum.

    Parameters
    ----------
    T : sequence
        Periods (s).
    Sa : sequence
        Spectral accelerations (g).
    site_class : str, optional
        Label printed on the figure (e.g. 'Site Class D').
    ax : matplotlib axis, optional
        Embed in an existing axis.
    save_as : str, optional
        Path to save the figure.
    **kwargs
        Passed to Matplotlib ``plot``.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    T = np.asarray(T)
    Sa = np.asarray(Sa)
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))
    else:
        fig = ax.figure
    style = {"color": "#1f3a93", "linewidth": 1.8}
    style.update(kwargs)
    ax.plot(T, Sa, **style)
    ax.set_xlabel("Period $T$ (s)")
    ax.set_ylabel("Spectral acceleration $S_a$ (g)")
    title = "Response spectrum"
    if site_class:
        title += f" ({site_class})"
    ax.set_title(title)
    ax.set_xscale("log")
    ax.grid(True, which="both", alpha=0.3)
    if save_as:
        fig.savefig(save_as, dpi=300, bbox_inches="tight")
    return fig
