"""
SoilProfile: layered geotechnical profile with stress calculations.

A profile is an ordered list of ``(top, bot, Soil)`` tuples plus an
optional water table depth. Coordinates increase downward from the
ground surface (z = 0).

Stress calculation rules (Das 2010, Ch. 5)
------------------------------------------
For a depth z:

1. Total stress sigma = integral of gamma * dz from 0 to z, using
   ``gamma`` above the water table and ``gamma_sat`` below it.
2. Pore pressure u = gamma_w * max(0, z - z_w).
3. Effective stress sigma' = sigma - u.

If the water table is above the ground surface (z_w < 0) it is treated
as flooding at the surface, contributing extra hydrostatic pressure.
If z_w is None or +inf the profile is assumed dry throughout.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np

from geoeq.core.constants import GAMMA_WATER
from geoeq.profile.soil import Soil

# A layer is (top, bottom, Soil) -- depths in m, positive downward.
LayerInput = Tuple[float, float, Soil]


@dataclass
class _Layer:
    top: float
    bot: float
    soil: Soil

    @property
    def thickness(self) -> float:
        return self.bot - self.top


class SoilProfile:
    """Layered soil profile with stress computations and plotting.

    Parameters
    ----------
    layers : sequence of (top, bot, Soil)
        Layers in order of increasing depth. Layers must be contiguous
        and non-overlapping. Tops and bottoms in metres.
    water_table : float, optional
        Depth of the phreatic surface (m, positive downward).
        ``None`` or ``np.inf`` means no water table (dry profile).

    Examples
    --------
    >>> from geoeq.profile import Soil, SoilProfile
    >>> p = SoilProfile([
    ...     (0, 2, Soil("Fill",       gamma=18)),
    ...     (2, 8, Soil("Soft Clay",  gamma=17, gamma_sat=18.5)),
    ...     (8, 20, Soil("Dense Sand", gamma=19, gamma_sat=20.5)),
    ... ], water_table=2.0)
    >>> round(p.effective_stress(10), 1)
    133.7
    """

    def __init__(
        self,
        layers: Sequence[LayerInput],
        water_table: Optional[float] = None,
    ):
        if not layers:
            raise ValueError("SoilProfile needs at least one layer.")
        self._layers: List[_Layer] = []
        prev_bot = layers[0][0]
        for i, (top, bot, soil) in enumerate(layers):
            if bot <= top:
                raise ValueError(
                    f"Layer {i}: bottom ({bot}) must be > top ({top}).")
            if not np.isclose(top, prev_bot):
                raise ValueError(
                    f"Layer {i} starts at {top} m but previous layer ended "
                    f"at {prev_bot} m -- profile must be contiguous.")
            if not isinstance(soil, Soil):
                raise TypeError(
                    f"Layer {i}: third element must be a Soil instance.")
            self._layers.append(_Layer(float(top), float(bot), soil))
            prev_bot = bot
        self.water_table = (
            float("inf") if water_table is None else float(water_table))

    # ---------------------------------------------------------------
    # Introspection
    # ---------------------------------------------------------------
    @property
    def top(self) -> float:
        return self._layers[0].top

    @property
    def bottom(self) -> float:
        return self._layers[-1].bot

    @property
    def n_layers(self) -> int:
        return len(self._layers)

    def layers(self) -> List[Tuple[float, float, Soil]]:
        """Return layers as a list of ``(top, bot, Soil)`` tuples."""
        return [(L.top, L.bot, L.soil) for L in self._layers]

    def layer_at(self, z: float) -> Soil:
        """Return the ``Soil`` instance at depth z (m)."""
        z = float(z)
        if z < self.top or z > self.bottom:
            raise ValueError(
                f"Depth {z} m is outside the profile "
                f"({self.top}..{self.bottom}).")
        for L in self._layers:
            if L.top <= z <= L.bot:
                return L.soil
        raise ValueError(f"Depth {z} m not found in any layer.")

    def add_layer(self, top: float, bot: float, soil: Soil) -> None:
        """Append a layer at the bottom of the profile."""
        if not np.isclose(top, self.bottom):
            raise ValueError(
                f"New layer must start at current bottom "
                f"({self.bottom} m), got {top} m.")
        if bot <= top:
            raise ValueError(f"bot ({bot}) must be > top ({top}).")
        self._layers.append(_Layer(float(top), float(bot), soil))

    # ---------------------------------------------------------------
    # Stress calculations  (Das Ch. 5, Terzaghi 1943)
    # ---------------------------------------------------------------
    def total_stress(self, z: Union[float, Iterable[float]]) -> Union[float, np.ndarray]:
        """Total vertical stress sigma_v at depth z (kPa).

        Notes
        -----
        Integrates ``gamma`` above water table and ``gamma_sat`` below.
        If the water table is above the ground surface, hydrostatic
        pressure of the standing water is added.
        """
        z_arr = np.atleast_1d(np.asarray(z, dtype=float))
        out = np.zeros_like(z_arr)
        for i, zi in enumerate(z_arr):
            out[i] = self._sigma_at(float(zi))
        return float(out[0]) if np.isscalar(z) else out

    def _sigma_at(self, z: float) -> float:
        if z < self.top:
            raise ValueError(
                f"Depth {z} m is above profile top ({self.top} m).")
        if z > self.bottom + 1e-9:
            raise ValueError(
                f"Depth {z} m is below profile bottom ({self.bottom} m).")
        zw = self.water_table
        sigma = 0.0
        # Surface flooding contribution.
        if zw < self.top:
            sigma += GAMMA_WATER * (self.top - zw)
        for L in self._layers:
            if z <= L.top:
                break
            z_in_layer = min(z, L.bot) - L.top
            # Portion of this layer above water table:
            above = max(0.0, min(zw, L.bot) - L.top)
            below = max(0.0, z_in_layer - above)
            # Cap by depth-in-layer if z lands inside this layer.
            above = min(above, z_in_layer)
            sigma += L.soil.gamma * above + float(L.soil.gamma_sat) * below
            if z <= L.bot:
                break
        return sigma

    def pore_pressure(
        self, z: Union[float, Iterable[float]]
    ) -> Union[float, np.ndarray]:
        """Hydrostatic pore water pressure u at depth z (kPa)."""
        z_arr = np.atleast_1d(np.asarray(z, dtype=float))
        u = GAMMA_WATER * np.maximum(0.0, z_arr - self.water_table)
        return float(u[0]) if np.isscalar(z) else u

    def effective_stress(
        self, z: Union[float, Iterable[float]]
    ) -> Union[float, np.ndarray]:
        """Effective vertical stress sigma'_v = sigma_v - u (kPa)."""
        return self.total_stress(z) - self.pore_pressure(z)

    def stress_at(self, z: float) -> dict:
        """Return ``{'sigma': ..., 'u': ..., 'sigma_eff': ...}`` at depth z."""
        return {
            "sigma": float(self.total_stress(z)),
            "u": float(self.pore_pressure(z)),
            "sigma_eff": float(self.effective_stress(z)),
        }

    # ---------------------------------------------------------------
    # Export / plotting
    # ---------------------------------------------------------------
    def to_dataframe(self):
        """Export layer data to a pandas DataFrame (if pandas is installed)."""
        try:
            import pandas as pd
        except ImportError as e:  # pragma: no cover -- soft dep
            raise ImportError("pandas is required for to_dataframe()") from e
        rows = []
        for L in self._layers:
            rows.append({
                "name": L.soil.name,
                "top": L.top,
                "bot": L.bot,
                "thickness": L.thickness,
                "gamma": L.soil.gamma,
                "gamma_sat": L.soil.gamma_sat,
                "phi": L.soil.phi,
                "c": L.soil.c,
            })
        return pd.DataFrame(rows)

    def plot(self, dz: float = 0.1, ax=None, show: bool = False, save_as=None):
        """Plot sigma, u, sigma' vs depth.

        Returns the Matplotlib figure for further customization.
        """
        import matplotlib.pyplot as plt
        depths = np.arange(self.top, self.bottom + dz / 2, dz)
        sigma = self.total_stress(depths)
        u = self.pore_pressure(depths)
        sigma_eff = self.effective_stress(depths)

        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 8))
        else:
            fig = ax.figure
        ax.plot(sigma, depths, label=r"$\sigma$ (total)",
                color="#1f3a93", linewidth=1.8)
        ax.plot(u, depths, label=r"$u$ (pore water)",
                color="#2980b9", linewidth=1.8, linestyle="--")
        ax.plot(sigma_eff, depths, label=r"$\sigma'$ (effective)",
                color="#c0392b", linewidth=1.8)

        # Water table line.
        if np.isfinite(self.water_table):
            ax.axhline(self.water_table, color="#3498db",
                       linewidth=1.0, linestyle=":")
            ax.text(0.02, self.water_table, "  WT",
                    transform=ax.get_yaxis_transform(),
                    va="center", color="#2980b9", fontsize=9)
        # Layer boundaries.
        for L in self._layers[1:]:
            ax.axhline(L.top, color="0.7", linewidth=0.6)

        ax.invert_yaxis()
        ax.set_xlabel("Stress (kPa)")
        ax.set_ylabel("Depth (m)")
        ax.set_title("Stress profile")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="lower right")

        if save_as:
            fig.savefig(save_as, dpi=300, bbox_inches="tight")
        if show:  # pragma: no cover -- interactive
            plt.show()
        return fig

    # ---------------------------------------------------------------
    # Iteration
    # ---------------------------------------------------------------
    def __iter__(self):
        for L in self._layers:
            yield (L.top, L.bot, L.soil)

    def __len__(self) -> int:
        return len(self._layers)

    def __repr__(self) -> str:  # pragma: no cover -- cosmetic
        return (f"SoilProfile({self.n_layers} layers, "
                f"{self.top}..{self.bottom} m, "
                f"WT={self.water_table})")


# ---------------------------------------------------------------------
# Free functions
# ---------------------------------------------------------------------
def mesh(profile: SoilProfile, dz: float = 0.5) -> np.ndarray:
    """Calculation grid of depths from profile top to bottom at spacing dz."""
    if dz <= 0:
        raise ValueError("dz must be positive.")
    return np.arange(profile.top, profile.bottom + dz / 2, dz)


def log_plot(boreholes, save_as=None):  # pragma: no cover -- light stub
    """Multi-borehole log plot.

    Parameters
    ----------
    boreholes : dict[str, SoilProfile]
        Mapping of borehole label to profile.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(boreholes), figsize=(3 * len(boreholes), 8),
                             sharey=True)
    if len(boreholes) == 1:
        axes = [axes]
    for ax, (name, profile) in zip(axes, boreholes.items()):
        for top, bot, soil in profile:
            ax.fill_betweenx([top, bot], 0, 1, alpha=0.4,
                             label=soil.name)
            ax.text(0.5, (top + bot) / 2, soil.name,
                    ha="center", va="center", fontsize=8)
        ax.set_xlim(0, 1)
        ax.set_xticks([])
        ax.set_title(name)
        ax.invert_yaxis()
    axes[0].set_ylabel("Depth (m)")
    if save_as:
        fig.savefig(save_as, dpi=300, bbox_inches="tight")
    return fig
