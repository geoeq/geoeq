r"""
Permeability (hydraulic conductivity) laboratory tests for GeoEq.

Functions
---------
constant_head
    Constant-head permeability test → k (ASTM D2434).
falling_head
    Falling-head permeability test → k (ASTM D5084).
permeability_plot
    Plot flow rate vs head gradient or head vs time.

References
----------
Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed., Ch. 5.
ASTM D2434 (2006). Standard Test for Permeability of Granular Soils.
ASTM D5084 (2016). Standard Test for Hydraulic Conductivity.
"""

from typing import Dict, List, Optional, Union
import numpy as np
from geoeq.core.validation import check_positive


def constant_head(
    Q: Union[float, np.ndarray],
    L: float,
    A: float,
    h: float,
    t: float = 1.0,
) -> Union[float, np.ndarray]:
    r"""
    Compute hydraulic conductivity from a constant-head permeability test.

    .. math::

        k = \frac{Q \, L}{A \, h \, t}
        \qquad \text{[Das Eq.\,5.11]}

    Parameters
    ----------
    Q : float or array_like
        Volume of water collected (cm³).
    L : float
        Length of soil specimen (cm).
    A : float
        Cross-sectional area of specimen (cm²).
    h : float
        Constant head difference (cm).
    t : float, default 1.0
        Duration of flow collection (s).

    Returns
    -------
    float or ndarray
        Hydraulic conductivity *k* (cm/s).

    References
    ----------
    Das (2021), Eq. 5.11; ASTM D2434.

    Examples
    --------
    >>> from geoeq.lab.permeability import constant_head
    >>> round(constant_head(Q=500, L=15, A=30, h=50, t=120), 6)
    0.041667
    """
    Q_arr = np.asarray(Q, dtype=float)
    check_positive(L, "L")
    check_positive(A, "A")
    check_positive(h, "h")
    check_positive(t, "t")
    check_positive(Q_arr, "Q")

    k = (Q_arr * L) / (A * h * t)

    if np.ndim(Q) == 0:
        return float(k)
    return np.asarray(k)


def falling_head(
    a: float,
    L: float,
    A: float,
    h1: float,
    h2: float,
    t: float,
) -> float:
    r"""
    Compute hydraulic conductivity from a falling-head permeability test.

    .. math::

        k = \frac{a \, L}{A \, t} \ln\!\left(\frac{h_1}{h_2}\right)
        \qquad \text{[Das Eq.\,5.13]}

    Parameters
    ----------
    a : float
        Cross-sectional area of the standpipe (cm²).
    L : float
        Length of soil specimen (cm).
    A : float
        Cross-sectional area of specimen (cm²).
    h1 : float
        Initial head in the standpipe (cm).
    h2 : float
        Final head in the standpipe (cm), must be < h1.
    t : float
        Elapsed time (s).

    Returns
    -------
    float
        Hydraulic conductivity *k* (cm/s).

    References
    ----------
    Das (2021), Eq. 5.13; ASTM D5084.

    Examples
    --------
    >>> from geoeq.lab.permeability import falling_head
    >>> k = falling_head(a=1.0, L=15, A=30, h1=100, h2=50, t=600)
    >>> round(k, 7)
    5.78e-04
    """
    check_positive(a, "a")
    check_positive(L, "L")
    check_positive(A, "A")
    check_positive(h1, "h1")
    check_positive(h2, "h2")
    check_positive(t, "t")
    if h2 >= h1:
        raise ValueError(f"h2 must be less than h1, got h1={h1}, h2={h2}.")

    k = (a * L) / (A * t) * np.log(h1 / h2)
    return float(k)


def permeability_plot(
    Q_values: Union[List[float], np.ndarray, None] = None,
    head_gradient: Union[List[float], np.ndarray, None] = None,
    time_values: Union[List[float], np.ndarray, None] = None,
    head_values: Union[List[float], np.ndarray, None] = None,
    test_type: str = "constant",
    ax=None,
    save_as: Optional[str] = None,
) -> Dict:
    r"""
    Plot permeability test data.

    For constant-head: plots Q vs hydraulic gradient i.
    For falling-head: plots ln(h) vs time.

    Parameters
    ----------
    Q_values : array_like, optional
        Flow volumes (cm³) — for constant-head.
    head_gradient : array_like, optional
        Hydraulic gradients — for constant-head.
    time_values : array_like, optional
        Time readings (s) — for falling-head.
    head_values : array_like, optional
        Head readings (cm) — for falling-head.
    test_type : {'constant', 'falling'}, default ``'constant'``
    ax : matplotlib.axes.Axes, optional
    save_as : str, optional

    Returns
    -------
    dict
        ``'ax'`` : Axes, ``'k'`` : float (slope-derived k estimate).

    Examples
    --------
    >>> from geoeq.lab.permeability import permeability_plot
    >>> res = permeability_plot(Q_values=[10, 20, 30],
    ...     head_gradient=[1.0, 2.0, 3.0], test_type='constant')
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.get_figure()

    k_est = None

    if test_type == "constant":
        if Q_values is None or head_gradient is None:
            raise ValueError("Q_values and head_gradient required for constant-head plot.")
        Q = np.asarray(Q_values, dtype=float)
        i = np.asarray(head_gradient, dtype=float)
        ax.scatter(i, Q, s=80, color="steelblue", edgecolors="navy", zorder=5)

        coeffs = np.polyfit(i, Q, 1)
        i_fit = np.linspace(0, float(i.max()) * 1.2, 100)
        ax.plot(i_fit, np.polyval(coeffs, i_fit), "r-", linewidth=2)

        k_est = float(coeffs[0])
        ax.set_xlabel("Hydraulic Gradient i", fontweight="bold")
        ax.set_ylabel("Flow Rate Q (cm³/s)", fontweight="bold")
        ax.set_title("Constant-Head Permeability Test", fontweight="bold")

    elif test_type == "falling":
        if time_values is None or head_values is None:
            raise ValueError("time_values and head_values required for falling-head plot.")
        t = np.asarray(time_values, dtype=float)
        h = np.asarray(head_values, dtype=float)
        ln_h = np.log(h)

        ax.scatter(t, ln_h, s=80, color="steelblue", edgecolors="navy", zorder=5)
        coeffs = np.polyfit(t, ln_h, 1)
        t_fit = np.linspace(0, float(t.max()) * 1.1, 100)
        ax.plot(t_fit, np.polyval(coeffs, t_fit), "r-", linewidth=2)

        ax.set_xlabel("Time (s)", fontweight="bold")
        ax.set_ylabel("ln(h)", fontweight="bold")
        ax.set_title("Falling-Head Permeability Test", fontweight="bold")
        k_est = float(-coeffs[0])

    else:
        raise ValueError(f"test_type must be 'constant' or 'falling', got '{test_type}'.")

    ax.grid(True, alpha=0.3)

    if save_as:
        plt.savefig(save_as, bbox_inches="tight", dpi=300)

    return {"ax": ax, "k": k_est}
