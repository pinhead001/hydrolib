"""
hydrolib.freq_plot - Frequency curve plotting for Streamlit display.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

if TYPE_CHECKING:
    from .bulletin17c import Bulletin17C

logger = logging.getLogger(__name__)


def plot_frequency_curve_streamlit(
    b17c: "Bulletin17C",
    site_name: Optional[str] = None,
    site_no: Optional[str] = None,
    figsize: Tuple[int, int] = (10, 6),
) -> plt.Figure:
    """Plot a Bulletin 17C frequency curve suitable for Streamlit display.

    Parameters
    ----------
    b17c : Bulletin17C
        A Bulletin17C instance (analysis will be run if not already done).
    site_name : str, optional
        Station name for the title.
    site_no : str, optional
        USGS station number for the title.
    figsize : tuple
        Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure
        The frequency curve figure.
    """
    if b17c.results is None:
        b17c.run_analysis()

    r = b17c.results

    def aep_to_x(p: float) -> float:
        return norm.ppf(1 - p)

    fig, ax = plt.subplots(figsize=figsize, facecolor="white")
    ax.set_facecolor("white")

    # --- Observed peaks (Weibull plotting positions) ---
    peak_flows = b17c._peak_flows  # noqa: SLF001
    sorted_flows = np.sort(peak_flows)[::-1]
    n = len(sorted_flows)
    weibull_aep = np.arange(1, n + 1) / (n + 1)
    x_obs = [aep_to_x(p) for p in weibull_aep]

    ax.scatter(
        x_obs,
        sorted_flows,
        c="blue",
        s=40,
        zorder=5,
        label="Observed Annual Peaks",
        edgecolors="darkblue",
        linewidth=0.5,
    )

    # --- Fitted LP3 curve at fine spacing ---
    aep_fine = np.linspace(0.001, 0.999, 200)
    quantiles_df = b17c.compute_quantiles(aep_fine)
    x_curve = [aep_to_x(p) for p in aep_fine]

    ax.plot(
        x_curve,
        quantiles_df["flow_cfs"].values,
        "b-",
        linewidth=2,
        label="LP3 Fitted Curve",
        zorder=4,
    )

    # --- 90% confidence interval ---
    aep_ci = np.array(
        [0.999, 0.995, 0.99, 0.95, 0.90, 0.80, 0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002]
    )
    cl = b17c.compute_confidence_limits(aep_ci, confidence=0.90)
    x_cl = [aep_to_x(p) for p in cl["aep"]]
    ax.fill_between(
        x_cl,
        cl["lower_5pct"],
        cl["upper_5pct"],
        alpha=0.2,
        color="blue",
        label="90% Confidence Interval",
    )
    ax.plot(x_cl, cl["lower_5pct"], "b--", linewidth=0.8, alpha=0.6)
    ax.plot(x_cl, cl["upper_5pct"], "b--", linewidth=0.8, alpha=0.6)

    # --- Low outlier threshold ---
    if r.low_outlier_threshold > 0 and r.n_low_outliers > 0:
        ax.axhline(
            r.low_outlier_threshold,
            color="red",
            linestyle="--",
            alpha=0.7,
            label=f"MGBT Threshold ({r.n_low_outliers} low outliers)",
        )

    # --- Axes ---
    ax.set_yscale("log")
    ax.set_ylabel("Peak Discharge (cfs)", fontsize=11)
    ax.set_xlabel("Annual Exceedance Probability", fontsize=11)

    # X-axis probability ticks
    prob_ticks = [0.999, 0.99, 0.95, 0.90, 0.80, 0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002]
    x_ticks = [aep_to_x(p) for p in prob_ticks]
    ax.set_xticks(x_ticks)
    tick_labels = []
    for p in prob_ticks:
        pct = p * 100
        if pct >= 1:
            tick_labels.append(f"{pct:g}%")
        else:
            tick_labels.append(f"{pct:.1f}%")
    ax.set_xticklabels(tick_labels, fontsize=9, rotation=45, ha="right")
    ax.set_xlim(aep_to_x(0.999), aep_to_x(0.002))

    # --- Secondary axis: return periods ---
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    rp_ticks = [1.5, 2, 5, 10, 25, 50, 100, 200, 500]
    rp_x = [aep_to_x(1 / rp) for rp in rp_ticks]
    ax2.set_xticks(rp_x)
    ax2.set_xticklabels([str(int(rp)) if rp == int(rp) else str(rp) for rp in rp_ticks], fontsize=9)
    ax2.set_xlabel("Return Period (years)", fontsize=11)

    # --- Title ---
    if site_name and site_no:
        title = f"{site_name}\nUSGS {site_no}"
    elif site_no:
        title = f"USGS {site_no}"
    else:
        title = "Flood Frequency Curve (Bulletin 17C)"
    ax.set_title(title, fontsize=12, fontweight="bold", pad=35)

    # --- LP3 parameters text box ---
    stats_text = (
        f"LP3 Parameters: "
        f"\u03bc={r.mean_log:.3f}  "
        f"\u03c3={r.std_log:.3f}  "
        f"\u03b3={r.skew_used:.3f}"
    )
    ax.text(
        0.02,
        0.02,
        stats_text,
        transform=ax.transAxes,
        fontsize=9,
        va="bottom",
        ha="left",
        family="monospace",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.9, edgecolor="gray"),
    )

    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()

    return fig
