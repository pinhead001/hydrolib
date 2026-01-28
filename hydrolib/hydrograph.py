"""
hydrolib.hydrograph - Hydrograph analysis and plotting
"""

from __future__ import annotations

from functools import lru_cache
from typing import ClassVar, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.dates import DateFormatter


class Hydrograph:
    """Class for hydrograph analysis and plotting."""

    MONTH_STARTS: ClassVar[List[int]] = [1, 32, 62, 93, 123, 154, 184, 215, 246, 274, 305, 335]
    MONTH_LABELS: ClassVar[List[str]] = [
        "Oct",
        "Nov",
        "Dec",
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
    ]

    @staticmethod
    @lru_cache(maxsize=1024)
    def day_of_water_year(month: int, day: int) -> int:
        """Convert month/day to day of water year (1-366, starting Oct 1)."""
        days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        wy_month = (month - 10) % 12
        dowy = sum(days_in_month[(10 + i - 1) % 12] for i in range(wy_month)) + day
        return dowy

    @classmethod
    def _compute_dowy_series(cls, dates: pd.DatetimeIndex) -> np.ndarray:
        """Compute day of water year for a date series."""
        return np.array([cls.day_of_water_year(d.month, d.day) for d in dates])

    @classmethod
    def plot_daily_timeseries(
        cls,
        daily_data: pd.DataFrame,
        site_name: str = None,
        site_no: str = None,
        save_path: str = None,
        figsize: Tuple[int, int] = (12, 5),
        color: str = "steelblue",
    ) -> plt.Figure:
        """Plot daily flow time series."""
        fig, ax = plt.subplots(figsize=figsize)

        ax.plot(daily_data.index, daily_data["flow_cfs"], color=color, linewidth=0.5, alpha=0.8)
        ax.set_yscale("log")
        ax.set_ylabel("Discharge (cfs)", fontsize=11)
        ax.set_xlabel("Date", fontsize=11)

        title = "Mean Daily Streamflow"
        if site_name and site_no:
            title = f"Mean Daily Streamflow\nUSGS {site_no} - {site_name}"
        elif site_no:
            title = f"Mean Daily Streamflow - USGS {site_no}"
        ax.set_title(title, fontsize=12, fontweight="bold")

        ax.grid(True, which="both", alpha=0.3)
        ax.xaxis.set_major_formatter(DateFormatter("%Y"))

        start_yr = daily_data.index.min().year
        end_yr = daily_data.index.max().year
        ax.annotate(
            f"Period of Record: {start_yr}-{end_yr}",
            xy=(0.02, 0.98),
            xycoords="axes fraction",
            fontsize=9,
            ha="left",
            va="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches="tight")

        return fig

    @classmethod
    def plot_summary_hydrograph(
        cls,
        daily_data: pd.DataFrame,
        site_name: str = None,
        site_no: str = None,
        save_path: str = None,
        figsize: Tuple[int, int] = (10, 6),
        percentiles: List[int] = None,
    ) -> plt.Figure:
        """Plot summary hydrograph with day of water year on x-axis."""
        if percentiles is None:
            percentiles = [10, 25, 50, 75, 90]

        df = daily_data.copy()
        df["dowy"] = cls._compute_dowy_series(df.index)

        stats_df = df.groupby("dowy")["flow_cfs"].agg(
            ["mean", "min", "max"] + [lambda x, p=p: np.nanpercentile(x, p) for p in percentiles]
        )
        stats_df.columns = ["mean", "min", "max"] + [f"p{p}" for p in percentiles]

        fig, ax = plt.subplots(figsize=figsize)

        if len(percentiles) >= 4:
            ax.fill_between(
                stats_df.index,
                stats_df[f"p{percentiles[0]}"],
                stats_df[f"p{percentiles[-1]}"],
                alpha=0.3,
                color="blue",
                label=f"{percentiles[0]}-{percentiles[-1]}th percentile",
            )
            ax.fill_between(
                stats_df.index,
                stats_df[f"p{percentiles[1]}"],
                stats_df[f"p{percentiles[-2]}"],
                alpha=0.4,
                color="blue",
                label=f"{percentiles[1]}-{percentiles[-2]}th percentile",
            )

        ax.plot(stats_df.index, stats_df["p50"], "b-", linewidth=2, label="Median")
        ax.plot(stats_df.index, stats_df["mean"], "k--", linewidth=1.5, label="Mean")

        ax.set_yscale("log")
        ax.set_ylabel("Discharge (cfs)", fontsize=11)
        ax.set_xlabel("Day of Water Year", fontsize=11)

        ax.set_xticks(cls.MONTH_STARTS)
        ax.set_xticklabels(cls.MONTH_LABELS)
        ax.set_xlim(1, 366)

        title = "Summary Hydrograph"
        if site_name and site_no:
            title = f"Summary Hydrograph\nUSGS {site_no} - {site_name}"
        elif site_no:
            title = f"Summary Hydrograph - USGS {site_no}"
        ax.set_title(title, fontsize=12, fontweight="bold")

        ax.legend(loc="upper right", fontsize=9)
        ax.grid(True, which="both", alpha=0.3)

        start_yr = daily_data.index.min().year
        end_yr = daily_data.index.max().year
        n_years = len(df.index.year.unique())
        ax.annotate(
            f"Period of Record: {start_yr}-{end_yr} ({n_years} years)",
            xy=(0.02, 0.02),
            xycoords="axes fraction",
            fontsize=9,
            ha="left",
            va="bottom",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches="tight")

        return fig
