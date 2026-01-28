"""
hydrolib.bulletin17c - Bulletin 17C Flood Frequency Analysis

Implements both Method of Moments (MOM) and Expected Moments Algorithm (EMA).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from datetime import datetime
from functools import cached_property
from typing import ClassVar, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.integrate import quad
from scipy.special import gammaln, ndtri

from .core import (
    AnalysisMethod,
    EMAParameters,
    FlowInterval,
    FrequencyResults,
    grubbs_beck_critical_value,
    kfactor,
    kfactor_array,
    log_pearson3_cdf,
    log_pearson3_pdf,
)


class FloodFrequencyAnalysis(ABC):
    """Abstract base class for flood frequency analysis methods."""

    STANDARD_AEP: ClassVar[np.ndarray] = np.array(
        [0.995, 0.99, 0.95, 0.90, 0.80, 0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002]
    )

    def __init__(
        self, peak_flows: np.ndarray, regional_skew: float = None, regional_skew_mse: float = None
    ):
        self._peak_flows = np.array(peak_flows)
        self._peak_flows = self._peak_flows[~np.isnan(self._peak_flows)]
        self._peak_flows = self._peak_flows[self._peak_flows > 0]

        self._regional_skew = regional_skew
        self._regional_skew_mse = regional_skew_mse
        self._results: Optional[FrequencyResults] = None

    @property
    def n(self) -> int:
        return len(self._peak_flows)

    @property
    def peak_flows(self) -> np.ndarray:
        return self._peak_flows.copy()

    @cached_property
    def log_flows(self) -> np.ndarray:
        return np.log10(self._peak_flows)

    @property
    def results(self) -> Optional[FrequencyResults]:
        return self._results

    @abstractmethod
    def run_analysis(self) -> FrequencyResults:
        pass

    def _compute_weighted_skew(self, station_skew: float) -> Optional[float]:
        """Compute weighted skew from station and regional values."""
        if self._regional_skew is None or self._regional_skew_mse is None:
            return None

        n = self.n
        mse_station = (6 / n) * (1 + (6 / n) * station_skew**2 + (15 / (n**2)) * station_skew**4)

        w_regional = mse_station / (mse_station + self._regional_skew_mse)
        w_station = self._regional_skew_mse / (mse_station + self._regional_skew_mse)

        return w_station * station_skew + w_regional * self._regional_skew

    def compute_quantiles(self, aep: np.ndarray = None) -> pd.DataFrame:
        """Compute flood frequency quantiles."""
        if self._results is None:
            self.run_analysis()

        if aep is None:
            aep = self.STANDARD_AEP

        K = kfactor_array(self._results.skew_used, aep)
        log_Q = self._results.mean_log + K * self._results.std_log
        Q = 10**log_Q

        return pd.DataFrame(
            {"aep": aep, "return_period": 1 / aep, "flow_cfs": Q, "log_flow": log_Q, "K_factor": K}
        )

    def compute_confidence_limits(
        self, aep: np.ndarray = None, confidence: float = 0.90
    ) -> pd.DataFrame:
        """Compute confidence limits for quantile estimates."""
        if self._results is None:
            self.run_analysis()

        if aep is None:
            aep = self.STANDARD_AEP

        alpha = 1 - confidence
        z_alpha = ndtri(1 - alpha / 2)

        n = self._results.n_systematic or self.n
        K = kfactor_array(self._results.skew_used, aep)
        G = self._results.skew_used

        var_factor = 1 / n + K**2 * (1 + 0.75 * G**2) / (2 * (n - 1))
        se_log = self._results.std_log * np.sqrt(var_factor)

        log_Q = self._results.mean_log + K * self._results.std_log
        log_lower = log_Q - z_alpha * se_log
        log_upper = log_Q + z_alpha * se_log

        return pd.DataFrame(
            {
                "aep": aep,
                "return_period": 1 / aep,
                "flow_cfs": 10**log_Q,
                "lower_5pct": 10**log_lower,
                "upper_5pct": 10**log_upper,
            }
        )

    def plot_frequency_curve(
        self,
        site_name: str = None,
        site_no: str = None,
        save_path: str = None,
        figsize: Tuple[int, int] = (10, 7),
        show_confidence: bool = True,
    ) -> plt.Figure:
        """Plot professional flood frequency curve."""
        if self._results is None:
            self.run_analysis()

        fig, ax = plt.subplots(figsize=figsize)

        def prob_to_x(p):
            return ndtri(1 - p)

        # Plot observed data
        sorted_flows = np.sort(self._peak_flows)[::-1]
        n = len(sorted_flows)
        plotting_position = np.arange(1, n + 1) / (n + 1)

        x_data = [prob_to_x(p) for p in plotting_position]
        ax.scatter(
            x_data,
            sorted_flows,
            c="blue",
            s=40,
            zorder=5,
            label="Observed Annual Peaks",
            edgecolors="darkblue",
            linewidth=0.5,
        )

        # Plot fitted curve
        aep_curve = np.array(
            [
                0.999,
                0.995,
                0.99,
                0.95,
                0.90,
                0.80,
                0.50,
                0.20,
                0.10,
                0.04,
                0.02,
                0.01,
                0.005,
                0.002,
                0.001,
            ]
        )
        K_curve = kfactor_array(self._results.skew_used, aep_curve)
        Q_curve = 10 ** (self._results.mean_log + K_curve * self._results.std_log)
        x_curve = [prob_to_x(p) for p in aep_curve]

        method_label = "EMA" if self._results.method == AnalysisMethod.EMA else "MOM"
        ax.plot(
            x_curve,
            Q_curve,
            "b-",
            linewidth=2,
            label=f"LP3 Fitted Curve ({method_label})",
            zorder=4,
        )

        # Confidence limits
        if show_confidence:
            cl = self.compute_confidence_limits(aep_curve)
            x_cl = [prob_to_x(p) for p in cl["aep"]]
            ax.fill_between(
                x_cl,
                cl["lower_5pct"],
                cl["upper_5pct"],
                alpha=0.2,
                color="blue",
                label="90% Confidence Interval",
            )
            ax.plot(x_cl, cl["lower_5pct"], "b--", linewidth=1, alpha=0.7)
            ax.plot(x_cl, cl["upper_5pct"], "b--", linewidth=1, alpha=0.7)

        # Low outliers
        if self._results.low_outlier_threshold > 0 and self._results.n_low_outliers > 0:
            threshold = self._results.low_outlier_threshold
            low_outliers = self._peak_flows[self._peak_flows < threshold]
            for lo in low_outliers:
                idx = np.where(sorted_flows == lo)[0]
                if len(idx) > 0:
                    ax.scatter(
                        [x_data[i] for i in idx],
                        [lo] * len(idx),
                        c="red",
                        s=60,
                        marker="x",
                        zorder=6,
                    )
            ax.axhline(
                threshold,
                color="red",
                linestyle=":",
                alpha=0.7,
                label=f"Low Outlier Threshold ({self._results.n_low_outliers})",
            )

        ax.set_yscale("log")
        ax.set_ylabel("Peak Discharge (cfs)", fontsize=11)
        ax.set_xlabel("Annual Exceedance Probability", fontsize=11)

        # X-axis ticks
        prob_ticks = [0.99, 0.95, 0.90, 0.80, 0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002]
        x_ticks = [prob_to_x(p) for p in prob_ticks]
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(
            [f"{p*100:.1f}%" if p >= 0.01 else f"{p*100:.2f}%" for p in prob_ticks],
            fontsize=9,
            rotation=45,
            ha="right",
        )

        # Secondary axis for return period
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        rp_ticks = [1.01, 2, 5, 10, 25, 50, 100, 200, 500]
        rp_x = [prob_to_x(1 / rp) for rp in rp_ticks]
        ax2.set_xticks(rp_x)
        ax2.set_xticklabels(
            [str(int(rp)) if rp >= 1 else f"{rp:.2f}" for rp in rp_ticks], fontsize=9
        )
        ax2.set_xlabel("Return Period (years)", fontsize=11)

        # Title
        title = f"Flood Frequency Analysis (Bulletin 17C - {method_label})"
        if site_name and site_no:
            title = f"Flood Frequency Analysis (Bulletin 17C - {method_label})\nUSGS {site_no} - {site_name}"
        ax.set_title(title, fontsize=12, fontweight="bold", pad=35)

        # Stats annotation
        r = self._results
        stats_text = f"n = {r.n_peaks}"
        if r.n_historical > 0:
            stats_text += f" ({r.n_systematic} sys + {r.n_historical} hist)"
        stats_text += f"\nMean(log Q) = {r.mean_log:.4f}"
        stats_text += f"\nStd(log Q) = {r.std_log:.4f}"
        stats_text += f"\nStation Skew = {r.skew_station:.3f}"
        if r.skew_weighted is not None:
            stats_text += f"\nWeighted Skew = {r.skew_weighted:.3f}"

        ax.annotate(
            stats_text,
            xy=(0.02, 0.98),
            xycoords="axes fraction",
            fontsize=9,
            ha="left",
            va="top",
            family="monospace",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.9),
        )

        ax.legend(loc="lower right", fontsize=9)
        ax.grid(True, which="both", alpha=0.3)
        ax.set_xlim(prob_to_x(0.999), prob_to_x(0.001))

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches="tight")

        return fig


class MethodOfMoments(FloodFrequencyAnalysis):
    """Traditional Method of Moments flood frequency analysis."""

    def run_analysis(self) -> FrequencyResults:
        """Run Method of Moments analysis."""
        log_flows = self.log_flows
        n = self.n

        mean_log = np.mean(log_flows)
        std_log = np.std(log_flows, ddof=1)

        skew_station = n * np.sum((log_flows - mean_log) ** 3) / ((n - 1) * (n - 2) * std_log**3)

        skew_weighted = self._compute_weighted_skew(skew_station)
        skew_used = skew_weighted if skew_weighted is not None else skew_station

        k_n = grubbs_beck_critical_value(n)
        threshold_log = mean_log - k_n * std_log
        low_outlier_threshold = 10**threshold_log
        n_low_outliers = int(np.sum(self._peak_flows < low_outlier_threshold))

        self._results = FrequencyResults(
            n_peaks=n,
            n_systematic=n,
            n_historical=0,
            n_censored=0,
            n_low_outliers=n_low_outliers,
            mean_log=mean_log,
            std_log=std_log,
            skew_station=skew_station,
            skew_regional=self._regional_skew,
            skew_weighted=skew_weighted,
            skew_used=skew_used,
            low_outlier_threshold=low_outlier_threshold,
            mgb_critical_value=k_n,
            method=AnalysisMethod.MOM,
        )

        self._results.quantiles = self.compute_quantiles()
        self._results.confidence_limits = self.compute_confidence_limits()

        return self._results


class ExpectedMomentsAlgorithm(FloodFrequencyAnalysis):
    """
    Bulletin 17C Expected Moments Algorithm (EMA) implementation.

    Handles systematic record, historical floods, censored observations,
    and low outliers identified by Multiple Grubbs-Beck test.
    """

    def __init__(
        self,
        peak_flows: np.ndarray,
        water_years: np.ndarray = None,
        regional_skew: float = None,
        regional_skew_mse: float = None,
        ema_params: EMAParameters = None,
        historical_peaks: List[Tuple[int, float]] = None,
        perception_thresholds: Dict[Tuple[int, int], float] = None,
    ):
        super().__init__(peak_flows, regional_skew, regional_skew_mse)

        if water_years is not None:
            self._water_years = np.array(water_years)
        else:
            end_year = datetime.now().year
            self._water_years = np.arange(end_year - len(peak_flows) + 1, end_year + 1)

        self._historical_peaks = historical_peaks or []
        self._perception_thresholds = perception_thresholds or {}
        self._ema_params = ema_params or self._auto_configure_ema_params()
        self._intervals: List[FlowInterval] = []

    def _auto_configure_ema_params(self) -> EMAParameters:
        """Auto-configure EMA parameters from data."""
        sys_start = int(self._water_years.min())
        sys_end = int(self._water_years.max())

        all_years = set(range(sys_start, sys_end + 1))
        recorded_years = set(self._water_years.astype(int))
        gaps = sorted(all_years - recorded_years)

        hist_start = None
        hist_end = None
        hist_threshold = None

        if self._historical_peaks:
            hist_years = [h[0] for h in self._historical_peaks]
            hist_start = min(hist_years)
            hist_end = max(max(hist_years), sys_start - 1)
            hist_threshold = max(h[1] for h in self._historical_peaks)
        elif gaps:
            first_gap = gaps[0]
            pre_gap_mask = self._water_years < first_gap
            if np.any(pre_gap_mask):
                pre_gap_flows = self._peak_flows[pre_gap_mask]
                hist_threshold = float(np.max(pre_gap_flows))
                hist_start = sys_start
                hist_end = first_gap - 1

        return EMAParameters(
            systematic_start=sys_start,
            systematic_end=sys_end,
            historical_start=hist_start,
            historical_end=hist_end,
            historical_threshold=hist_threshold,
        )

    def _build_flow_intervals(self, low_threshold: float = 0.0) -> List[FlowInterval]:
        """Build flow intervals for EMA analysis."""
        intervals = []

        for flow, year in zip(self._peak_flows, self._water_years):
            year = int(year)

            if flow < low_threshold:
                intervals.append(
                    FlowInterval.from_censored(
                        lower=0, upper=low_threshold, year=year, perception_threshold=0.0
                    )
                )
            else:
                intervals.append(FlowInterval.from_peak(flow, year))

        for year, flow in self._historical_peaks:
            threshold = self._ema_params.historical_threshold or flow
            intervals.append(
                FlowInterval.from_historical(flow, year, perception_threshold=threshold)
            )

        if self._ema_params.historical_start and self._ema_params.historical_threshold:
            hist_recorded_years = {h[0] for h in self._historical_peaks}
            sys_years = set(self._water_years.astype(int))

            for year in range(
                self._ema_params.historical_start, self._ema_params.historical_end + 1
            ):
                if year not in hist_recorded_years and year not in sys_years:
                    intervals.append(
                        FlowInterval.from_censored(
                            lower=0,
                            upper=self._ema_params.historical_threshold,
                            year=year,
                            perception_threshold=self._ema_params.historical_threshold,
                        )
                    )

        self._intervals = sorted(intervals, key=lambda x: x.year)
        return self._intervals

    def _compute_ema_moments(
        self, mean_log: float, std_log: float, skew: float
    ) -> Tuple[float, float, float]:
        """Compute expected moments given current parameter estimates."""
        sum_x = 0.0
        sum_x2 = 0.0
        sum_x3 = 0.0

        for interval in self._intervals:
            if not interval.is_censored:
                x = np.log10(interval.lower)
                sum_x += x
                sum_x2 += x**2
                sum_x3 += x**3
            else:
                lower = interval.lower if interval.lower > 0 else 1e-10
                upper = interval.upper if np.isfinite(interval.upper) else 1e10

                p_lower = log_pearson3_cdf(lower, mean_log, std_log, skew)
                p_upper = log_pearson3_cdf(upper, mean_log, std_log, skew)
                p_interval = p_upper - p_lower

                if p_interval < 1e-10:
                    x = (np.log10(lower) + np.log10(upper)) / 2
                    sum_x += x
                    sum_x2 += x**2
                    sum_x3 += x**3
                else:
                    log_lower = np.log10(lower)
                    log_upper = np.log10(upper)

                    def pdf(log_q):
                        return log_pearson3_pdf(log_q, mean_log, std_log, skew)

                    try:
                        ex, _ = quad(lambda x: x * pdf(x), log_lower, log_upper)
                        ex /= p_interval
                        ex2, _ = quad(lambda x: x**2 * pdf(x), log_lower, log_upper)
                        ex2 /= p_interval
                        ex3, _ = quad(lambda x: x**3 * pdf(x), log_lower, log_upper)
                        ex3 /= p_interval
                    except:
                        x = (log_lower + log_upper) / 2
                        ex, ex2, ex3 = x, x**2, x**3

                    sum_x += ex
                    sum_x2 += ex2
                    sum_x3 += ex3

        return sum_x, sum_x2, sum_x3

    def _ema_iteration(
        self, mean_log: float, std_log: float, skew: float
    ) -> Tuple[float, float, float]:
        """Perform one EMA iteration."""
        n = len(self._intervals)

        sum_x, sum_x2, sum_x3 = self._compute_ema_moments(mean_log, std_log, skew)

        new_mean = sum_x / n
        var = (sum_x2 - n * new_mean**2) / (n - 1)
        new_std = np.sqrt(max(var, 1e-10))

        if n > 2:
            m3 = sum_x3 - 3 * new_mean * sum_x2 + 2 * n * new_mean**3
            new_skew = (n * m3) / ((n - 1) * (n - 2) * new_std**3)
        else:
            new_skew = skew

        return new_mean, new_std, new_skew

    def _multiple_grubbs_beck(self) -> Tuple[float, int]:
        """Perform Multiple Grubbs-Beck test for low outliers."""
        log_flows = self.log_flows
        current_flows = np.sort(log_flows)
        n_outliers = 0
        threshold_log = -np.inf

        while len(current_flows) > 10:
            mean_log = np.mean(current_flows)
            std_log = np.std(current_flows, ddof=1)

            k_n = grubbs_beck_critical_value(len(current_flows))
            test_threshold = mean_log - k_n * std_log

            if current_flows[0] < test_threshold:
                threshold_log = test_threshold
                n_outliers += 1
                current_flows = current_flows[1:]
            else:
                break

        threshold = 10**threshold_log if threshold_log > -np.inf else 0.0
        return threshold, n_outliers

    def run_analysis(self) -> FrequencyResults:
        """Run EMA flood frequency analysis."""
        log_flows = self.log_flows
        n = self.n

        mean_log = np.mean(log_flows)
        std_log = np.std(log_flows, ddof=1)
        skew_station = n * np.sum((log_flows - mean_log) ** 3) / ((n - 1) * (n - 2) * std_log**3)

        low_threshold, n_low_outliers = self._multiple_grubbs_beck()
        self._ema_params.low_outlier_threshold = low_threshold

        self._build_flow_intervals(low_threshold)

        converged = False
        iteration = 0

        for iteration in range(self._ema_params.max_iterations):
            new_mean, new_std, new_skew = self._ema_iteration(mean_log, std_log, skew_station)

            if (
                abs(new_mean - mean_log) < self._ema_params.tolerance
                and abs(new_std - std_log) < self._ema_params.tolerance
                and abs(new_skew - skew_station) < self._ema_params.tolerance
            ):
                converged = True
                mean_log = new_mean
                std_log = new_std
                skew_station = new_skew
                break

            mean_log = new_mean
            std_log = new_std
            skew_station = new_skew

        skew_weighted = self._compute_weighted_skew(skew_station)
        skew_used = skew_weighted if skew_weighted is not None else skew_station

        n_systematic = sum(1 for i in self._intervals if i.is_systematic and not i.is_censored)
        n_historical = sum(1 for i in self._intervals if i.is_historical)
        n_censored = sum(1 for i in self._intervals if i.is_censored)

        self._results = FrequencyResults(
            n_peaks=len(self._intervals),
            n_systematic=n_systematic,
            n_historical=n_historical,
            n_censored=n_censored,
            n_low_outliers=n_low_outliers,
            mean_log=mean_log,
            std_log=std_log,
            skew_station=skew_station,
            skew_regional=self._regional_skew,
            skew_weighted=skew_weighted,
            skew_used=skew_used,
            low_outlier_threshold=low_threshold,
            mgb_critical_value=grubbs_beck_critical_value(n),
            method=AnalysisMethod.EMA,
            ema_iterations=iteration + 1,
            ema_converged=converged,
        )

        self._results.quantiles = self.compute_quantiles()
        self._results.confidence_limits = self.compute_confidence_limits()

        return self._results

    @property
    def intervals(self) -> List[FlowInterval]:
        return self._intervals.copy()


class Bulletin17C:
    """
    Unified interface for Bulletin 17C flood frequency analysis.

    Allows selection between Method of Moments (MOM) and Expected
    Moments Algorithm (EMA) methods.

    Examples
    --------
    >>> b17c = Bulletin17C(peak_flows)
    >>> results = b17c.run_analysis(method='mom')

    >>> b17c = Bulletin17C(peak_flows, regional_skew=-0.05, regional_skew_mse=0.12)
    >>> results = b17c.run_analysis(method='ema')
    """

    def __init__(
        self,
        peak_flows: np.ndarray,
        water_years: np.ndarray = None,
        regional_skew: float = None,
        regional_skew_mse: float = None,
        historical_peaks: List[Tuple[int, float]] = None,
        perception_thresholds: Dict[Tuple[int, int], float] = None,
        ema_params: EMAParameters = None,
    ):
        self._peak_flows = np.array(peak_flows)
        self._water_years = water_years
        self._regional_skew = regional_skew
        self._regional_skew_mse = regional_skew_mse
        self._historical_peaks = historical_peaks
        self._perception_thresholds = perception_thresholds
        self._ema_params = ema_params

        self._analyzer: Optional[FloodFrequencyAnalysis] = None
        self._results: Optional[FrequencyResults] = None

    @property
    def results(self) -> Optional[FrequencyResults]:
        return self._results

    @property
    def quantiles(self) -> Optional[pd.DataFrame]:
        return self._results.quantiles if self._results else None

    @property
    def confidence_limits(self) -> Optional[pd.DataFrame]:
        return self._results.confidence_limits if self._results else None

    @property
    def n(self) -> int:
        return len(self._peak_flows)

    @property
    def mean_log(self) -> Optional[float]:
        return self._results.mean_log if self._results else None

    @property
    def std_log(self) -> Optional[float]:
        return self._results.std_log if self._results else None

    @property
    def skew_station(self) -> Optional[float]:
        return self._results.skew_station if self._results else None

    @property
    def skew_weighted(self) -> Optional[float]:
        return self._results.skew_weighted if self._results else None

    @property
    def skew_used(self) -> Optional[float]:
        return self._results.skew_used if self._results else None

    @property
    def low_outlier_threshold(self) -> Optional[float]:
        return self._results.low_outlier_threshold if self._results else None

    @property
    def n_low_outliers(self) -> Optional[int]:
        return self._results.n_low_outliers if self._results else None

    def run_analysis(self, method: Union[str, AnalysisMethod] = "ema") -> FrequencyResults:
        """
        Run flood frequency analysis.

        Parameters
        ----------
        method : str or AnalysisMethod
            'mom' for Method of Moments, 'ema' for Expected Moments Algorithm
        """
        if isinstance(method, str):
            method = AnalysisMethod[method.upper()]

        if method == AnalysisMethod.MOM:
            self._analyzer = MethodOfMoments(
                self._peak_flows,
                regional_skew=self._regional_skew,
                regional_skew_mse=self._regional_skew_mse,
            )
        else:
            self._analyzer = ExpectedMomentsAlgorithm(
                self._peak_flows,
                water_years=self._water_years,
                regional_skew=self._regional_skew,
                regional_skew_mse=self._regional_skew_mse,
                ema_params=self._ema_params,
                historical_peaks=self._historical_peaks,
                perception_thresholds=self._perception_thresholds,
            )

        self._results = self._analyzer.run_analysis()
        return self._results

    def compute_quantiles(self, aep: np.ndarray = None) -> pd.DataFrame:
        if self._analyzer is None:
            self.run_analysis()
        return self._analyzer.compute_quantiles(aep)

    def compute_confidence_limits(
        self, aep: np.ndarray = None, confidence: float = 0.90
    ) -> pd.DataFrame:
        if self._analyzer is None:
            self.run_analysis()
        return self._analyzer.compute_confidence_limits(aep, confidence)

    def plot_frequency_curve(
        self,
        site_name: str = None,
        site_no: str = None,
        save_path: str = None,
        figsize: Tuple[int, int] = (10, 7),
        show_confidence: bool = True,
    ) -> plt.Figure:
        if self._analyzer is None:
            self.run_analysis()
        return self._analyzer.plot_frequency_curve(
            site_name, site_no, save_path, figsize, show_confidence
        )
