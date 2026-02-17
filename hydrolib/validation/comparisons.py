"""
Comparison engine for native vs reference frequency analysis results.

Compares HydroLib native EMA output against PeakfqSA/reference results
with configurable tolerance thresholds per output category.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

from hydrolib.peakfqsa.parsers import PeakfqSAResult

logger = logging.getLogger(__name__)


@dataclass
class ComparisonResult:
    """Result of comparing native and reference frequency analyses.

    Parameters
    ----------
    passed : bool
        True if all differences are within tolerance.
    tolerance_pct : float
        Tolerance threshold used for the comparison.
    parameter_diffs : dict[str, float]
        Parameter name to percent difference mapping.
    quantile_diffs : dict[float, float]
        AEP to percent difference mapping for quantiles.
    ci_diffs : dict[float, float]
        AEP to percent difference mapping for CI bounds.
    max_diff_pct : float
        Maximum percent difference across all comparisons.
    summary : str
        Human-readable one-line summary of comparison.
    """

    passed: bool = False
    tolerance_pct: float = 1.0
    parameter_diffs: dict[str, float] = field(default_factory=dict)
    quantile_diffs: dict[float, float] = field(default_factory=dict)
    ci_diffs: dict[float, float] = field(default_factory=dict)
    max_diff_pct: float = 0.0
    summary: str = ""


class FrequencyComparator:
    """Compare native HydroLib analysis against a reference result.

    Parameters
    ----------
    tolerance_pct : float
        Default tolerance for quantile comparisons.
    parameter_tolerance_pct : float
        Tolerance for LP3 parameter comparisons.
    ci_tolerance_pct : float
        Tolerance for confidence interval comparisons.
    """

    def __init__(
        self,
        tolerance_pct: float = 1.0,
        parameter_tolerance_pct: float = 0.5,
        ci_tolerance_pct: float = 2.0,
    ) -> None:
        self.tolerance_pct = tolerance_pct
        self.parameter_tolerance_pct = parameter_tolerance_pct
        self.ci_tolerance_pct = ci_tolerance_pct

    def compare(
        self,
        native: dict[str, Any],
        reference: PeakfqSAResult,
    ) -> ComparisonResult:
        """Compare native analysis output against a reference result.

        Parameters
        ----------
        native : dict
            Output from HydroLib native analysis. Expected keys:
            'parameters', 'quantiles', 'confidence_intervals'.
        reference : PeakfqSAResult
            Reference result to compare against.

        Returns
        -------
        ComparisonResult
            Detailed comparison with per-field differences.
        """
        # TODO: Implement full comparison
        raise NotImplementedError

    def compare_parameters(self, native: dict[str, Any], ref: PeakfqSAResult) -> dict[str, float]:
        """Compare LP3 parameters between native and reference.

        Returns
        -------
        dict[str, float]
            Parameter name to percent difference.
        """
        # TODO: Implement parameter comparison
        raise NotImplementedError

    def compare_quantiles(self, native: dict[str, Any], ref: PeakfqSAResult) -> dict[float, float]:
        """Compare quantile estimates between native and reference.

        Returns
        -------
        dict[float, float]
            AEP to percent difference.
        """
        # TODO: Implement quantile comparison
        raise NotImplementedError

    def compare_ci(self, native: dict[str, Any], ref: PeakfqSAResult) -> dict[float, float]:
        """Compare confidence intervals between native and reference.

        Returns
        -------
        dict[float, float]
            AEP to percent difference on bounds.
        """
        # TODO: Implement CI comparison
        raise NotImplementedError
