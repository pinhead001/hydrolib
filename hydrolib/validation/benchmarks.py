"""
Bulletin 17C benchmark test cases.

Each benchmark is a named, self-contained test scenario with known inputs
and expected outputs from official USGS examples. Used for validating
HydroLib native EMA against published reference values.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Optional

from hydrolib.peakfqsa.config import PeakfqSAConfig
from hydrolib.peakfqsa.parsers import PeakfqSAResult
from hydrolib.validation.comparisons import ComparisonResult

logger = logging.getLogger(__name__)


@dataclass
class Benchmark:
    """A single benchmark test case.

    Parameters
    ----------
    name : str
        Short identifier for the benchmark.
    description : str
        Human-readable description.
    peaks : dict[int, float]
        Systematic peak flows {water_year: discharge}.
    historical : dict[int, float]
        Historical peaks {water_year: discharge}.
    thresholds : list[dict]
        Perception thresholds.
    expected_quantiles : dict[float, float]
        Expected AEP to discharge mapping.
    tolerance_pct : float
        Tolerance for passing.
    regional_skew : float
        Regional skew coefficient.
    regional_skew_sd : float
        Standard deviation of regional skew.
    begyear : int
        Analysis start year.
    endyear : int
        Analysis end year.
    """

    name: str = ""
    description: str = ""
    peaks: dict[int, float] = field(default_factory=dict)
    historical: dict[int, float] = field(default_factory=dict)
    thresholds: list[dict[str, Any]] = field(default_factory=list)
    expected_quantiles: dict[float, float] = field(default_factory=dict)
    tolerance_pct: float = 1.0
    regional_skew: float = -0.302
    regional_skew_sd: float = 0.55
    begyear: int = 0
    endyear: int = 0

    def run_native(self) -> Any:
        """Run benchmark using HydroLib native EMA.

        Returns
        -------
        AnalysisResult
            Native analysis result.
        """
        # TODO: Implement native benchmark run
        raise NotImplementedError

    def run_peakfqsa(self, config: PeakfqSAConfig) -> PeakfqSAResult:
        """Run benchmark using PeakfqSA subprocess.

        Parameters
        ----------
        config : PeakfqSAConfig
            PeakfqSA configuration.

        Returns
        -------
        PeakfqSAResult
            Reference result.
        """
        # TODO: Implement PeakfqSA benchmark run
        raise NotImplementedError

    def validate(self, config: Optional[PeakfqSAConfig] = None) -> ComparisonResult:
        """Run both native and reference, return comparison.

        Parameters
        ----------
        config : PeakfqSAConfig or None
            PeakfqSA config. If None, only native results are validated
            against expected_quantiles.

        Returns
        -------
        ComparisonResult
            Comparison between native and reference (or expected) results.
        """
        # TODO: Implement benchmark validation
        raise NotImplementedError


# Registry of available benchmarks
BENCHMARKS: dict[str, Benchmark] = {}


def run_all_benchmarks(
    peakfqsa_config: Optional[PeakfqSAConfig] = None,
) -> dict[str, ComparisonResult]:
    """Run all registered benchmarks.

    Parameters
    ----------
    peakfqsa_config : PeakfqSAConfig or None
        Configuration for PeakfqSA. If None, only native validation is run.

    Returns
    -------
    dict[str, ComparisonResult]
        Benchmark name to comparison result mapping.
    """
    # TODO: Implement run_all
    raise NotImplementedError


def print_benchmark_report(results: dict[str, ComparisonResult]) -> None:
    """Print a formatted report of benchmark results.

    Parameters
    ----------
    results : dict[str, ComparisonResult]
        Results from run_all_benchmarks.
    """
    # TODO: Implement report printer
    raise NotImplementedError
