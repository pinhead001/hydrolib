"""
hydrolib.regression - Regional regression and ROI peak flow estimation.

Implements the GLS/ROI comparison workflow for Tennessee streams per
SIR 2024-5130 (Ladd & Ensminger, 2025).

Typical usage
-------------
>>> from hydrolib.regression import BasinCharacteristics, HydrologicArea
>>> from hydrolib.regression import SIR2024_5130
>>> from hydrolib.regression import RoiAnalysis
>>> from hydrolib.regression import PeakFlowComparator

Workflow
--------
1. Define basin characteristics (DA, S1085) from StreamStats.
2. Compute GLS regression estimates via SIR2024_5130.
3. Run ROI analysis to fetch NWIS peak statistics for nearby gauges.
4. Compare GLS and ROI estimates; compute inverse-variance-weighted average
   per Bulletin 17C Appendix 8 framework.
"""

from hydrolib.regression.basin_chars import BasinCharacteristics, HydrologicArea
from hydrolib.regression.comparator import (
    ComparisonResult,
    KarstFlag,
    PeakFlowComparator,
    SiteAssessment,
    WeightedEstimate,
)
from hydrolib.regression.roi_analysis import RoiAnalysis, RoiSite
from hydrolib.regression.sir2024_5130 import SIR2024_5130, GlsEquation

__all__ = [
    "BasinCharacteristics",
    "HydrologicArea",
    "GlsEquation",
    "SIR2024_5130",
    "RoiSite",
    "RoiAnalysis",
    "ComparisonResult",
    "KarstFlag",
    "SiteAssessment",
    "WeightedEstimate",
    "PeakFlowComparator",
]
