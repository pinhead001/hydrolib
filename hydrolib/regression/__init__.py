"""
hydrolib.regression - Regional regression and ROI peak flow estimation.

Provides a publication-agnostic workflow for comparing regional regression
equation estimates with Region of Influence (ROI) estimates at ungaged sites.

Core generic classes
--------------------
:class:`HydrologicRegion`
    Describes a hydrologic region (area / zone) for which GLS regression
    equations exist.  Hashable; usable as a dict key.

:class:`BasinCharacteristics`
    Stores basin attributes in a flexible ``predictors`` dict keyed by
    USGS StreamStats parameter codes (``DRNAREA``, ``CSL1085LFP``, etc.).
    Works for any state and any combination of predictor variables.

:class:`GlsEquation`
    Single regression equation: ``log10(Q_T) = intercept + Σ coeff·log10(predictor)``.

:class:`RegressionTable`
    Publication-agnostic container for a full set of GLS equations indexed
    by (region, AEP).  Loads from dict or wide-format CSV.

:class:`RoiAnalysis`
    Fetches NWIS peak records, fits B17C EMA, computes characteristic-space
    weights, and returns ROI-weighted quantile estimates.

:class:`PeakFlowComparator`
    Orchestrates GLS vs ROI comparison with B17C Appendix 8 inverse-variance
    weighting and karst-aware design estimate selection.

Tennessee-specific helpers
--------------------------
``from hydrolib.regression.sir2024_5130 import (
    HydrologicArea, SIR2024_5130, TN_AREA1, TN_AREA2, TN_AREA3, TN_AREA4
)``

Typical usage
-------------
>>> from hydrolib.regression import (
...     HydrologicRegion, BasinCharacteristics,
...     RegressionTable, RoiAnalysis, RoiSite,
...     PeakFlowComparator, SiteAssessment, KarstFlag,
... )
>>> from hydrolib.regression.sir2024_5130 import HydrologicArea, SIR2024_5130
"""

from hydrolib.regression.basin_chars import BasinCharacteristics
from hydrolib.regression.comparator import (
    ComparisonResult,
    KarstFlag,
    PeakFlowComparator,
    SiteAssessment,
    WeightedEstimate,
)
from hydrolib.regression.region import HydrologicRegion, RegionRegistry
from hydrolib.regression.regression_table import GlsEquation, RegressionTable
from hydrolib.regression.roi_analysis import RoiAnalysis, RoiSite
from hydrolib.regression.sir2024_5130 import (
    SIR2024_5130,
    TN_AREA1,
    TN_AREA2,
    TN_AREA3,
    TN_AREA4,
    TN_REGIONS,
    HydrologicArea,
)

__all__ = [
    # Generic
    "HydrologicRegion",
    "RegionRegistry",
    "BasinCharacteristics",
    "GlsEquation",
    "RegressionTable",
    "RoiSite",
    "RoiAnalysis",
    "ComparisonResult",
    "KarstFlag",
    "SiteAssessment",
    "WeightedEstimate",
    "PeakFlowComparator",
    # Tennessee / SIR 2024-5130 convenience exports
    "HydrologicArea",
    "SIR2024_5130",
    "TN_AREA1",
    "TN_AREA2",
    "TN_AREA3",
    "TN_AREA4",
    "TN_REGIONS",
]
