"""
hydrolib.regression - Nationwide regional regression and ROI peak flow estimation.

Provides a publication-agnostic workflow for comparing regional GLS regression
equation estimates with Region of Influence (ROI) estimates at ungaged sites,
supporting any state and any combination of StreamStats predictor variables.

Core classes
------------
:class:`HydrologicRegion`
    Describes one hydrologic region (area / zone) for which GLS regression
    equations exist.  Frozen dataclass — hashable and usable as dict key.

:class:`BasinCharacteristics`
    Stores basin attributes in a ``predictors`` dict keyed by USGS StreamStats
    parameter codes (``DRNAREA``, ``CSL1085LFP``, ``ELEV``, ``PRECIP``, …).
    Works for any state and any combination of predictor variables.

:class:`GlsEquation`
    One regression equation: ``log10(Q_T) = intercept + Σ coeff·log10(predictor)``.

:class:`RegressionTable`
    Multi-region, multi-AEP equation table.  Supports dict/CSV loading,
    nationwide ``merge()``, ``batch_estimate()``, and ``filter_by_state()``.

:class:`RoiAnalysis`
    Fetches NWIS peak records, fits B17C EMA, computes characteristic-space
    weights, and returns ROI-weighted quantile estimates.

:class:`PeakFlowComparator`
    Orchestrates GLS vs ROI comparison with B17C Appendix 8 inverse-variance
    weighting and karst-aware design estimate selection.

Predictor code constants
------------------------
``DRNAREA``, ``CSL1085LFP``, ``ELEV``, ``PRECIP``, ``I2``,
``IMPERV``, ``FOREST``, ``BFIPCT``, ``SLOPE``, ``LC06IMP``

State-specific region constants and table builders
--------------------------------------------------
``from hydrolib.regression.states import (``
``    TN_AREA1, TN_AREA2, TN_AREA3, TN_AREA4, build_tennessee_table,``
``    GA_BLUE_RIDGE, GA_VALLEY_RIDGE, GA_PIEDMONT, GA_COASTAL, build_georgia_table,``
``    MT_MOUNTAIN, MT_FOOTHILLS, MT_PLAINS, build_montana_table,``
``)``

Typical usage
-------------
::

    from hydrolib.regression import (
        BasinCharacteristics, RegressionTable, RoiAnalysis,
        PeakFlowComparator, SiteAssessment, KarstFlag,
        DRNAREA, CSL1085LFP,
    )
    from hydrolib.regression.states.tennessee import TN_AREA2, build_tennessee_table

    table = build_tennessee_table()
    basin = BasinCharacteristics(
        site_no="03606500",
        site_name="Big Sandy River at Bruceton, TN",
        region=TN_AREA2,
        predictors={DRNAREA: 779.0, CSL1085LFP: 2.4},
    )
    for aep, q in table.estimate_all_aeps(basin).items():
        print(f"Q{1/aep:.0f}: {q:,.0f} cfs")

Nationwide usage::

    from hydrolib.regression.regression_table import RegressionTable
    from hydrolib.regression.states import (
        build_tennessee_table, build_georgia_table, build_montana_table,
    )
    national = RegressionTable.merge(
        build_tennessee_table(), build_georgia_table(), build_montana_table(),
    )
    results = national.batch_estimate(basin_list, aep=0.01)
"""

from hydrolib.regression.basin_chars import (
    BFIPCT,
    CSL1085LFP,
    DRNAREA,
    ELEV,
    FOREST,
    I2,
    IMPERV,
    LC06IMP,
    PRECIP,
    SLOPE,
    BasinCharacteristics,
)
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

__all__ = [
    # Generic framework
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
    # StreamStats predictor code constants
    "DRNAREA",
    "CSL1085LFP",
    "ELEV",
    "PRECIP",
    "I2",
    "IMPERV",
    "FOREST",
    "BFIPCT",
    "SLOPE",
    "LC06IMP",
]
