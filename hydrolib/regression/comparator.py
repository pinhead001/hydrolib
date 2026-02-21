"""
hydrolib.regression.comparator - GLS regression vs ROI comparison engine.

Combines estimates from two sources:

1. **GLS regression** (SIR 2024-5130) — regional power-law equation applied
   to StreamStats basin characteristics.
2. **ROI analysis** — weighted average of at-site B17C EMA estimates from
   nearby NWIS-gauged sites.

Then computes an inverse-variance-weighted average per the framework
described in Bulletin 17C Appendix 8 (England et al., 2019, p. A8-1—A8-5).

References
----------
England, J.F., Jr., Cohn, T.A., Faber, B.A., Stedinger, J.R.,
Thomas, W.O., Jr., Veilleux, A.G., Kiang, J.E., and Mason, R.R., Jr.,
2018, Guidelines for determining flood flow frequency—Bulletin 17C:
U.S. Geological Survey Techniques and Methods, book 4, chap. B5, 148 p.,
https://doi.org/10.3133/tm4B5.

Tasker, G.D., and Stedinger, J.R., 1989, An operational GLS model for
hydrologic regression: Journal of Hydrology, v. 111, p. 361–375.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from hydrolib.regression.basin_chars import BasinCharacteristics
from hydrolib.regression.regression_table import RegressionTable
from hydrolib.regression.roi_analysis import STANDARD_AEPS, RoiAnalysis

log = logging.getLogger(__name__)

STANDARD_RETURN_PERIODS = tuple(1.0 / a for a in STANDARD_AEPS)


# ---------------------------------------------------------------------------
# Karst / site reconnaissance flags
# ---------------------------------------------------------------------------


class KarstFlag(Enum):
    """
    Site reconnaissance assessment for karst or losing-stream behaviour.

    Design implications (per Bulletin 17C and TDOT practice)
    ----------------------------------------------------------
    * ``NO_KARST``  — no evidence; use the inverse-variance-weighted average.
    * ``POSSIBLE``  — some indicators; weight ROI estimates more heavily.
    * ``CONFIRMED`` — confirmed karst or losing stream; favour the **lower**
      (ROI) estimate for conservative flood-frequency design, per standard
      TDOT / USACE practice for karst terrains.

    Notes
    -----
    Bulletin 17C does not prescribe a karst-specific procedure; the
    ``CONFIRMED`` recommendation follows common engineering judgment that
    karst drainage losses cause USGS regression equations calibrated on
    gauged sites (which may include karst losses) to overpredict flows for
    sites with more severe losing-stream behaviour.
    """

    NO_KARST = auto()
    POSSIBLE = auto()
    CONFIRMED = auto()


@dataclass
class SiteAssessment:
    """
    Field or GIS-based site assessment used to inform estimate selection.

    Parameters
    ----------
    karst_flag : KarstFlag
        Karst/losing-stream classification.
    notes : str
        Free-text reconnaissance notes (sinkholes, springs, geology, etc.).
    karst_score : float
        Qualitative karst severity score (0.0 = none, 1.0 = severe).
        Used to blend GLS and ROI estimates when karst is possible but
        not confirmed.
    """

    karst_flag: KarstFlag = KarstFlag.NO_KARST
    notes: str = ""
    karst_score: float = 0.0

    def __post_init__(self) -> None:
        if not (0.0 <= self.karst_score <= 1.0):
            raise ValueError(f"karst_score must be in [0, 1], got {self.karst_score}")


# ---------------------------------------------------------------------------
# Per-AEP comparison result
# ---------------------------------------------------------------------------


@dataclass
class ComparisonResult:
    """
    Peak flow comparison result for a single AEP.

    Attributes
    ----------
    aep : float
        Annual exceedance probability.
    return_period : float
        Return period in years (1/AEP).
    gls_q : float or None
        GLS regression peak flow estimate (cfs).
    roi_q : float or None
        ROI-weighted peak flow estimate (cfs).
    gls_variance : float or None
        Variance of GLS estimate in log10(cfs) space.
    roi_variance : float or None
        Variance of ROI estimate in log10(cfs) space.
    weighted_q : float or None
        Inverse-variance-weighted peak flow estimate (cfs).
    weighted_variance : float or None
        Variance of the weighted estimate (log10 space).
    gls_sep_pct : float or None
        Average standard error of prediction (%) for GLS equation.
    ratio_roi_to_gls : float or None
        Ratio ROI / GLS (< 1 → ROI predicts lower flows than GLS).
    diff_pct : float or None
        Percent difference: (GLS - ROI) / ((GLS + ROI)/2) * 100.
    """

    aep: float
    return_period: float
    gls_q: Optional[float] = None
    roi_q: Optional[float] = None
    gls_variance: Optional[float] = None
    roi_variance: Optional[float] = None
    weighted_q: Optional[float] = None
    weighted_variance: Optional[float] = None
    gls_sep_pct: Optional[float] = None
    ratio_roi_to_gls: Optional[float] = None
    diff_pct: Optional[float] = None

    def to_dict(self) -> dict:
        """Return as a plain dictionary."""
        return {k: v for k, v in self.__dict__.items()}


# ---------------------------------------------------------------------------
# Weighted estimate (B17C Appendix 8)
# ---------------------------------------------------------------------------


@dataclass
class WeightedEstimate:
    """
    Inverse-variance-weighted combination of GLS and ROI estimates.

    Implements Equation A8-1 of Bulletin 17C (England et al., 2018, p. A8-2):

    .. math::

        \\log_{10}(\\hat{Q}) = \\frac{
            \\log_{10}(Q_{\\text{GLS}}) / V_{\\text{GLS}}
            + \\log_{10}(Q_{\\text{ROI}}) / V_{\\text{ROI}}
        }{
            1 / V_{\\text{GLS}} + 1 / V_{\\text{ROI}}
        }

    where :math:`V` denotes variance in log10 space.

    The combined variance is:

    .. math::

        V_{\\text{comb}} = \\frac{1}{1/V_{\\text{GLS}} + 1/V_{\\text{ROI}}}

    Parameters
    ----------
    comparison_results : list of ComparisonResult
        One per AEP.
    site_assessment : SiteAssessment
        Site reconnaissance information.
    """

    comparison_results: List[ComparisonResult] = field(default_factory=list)
    site_assessment: SiteAssessment = field(default_factory=SiteAssessment)

    def design_estimate(self, aep: float) -> Optional[float]:
        """
        Return the recommended design estimate for a given AEP.

        Decision logic:

        * ``KarstFlag.CONFIRMED`` → return the *lower* of GLS and ROI
          (conservative for losing-stream basins).
        * ``KarstFlag.POSSIBLE`` → return the inverse-variance-weighted
          average, then blend toward the lower bound by ``karst_score``.
        * ``KarstFlag.NO_KARST`` → return the inverse-variance-weighted
          average directly.

        Parameters
        ----------
        aep : float
            Annual exceedance probability.

        Returns
        -------
        float or None
        """
        cr = self._get_cr(aep)
        if cr is None:
            return None

        flag = self.site_assessment.karst_flag
        if flag == KarstFlag.CONFIRMED:
            # Use lower of the two estimates (conservative for karst)
            if cr.gls_q is not None and cr.roi_q is not None:
                lower = min(cr.gls_q, cr.roi_q)
                log.info(
                    "AEP=%.4f: CONFIRMED karst → using lower estimate %.1f cfs "
                    "(GLS=%.1f, ROI=%.1f)",
                    aep,
                    lower,
                    cr.gls_q,
                    cr.roi_q,
                )
                return lower
            return cr.roi_q or cr.gls_q

        elif flag == KarstFlag.POSSIBLE:
            w_avg = cr.weighted_q
            if w_avg is None:
                return cr.roi_q or cr.gls_q
            lower = (
                min(v for v in [cr.gls_q, cr.roi_q] if v is not None)
                if cr.gls_q and cr.roi_q
                else w_avg
            )
            ks = self.site_assessment.karst_score
            blended = w_avg * (1 - ks) + lower * ks
            log.info(
                "AEP=%.4f: POSSIBLE karst (score=%.2f) → blended %.1f cfs",
                aep,
                ks,
                blended,
            )
            return blended

        else:
            return cr.weighted_q or cr.roi_q or cr.gls_q

    def _get_cr(self, aep: float) -> Optional[ComparisonResult]:
        for cr in self.comparison_results:
            if math.isclose(cr.aep, aep, rel_tol=1e-6):
                return cr
        return None

    def to_dataframe(self) -> pd.DataFrame:
        """Return all comparison results as a pandas DataFrame."""
        if not self.comparison_results:
            return pd.DataFrame()
        return pd.DataFrame([cr.to_dict() for cr in self.comparison_results])


# ---------------------------------------------------------------------------
# PeakFlowComparator – the main orchestrator
# ---------------------------------------------------------------------------


class PeakFlowComparator:
    """
    Orchestrates the GLS-regression vs ROI comparison workflow.

    Parameters
    ----------
    basin : BasinCharacteristics
        Basin characteristics for the target site.
    gls_table : RegressionTable
        Populated GLS equation table (any state / publication).
    roi_analysis : RoiAnalysis
        Fully-run ROI analysis (peaks fetched, EMA fitted, weights computed).
    site_assessment : SiteAssessment, optional
        Site reconnaissance information (karst flags, notes).

    Examples
    --------
    >>> from hydrolib.regression import (
    ...     BasinCharacteristics, HydrologicRegion,
    ...     RegressionTable, RoiAnalysis, PeakFlowComparator,
    ...     SiteAssessment, KarstFlag,
    ... )
    >>> from hydrolib.regression.sir2024_5130 import HydrologicArea, SIR2024_5130
    >>>
    >>> basin = BasinCharacteristics(
    ...     site_no="UNGAGED01", site_name="Example Creek",
    ...     region=HydrologicArea.AREA2,
    ...     predictors={"DRNAREA": 85.3, "CSL1085LFP": 4.7},
    ... )
    >>> assessment = SiteAssessment(
    ...     karst_flag=KarstFlag.POSSIBLE,
    ...     notes="Several sinkholes observed within 2 miles of channel",
    ...     karst_score=0.4,
    ... )
    >>> comparator = PeakFlowComparator(basin, gls_table, roi, assessment)
    >>> result = comparator.compare()
    >>> print(result.to_dataframe())
    """

    def __init__(
        self,
        basin: BasinCharacteristics,
        gls_table: RegressionTable,
        roi_analysis: RoiAnalysis,
        site_assessment: Optional[SiteAssessment] = None,
    ) -> None:
        self.basin = basin
        self.gls_table = gls_table
        self.roi_analysis = roi_analysis
        self.site_assessment = site_assessment or SiteAssessment()

    # ------------------------------------------------------------------
    # Main comparison method
    # ------------------------------------------------------------------

    def compare(
        self,
        aeps: Optional[Tuple[float, ...]] = None,
    ) -> WeightedEstimate:
        """
        Run the full GLS vs ROI comparison for all (or specified) AEPs.

        Sequence:

        1. Compute GLS estimates from the SIR 2024-5130 regression table.
        2. Retrieve ROI-weighted quantile estimates.
        3. Compute inverse-variance-weighted average per B17C Appendix 8.
        4. Return a :class:`WeightedEstimate` containing per-AEP results.

        Parameters
        ----------
        aeps : tuple of float, optional
            AEPs to compute (default: all standard AEPs).

        Returns
        -------
        WeightedEstimate
        """
        aeps = aeps or STANDARD_AEPS

        # GLS estimates
        gls_estimates = self._compute_gls(aeps)

        # ROI estimates (compute_roi_quantiles already called internally)
        roi_estimates = self.roi_analysis.compute_roi_quantiles()
        roi_variances = self.roi_analysis.roi_variance()

        results: List[ComparisonResult] = []
        for aep in aeps:
            gls_q = gls_estimates.get(aep)
            roi_q = roi_estimates.get(aep)
            gls_var = self._gls_variance(aep)
            roi_var = roi_variances.get(aep)

            weighted_q, weighted_var = self._b17c_weighted_average(gls_q, roi_q, gls_var, roi_var)

            ratio = (roi_q / gls_q) if (gls_q and roi_q and gls_q > 0) else None
            diff = None
            if gls_q is not None and roi_q is not None and (gls_q + roi_q) > 0:
                diff = (gls_q - roi_q) / ((gls_q + roi_q) / 2.0) * 100.0

            sep_pct = None
            try:
                eq = self.gls_table.get_equation(self.basin.region, aep)
                sep_pct = eq.sep_pct
            except KeyError:
                pass

            results.append(
                ComparisonResult(
                    aep=aep,
                    return_period=1.0 / aep,
                    gls_q=gls_q,
                    roi_q=roi_q,
                    gls_variance=gls_var,
                    roi_variance=(
                        roi_var if (roi_var is not None and not math.isnan(roi_var)) else None
                    ),
                    weighted_q=weighted_q,
                    weighted_variance=weighted_var,
                    gls_sep_pct=sep_pct,
                    ratio_roi_to_gls=ratio,
                    diff_pct=diff,
                )
            )

        return WeightedEstimate(
            comparison_results=results,
            site_assessment=self.site_assessment,
        )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _compute_gls(self, aeps: Tuple[float, ...]) -> Dict[float, float]:
        """Compute GLS regression estimates for the target basin."""
        estimates: Dict[float, float] = {}
        for aep in aeps:
            try:
                estimates[aep] = self.gls_table.estimate(self.basin, aep)
            except (KeyError, ValueError) as exc:
                log.warning("GLS estimate failed at AEP=%.4f: %s", aep, exc)
        return estimates

    def _gls_variance(self, aep: float) -> Optional[float]:
        """Return GLS prediction variance (log10 space) for the given AEP."""
        try:
            return self.gls_table.get_variance(self.basin.region, aep)
        except KeyError:
            return None

    @staticmethod
    def _b17c_weighted_average(
        gls_q: Optional[float],
        roi_q: Optional[float],
        gls_var: Optional[float],
        roi_var: Optional[float],
    ) -> Tuple[Optional[float], Optional[float]]:
        """
        Inverse-variance-weighted average in log10 space (B17C Appendix 8).

        Returns (weighted_q_cfs, weighted_variance_log10).

        Falls back to the available single estimate if one is missing,
        or None if both are missing.
        """
        have_gls = gls_q is not None and gls_q > 0 and gls_var is not None and gls_var > 0
        have_roi = roi_q is not None and roi_q > 0 and roi_var is not None and roi_var > 0

        if have_gls and have_roi:
            log_gls = math.log10(gls_q)
            log_roi = math.log10(roi_q)
            w_gls = 1.0 / gls_var
            w_roi = 1.0 / roi_var
            w_total = w_gls + w_roi
            log_weighted = (w_gls * log_gls + w_roi * log_roi) / w_total
            var_weighted = 1.0 / w_total
            return 10**log_weighted, var_weighted

        elif have_gls and not have_roi:
            return gls_q, gls_var
        elif have_roi and not have_gls:
            return roi_q, roi_var
        elif gls_q is not None and gls_q > 0:
            # GLS available but no variance info
            return gls_q, None
        elif roi_q is not None and roi_q > 0:
            return roi_q, None
        return None, None

    # ------------------------------------------------------------------
    # Reporting
    # ------------------------------------------------------------------

    def report(self, weighted: WeightedEstimate) -> str:
        """
        Generate a text comparison report.

        Parameters
        ----------
        weighted : WeightedEstimate
            Result from :meth:`compare`.

        Returns
        -------
        str
            Multi-line formatted report.
        """
        lines = [
            "=" * 70,
            "GLS Regression vs ROI Analysis — Peak Flow Comparison",
            f"Site  : {self.basin.site_no} — {self.basin.site_name}",
            f"Area  : {self.basin.region.label}",
            f"DA    : {self.basin.predictors.get('DRNAREA', float('nan')):.2f} sq mi",
        ]
        slope = self.basin.predictors.get("CSL1085LFP")
        if slope is not None:
            lines.append(f"S1085 : {slope:.3f} ft/mi")

        karst = weighted.site_assessment.karst_flag.name
        lines.append(f"Karst : {karst}  ({weighted.site_assessment.notes or 'no notes'})")
        lines.append("-" * 70)

        header = (
            f"{'AEP':>7}  {'T(yr)':>6}  {'GLS(cfs)':>10}  {'ROI(cfs)':>10}  "
            f"{'Weighted':>10}  {'Diff%':>7}  {'SEP%':>6}  {'Design':>10}"
        )
        lines.append(header)
        lines.append("-" * 70)

        for cr in sorted(weighted.comparison_results, key=lambda x: x.aep, reverse=True):
            design = weighted.design_estimate(cr.aep)
            gls_str = f"{cr.gls_q:.0f}" if cr.gls_q else "N/A"
            roi_str = f"{cr.roi_q:.0f}" if cr.roi_q else "N/A"
            wt_str = f"{cr.weighted_q:.0f}" if cr.weighted_q else "N/A"
            diff_str = f"{cr.diff_pct:.1f}" if cr.diff_pct is not None else "N/A"
            sep_str = f"{cr.gls_sep_pct:.1f}" if cr.gls_sep_pct is not None else "N/A"
            des_str = f"{design:.0f}" if design is not None else "N/A"

            lines.append(
                f"{cr.aep:>7.4f}  {cr.return_period:>6.0f}  "
                f"{gls_str:>10}  {roi_str:>10}  "
                f"{wt_str:>10}  {diff_str:>7}  {sep_str:>6}  {des_str:>10}"
            )

        lines.append("=" * 70)

        # Notes
        if weighted.site_assessment.karst_flag == KarstFlag.CONFIRMED:
            lines.append(
                "NOTE: Karst confirmed.  Design estimates use the LOWER of GLS "
                "and ROI per standard practice for losing-stream basins."
            )
        elif weighted.site_assessment.karst_flag == KarstFlag.POSSIBLE:
            ks = weighted.site_assessment.karst_score
            lines.append(
                f"NOTE: Karst possible (score={ks:.2f}).  Design estimates blend "
                "inverse-variance-weighted average toward lower bound."
            )

        lines.append("REF: B17C Appendix 8 inverse-variance weighting | SIR 2024-5130 GLS")
        return "\n".join(lines)

    def to_dataframe(self, weighted: WeightedEstimate) -> pd.DataFrame:
        """
        Return comparison results as a tidy pandas DataFrame with design column.

        Parameters
        ----------
        weighted : WeightedEstimate

        Returns
        -------
        pd.DataFrame
        """
        df = weighted.to_dataframe()
        if df.empty:
            return df
        df["design_q_cfs"] = df["aep"].apply(lambda a: weighted.design_estimate(a))
        return df.sort_values("aep", ascending=False).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Convenience: re-export KarstFlag for backwards compatibility
# ---------------------------------------------------------------------------

__all__ = [
    "KarstFlag",
    "SiteAssessment",
    "ComparisonResult",
    "WeightedEstimate",
    "PeakFlowComparator",
]
