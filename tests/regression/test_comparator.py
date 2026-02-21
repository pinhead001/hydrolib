"""
Tests for hydrolib.regression.comparator.

Exercises the GLS-vs-ROI comparison logic, inverse-variance weighting,
karst assessment logic, and report generation.
"""

from __future__ import annotations

import math
from typing import Dict, List
from unittest.mock import MagicMock, patch

import pytest

from hydrolib.core import log_pearson3_ppf
from hydrolib.regression.basin_chars import (
    CSL1085LFP,
    DRNAREA,
    BasinCharacteristics,
)
from hydrolib.regression.comparator import (
    ComparisonResult,
    KarstFlag,
    PeakFlowComparator,
    SiteAssessment,
    WeightedEstimate,
)
from hydrolib.regression.regression_table import RegressionTable
from hydrolib.regression.roi_analysis import STANDARD_AEPS, RoiAnalysis, RoiSite
from hydrolib.regression.states.tennessee import TN_AREA2

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def basin() -> BasinCharacteristics:
    return BasinCharacteristics(
        site_no="SITE01",
        site_name="Test Creek at Test, TN",
        region=TN_AREA2,
        predictors={DRNAREA: 85.3, CSL1085LFP: 4.7},
    )


@pytest.fixture
def gls_table(basin) -> RegressionTable:
    """Synthetic GLS table for Area 2 across all standard AEPs."""
    b0_by_aep = {
        0.5: 1.80,
        0.2: 2.10,
        0.1: 2.30,
        0.04: 2.55,
        0.02: 2.70,
        0.01: 2.90,
        0.005: 3.05,
        0.002: 3.25,
    }
    sep_by_aep = {
        0.5: 44.0,
        0.2: 40.0,
        0.1: 38.0,
        0.04: 37.0,
        0.02: 36.0,
        0.01: 35.0,
        0.005: 36.0,
        0.002: 38.0,
    }
    coeff_dict = {
        (TN_AREA2, aep): {
            "intercept": b0,
            "coefficients": {"DRNAREA": 0.75, "CSL1085LFP": 0.38},
            "sep_pct": sep_by_aep[aep],
            "pseudo_r2": 0.92,
            "eyr": 18.0,
        }
        for aep, b0 in b0_by_aep.items()
    }
    return RegressionTable.load_from_dict(coeff_dict)


def make_candidate_pool(target_basin: BasinCharacteristics) -> List[RoiSite]:
    """Create 8 candidate sites with pre-fitted LP3 parameters."""
    configs = [
        ("A", 70.0, 4.0, 2.70, 0.24, -0.08),
        ("B", 85.0, 4.9, 2.80, 0.26, -0.10),
        ("C", 100.0, 5.0, 2.85, 0.25, -0.09),
        ("D", 90.0, 4.5, 2.78, 0.27, -0.12),
        ("E", 120.0, 5.5, 2.90, 0.23, -0.07),
        ("F", 60.0, 3.8, 2.65, 0.28, -0.11),
        ("G", 80.0, 4.7, 2.75, 0.25, -0.09),
        ("H", 95.0, 4.8, 2.82, 0.26, -0.10),
    ]
    sites = []
    for sno, da, slope, mu, sig, skw in configs:
        b = BasinCharacteristics(
            site_no=sno,
            site_name=f"Site {sno}",
            region=TN_AREA2,
            predictors={DRNAREA: da, CSL1085LFP: slope},
        )
        site = RoiSite(
            site_no=sno,
            site_name=f"Site {sno}",
            basin=b,
            n_peaks=40,
            lp3_mean=mu,
            lp3_std=sig,
            lp3_skew=skw,
            record_source="nwis",
        )
        for aep in STANDARD_AEPS:
            site.at_site_quantiles[aep] = log_pearson3_ppf(1.0 - aep, mu, sig, skw)
        sites.append(site)
    return sites


@pytest.fixture
def roi_analysis(basin) -> RoiAnalysis:
    pool = make_candidate_pool(basin)
    roi = RoiAnalysis(target=basin, candidate_sites=pool, min_sites=3)
    roi.compute_weights()
    return roi


@pytest.fixture
def comparator(basin, gls_table, roi_analysis) -> PeakFlowComparator:
    return PeakFlowComparator(
        basin=basin,
        gls_table=gls_table,
        roi_analysis=roi_analysis,
        site_assessment=SiteAssessment(),
    )


# ---------------------------------------------------------------------------
# SiteAssessment tests
# ---------------------------------------------------------------------------


class TestSiteAssessment:
    def test_default_no_karst(self):
        sa = SiteAssessment()
        assert sa.karst_flag == KarstFlag.NO_KARST
        assert sa.karst_score == 0.0

    def test_karst_score_range(self):
        with pytest.raises(ValueError, match="karst_score"):
            SiteAssessment(karst_score=1.5)

    def test_confirmed_karst(self):
        sa = SiteAssessment(
            karst_flag=KarstFlag.CONFIRMED,
            notes="Sinkholes confirmed in basin",
            karst_score=1.0,
        )
        assert sa.karst_flag == KarstFlag.CONFIRMED


# ---------------------------------------------------------------------------
# ComparisonResult tests
# ---------------------------------------------------------------------------


class TestComparisonResult:
    def test_to_dict_contains_all_fields(self):
        cr = ComparisonResult(aep=0.01, return_period=100.0, gls_q=5000.0, roi_q=4200.0)
        d = cr.to_dict()
        assert "aep" in d
        assert "gls_q" in d
        assert "roi_q" in d

    def test_none_fields_allowed(self):
        cr = ComparisonResult(aep=0.01, return_period=100.0)
        assert cr.gls_q is None


# ---------------------------------------------------------------------------
# PeakFlowComparator.compare() tests
# ---------------------------------------------------------------------------


class TestPeakFlowComparatorCompare:
    def test_compare_returns_weighted_estimate(self, comparator):
        we = comparator.compare()
        assert isinstance(we, WeightedEstimate)

    def test_compare_all_aeps_present(self, comparator):
        we = comparator.compare()
        aeps_in_result = {cr.aep for cr in we.comparison_results}
        for aep in STANDARD_AEPS:
            assert aep in aeps_in_result

    def test_gls_estimates_positive(self, comparator):
        we = comparator.compare()
        for cr in we.comparison_results:
            assert cr.gls_q is not None
            assert cr.gls_q > 0

    def test_roi_estimates_positive(self, comparator):
        we = comparator.compare()
        for cr in we.comparison_results:
            assert cr.roi_q is not None
            assert cr.roi_q > 0

    def test_weighted_estimates_positive(self, comparator):
        we = comparator.compare()
        for cr in we.comparison_results:
            assert cr.weighted_q is not None
            assert cr.weighted_q > 0

    def test_100yr_greater_than_2yr(self, comparator):
        we = comparator.compare()
        q_map = {cr.aep: cr.weighted_q for cr in we.comparison_results if cr.weighted_q}
        assert q_map[0.01] > q_map[0.5]

    def test_diff_pct_computed(self, comparator):
        we = comparator.compare()
        for cr in we.comparison_results:
            assert cr.diff_pct is not None

    def test_sep_pct_populated(self, comparator):
        we = comparator.compare()
        for cr in we.comparison_results:
            assert cr.gls_sep_pct is not None
            assert cr.gls_sep_pct > 0

    def test_ratio_roi_to_gls_positive(self, comparator):
        we = comparator.compare()
        for cr in we.comparison_results:
            assert cr.ratio_roi_to_gls is not None
            assert cr.ratio_roi_to_gls > 0


# ---------------------------------------------------------------------------
# Inverse-variance-weighted average (B17C Appendix 8) unit tests
# ---------------------------------------------------------------------------


class TestB17CWeightedAverage:
    """Tests for PeakFlowComparator._b17c_weighted_average static method."""

    def test_equal_variance_is_arithmetic_mean_in_log(self):
        """With equal variances, weighted log10 average = arithmetic mean of log10."""
        q1, q2 = 5000.0, 3000.0
        v = 0.04
        w_q, w_v = PeakFlowComparator._b17c_weighted_average(q1, q2, v, v)
        expected_log = (math.log10(q1) + math.log10(q2)) / 2.0
        assert math.isclose(math.log10(w_q), expected_log, rel_tol=1e-9)

    def test_lower_variance_gets_more_weight(self):
        """Estimate with lower variance should pull the weighted average toward it."""
        q_low_var = 8000.0  # high-confidence estimate
        q_high_var = 4000.0  # low-confidence estimate
        v_low = 0.01
        v_high = 0.10
        w_q, _ = PeakFlowComparator._b17c_weighted_average(q_low_var, q_high_var, v_low, v_high)
        # Weighted average should be closer to q_low_var (low variance)
        assert math.log10(w_q) > (math.log10(q_low_var) + math.log10(q_high_var)) / 2.0

    def test_combined_variance_less_than_either(self):
        """Combined variance must be less than either input variance."""
        q1, q2 = 5000.0, 4000.0
        v1, v2 = 0.04, 0.06
        _, v_comb = PeakFlowComparator._b17c_weighted_average(q1, q2, v1, v2)
        assert v_comb < v1
        assert v_comb < v2

    def test_missing_gls_returns_roi(self):
        q_roi = 4000.0
        v_roi = 0.04
        w_q, _ = PeakFlowComparator._b17c_weighted_average(None, q_roi, None, v_roi)
        assert w_q == q_roi

    def test_missing_roi_returns_gls(self):
        q_gls = 5000.0
        v_gls = 0.04
        w_q, _ = PeakFlowComparator._b17c_weighted_average(q_gls, None, v_gls, None)
        assert w_q == q_gls

    def test_both_missing_returns_none(self):
        w_q, w_v = PeakFlowComparator._b17c_weighted_average(None, None, None, None)
        assert w_q is None
        assert w_v is None

    def test_zero_variance_not_used(self):
        """Zero variance would cause division by zero; fall back gracefully."""
        w_q, w_v = PeakFlowComparator._b17c_weighted_average(5000.0, 4000.0, 0.0, 0.0)
        # Both variances are 0 → no valid weighting; fall back to GLS (first available)
        assert w_q in (5000.0, 4000.0) or w_q is None


# ---------------------------------------------------------------------------
# Karst assessment / design estimate tests
# ---------------------------------------------------------------------------


class TestKarstDesignEstimate:
    def _make_we(
        self, gls_q: float, roi_q: float, karst: KarstFlag, ks: float = 0.0
    ) -> WeightedEstimate:
        """Helper: create a WeightedEstimate with one ComparisonResult."""
        wq, wv = PeakFlowComparator._b17c_weighted_average(gls_q, roi_q, 0.04, 0.04)
        cr = ComparisonResult(
            aep=0.01,
            return_period=100.0,
            gls_q=gls_q,
            roi_q=roi_q,
            gls_variance=0.04,
            roi_variance=0.04,
            weighted_q=wq,
            weighted_variance=wv,
        )
        sa = SiteAssessment(karst_flag=karst, karst_score=ks)
        return WeightedEstimate(comparison_results=[cr], site_assessment=sa)

    def test_no_karst_returns_weighted(self):
        """No karst → design estimate equals inverse-variance-weighted Q."""
        gls, roi = 6000.0, 4000.0
        we = self._make_we(gls, roi, KarstFlag.NO_KARST)
        design = we.design_estimate(0.01)
        # Should be the weighted average (in log space, equal variances → geometric mean)
        expected_log = (math.log10(gls) + math.log10(roi)) / 2.0
        assert math.isclose(math.log10(design), expected_log, rel_tol=1e-6)

    def test_confirmed_karst_returns_lower(self):
        """Confirmed karst → design estimate is the lower of GLS and ROI."""
        gls, roi = 6000.0, 4000.0
        we = self._make_we(gls, roi, KarstFlag.CONFIRMED)
        design = we.design_estimate(0.01)
        assert design == min(gls, roi)

    def test_confirmed_karst_roi_higher_returns_gls(self):
        """Confirmed karst → still returns lower even if ROI > GLS."""
        gls, roi = 3500.0, 7000.0
        we = self._make_we(gls, roi, KarstFlag.CONFIRMED)
        design = we.design_estimate(0.01)
        assert design == min(gls, roi)

    def test_possible_karst_blends_toward_lower(self):
        """Possible karst with ks=0.5 should blend midway toward lower bound."""
        gls, roi = 6000.0, 4000.0
        we = self._make_we(gls, roi, KarstFlag.POSSIBLE, ks=0.5)
        design = we.design_estimate(0.01)
        w_avg, _ = PeakFlowComparator._b17c_weighted_average(gls, roi, 0.04, 0.04)
        lower = min(gls, roi)
        expected = w_avg * 0.5 + lower * 0.5
        assert math.isclose(design, expected, rel_tol=1e-6)

    def test_possible_karst_zero_score_equals_no_karst(self):
        """Possible karst with ks=0 should equal the no-karst weighted average."""
        gls, roi = 6000.0, 4000.0
        we_possible = self._make_we(gls, roi, KarstFlag.POSSIBLE, ks=0.0)
        we_none = self._make_we(gls, roi, KarstFlag.NO_KARST)
        assert math.isclose(
            we_possible.design_estimate(0.01),
            we_none.design_estimate(0.01),
            rel_tol=1e-6,
        )

    def test_missing_aep_returns_none(self):
        """Querying an AEP not in the result list should return None."""
        cr = ComparisonResult(aep=0.01, return_period=100.0, gls_q=5000.0, roi_q=4000.0)
        we = WeightedEstimate(comparison_results=[cr])
        assert we.design_estimate(0.002) is None


# ---------------------------------------------------------------------------
# PeakFlowComparator.report() and to_dataframe() tests
# ---------------------------------------------------------------------------


class TestComparatorOutputs:
    def test_report_contains_site_info(self, comparator):
        we = comparator.compare()
        report = comparator.report(we)
        assert "SITE01" in report
        assert "Area" in report

    def test_report_has_karst_note_for_confirmed(self, basin, gls_table, roi_analysis):
        sa = SiteAssessment(
            karst_flag=KarstFlag.CONFIRMED,
            notes="Sinkhole observed near intake",
        )
        c = PeakFlowComparator(basin, gls_table, roi_analysis, sa)
        we = c.compare()
        report = c.report(we)
        assert "karst" in report.lower() or "CONFIRMED" in report

    def test_to_dataframe_shape(self, comparator):
        we = comparator.compare()
        df = comparator.to_dataframe(we)
        assert len(df) == len(STANDARD_AEPS)
        assert "design_q_cfs" in df.columns

    def test_to_dataframe_sorted_descending_aep(self, comparator):
        we = comparator.compare()
        df = comparator.to_dataframe(we)
        aeps = df["aep"].tolist()
        assert aeps == sorted(aeps, reverse=True)

    def test_weighted_estimate_to_dataframe(self, comparator):
        we = comparator.compare()
        df = we.to_dataframe()
        assert not df.empty
        assert "gls_q" in df.columns
        assert "roi_q" in df.columns
