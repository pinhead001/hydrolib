"""
Tests for hydrolib.regression.roi_analysis.

Uses mocking for NWIS HTTP calls and a simplified EMA (via B17CEngine)
to validate the ROI distance/weighting and quantile computation logic.
"""

from __future__ import annotations

import math
import textwrap
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from hydrolib.core import log_pearson3_ppf
from hydrolib.regression.basin_chars import (
    CSL1085LFP,
    DRNAREA,
    BasinCharacteristics,
)
from hydrolib.regression.roi_analysis import (
    STANDARD_AEPS,
    RoiAnalysis,
    RoiSite,
    fetch_nwis_peak_data,
)
from hydrolib.regression.sir2024_5130 import TN_AREA2, TN_AREA3

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def target() -> BasinCharacteristics:
    """Target ungaged site in Area 2."""
    return BasinCharacteristics(
        site_no="UNGAGED",
        site_name="Example Creek",
        region=TN_AREA2,
        predictors={DRNAREA: 100.0, CSL1085LFP: 5.0},
        latitude=36.0,
        longitude=-86.5,
    )


def make_roi_site(
    site_no: str,
    da: float,
    slope: float,
    lp3_mean: float = 2.8,
    lp3_std: float = 0.25,
    lp3_skew: float = -0.1,
    lat: float = 36.0,
    lon: float = -86.5,
    n_peaks: int = 40,
) -> RoiSite:
    """Helper to create an RoiSite with pre-populated LP3 parameters."""
    basin = BasinCharacteristics(
        site_no=site_no,
        site_name=f"Site {site_no}",
        region=TN_AREA2,
        predictors={DRNAREA: da, CSL1085LFP: slope},
        latitude=lat,
        longitude=lon,
    )
    site = RoiSite(
        site_no=site_no,
        site_name=f"Site {site_no}",
        basin=basin,
        weight=0.0,
        n_peaks=n_peaks,
        lp3_mean=lp3_mean,
        lp3_std=lp3_std,
        lp3_skew=lp3_skew,
        record_source="nwis",
    )
    for aep in STANDARD_AEPS:
        site.at_site_quantiles[aep] = log_pearson3_ppf(1.0 - aep, lp3_mean, lp3_std, lp3_skew)
    return site


@pytest.fixture
def candidate_pool(target) -> list:
    """Pool of 8 candidate sites with varying DA and slope."""
    return [
        make_roi_site("A", da=80.0, slope=4.5, lat=35.9, lon=-86.4),
        make_roi_site("B", da=90.0, slope=5.1, lat=36.0, lon=-86.6),
        make_roi_site("C", da=110.0, slope=5.3, lat=36.1, lon=-86.5),
        make_roi_site("D", da=200.0, slope=3.0, lat=36.3, lon=-86.8),
        make_roi_site("E", da=50.0, slope=8.0, lat=35.8, lon=-86.3),
        make_roi_site("F", da=150.0, slope=4.0, lat=36.2, lon=-86.7),
        make_roi_site("G", da=120.0, slope=5.0, lat=36.0, lon=-86.5),
        make_roi_site("H", da=95.0, slope=4.8, lat=36.0, lon=-86.4),
    ]


@pytest.fixture
def roi(target, candidate_pool) -> RoiAnalysis:
    return RoiAnalysis(
        target=target,
        candidate_sites=candidate_pool,
        weight_method="char_space",
        min_sites=3,
        max_sites=8,
    )


# ---------------------------------------------------------------------------
# RoiSite tests
# ---------------------------------------------------------------------------


class TestRoiSite:
    def test_has_lp3_true(self):
        site = make_roi_site("X", 100.0, 5.0)
        assert site.has_lp3

    def test_has_lp3_false(self):
        basin = BasinCharacteristics(
            site_no="X",
            site_name="X",
            region=TN_AREA3,
            predictors={DRNAREA: 10.0},
        )
        site = RoiSite(site_no="X", site_name="X", basin=basin)
        assert not site.has_lp3

    def test_quantile_from_cache(self):
        site = make_roi_site("X", 100.0, 5.0)
        q = site.quantile(0.01)
        assert q is not None and q > 0

    def test_quantile_computed_from_lp3(self):
        """If AEP not in cache but LP3 params exist, it should compute on-the-fly."""
        basin = BasinCharacteristics(
            site_no="X",
            site_name="X",
            region=TN_AREA3,
            predictors={DRNAREA: 100.0},
        )
        site = RoiSite(
            site_no="X",
            site_name="X",
            basin=basin,
            lp3_mean=2.8,
            lp3_std=0.25,
            lp3_skew=-0.1,
        )
        # No cached quantiles
        q = site.quantile(0.01)
        assert q is not None and q > 0

    def test_quantile_none_without_lp3(self):
        basin = BasinCharacteristics(
            site_no="X",
            site_name="X",
            region=TN_AREA3,
            predictors={DRNAREA: 10.0},
        )
        site = RoiSite(site_no="X", site_name="X", basin=basin)
        assert site.quantile(0.01) is None


# ---------------------------------------------------------------------------
# fetch_nwis_peak_data tests (mocked HTTP)
# ---------------------------------------------------------------------------

MOCK_NWIS_RDB = textwrap.dedent("""\
    # --------------------------------
    # USGS Water Data
    # site_no  peak_dt         peak_va peak_cd
    # -------- --------------- ------- -------
    agency_cd\tsite_no\tpeak_dt\tpeak_va\tpeak_cd
    5s\t15s\t10d\t10d\t10d
    USGS\t03606500\t1970-01-01\t25000\t
    USGS\t03606500\t1971-01-01\t18000\t
    USGS\t03606500\t1972-01-01\t32000\tE
    """)


class TestFetchNwisPeakData:
    def test_success_returns_dataframe(self):
        mock_resp = MagicMock()
        mock_resp.text = MOCK_NWIS_RDB
        mock_resp.raise_for_status = MagicMock()

        with patch("hydrolib.regression.roi_analysis.requests.get", return_value=mock_resp):
            df = fetch_nwis_peak_data("03606500")

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 3
        assert "peak_va" in df.columns
        assert all(df["peak_va"] > 0)

    def test_network_error_returns_empty(self):
        import requests as req

        with patch(
            "hydrolib.regression.roi_analysis.requests.get",
            side_effect=req.ConnectionError("timeout"),
        ):
            df = fetch_nwis_peak_data("99999999")

        assert isinstance(df, pd.DataFrame)
        assert df.empty

    def test_http_error_returns_empty(self):
        import requests as req

        mock_resp = MagicMock()
        mock_resp.raise_for_status.side_effect = req.HTTPError("404")
        with patch("hydrolib.regression.roi_analysis.requests.get", return_value=mock_resp):
            df = fetch_nwis_peak_data("99999999")

        assert df.empty

    def test_zero_flows_filtered(self):
        rdb = MOCK_NWIS_RDB.replace("25000", "0")
        mock_resp = MagicMock()
        mock_resp.text = rdb
        mock_resp.raise_for_status = MagicMock()

        with patch("hydrolib.regression.roi_analysis.requests.get", return_value=mock_resp):
            df = fetch_nwis_peak_data("03606500")

        assert all(df["peak_va"] > 0)


# ---------------------------------------------------------------------------
# RoiAnalysis — distance and weighting tests
# ---------------------------------------------------------------------------


class TestRoiAnalysisWeighting:
    def test_char_distance_identical(self, target):
        """Distance from target to itself should be 0."""
        d = RoiAnalysis._char_distance(target, target)
        assert d == pytest.approx(0.0)

    def test_char_distance_da_only(self):
        """Verify distance formula for DA-only case."""
        t = BasinCharacteristics(
            site_no="T",
            site_name="T",
            region=TN_AREA3,
            predictors={DRNAREA: 100.0},
        )
        c = BasinCharacteristics(
            site_no="C",
            site_name="C",
            region=TN_AREA3,
            predictors={DRNAREA: 200.0},
        )
        expected = abs(math.log10(200.0) - math.log10(100.0))
        assert RoiAnalysis._char_distance(t, c) == pytest.approx(expected, rel=1e-9)

    def test_char_distance_with_slope(self, target):
        site = make_roi_site("X", da=200.0, slope=10.0)
        d = RoiAnalysis._char_distance(target, site.basin)
        assert d > 0

    def test_geo_distance_zero_coords(self, target):
        """Missing coordinates should return 0 (equal weights default)."""
        no_coord = BasinCharacteristics(
            site_no="X",
            site_name="X",
            region=TN_AREA2,
            predictors={DRNAREA: 100.0, CSL1085LFP: 5.0},
        )
        assert RoiAnalysis._geo_distance(target, no_coord) == 0.0

    def test_geo_distance_known_pair(self):
        """Known distance: Nashville to Memphis ~ 320 km."""
        nashville = BasinCharacteristics(
            site_no="A",
            site_name="A",
            region=TN_AREA3,
            predictors={DRNAREA: 1.0},
            latitude=36.1627,
            longitude=-86.7816,
        )
        memphis = BasinCharacteristics(
            site_no="B",
            site_name="B",
            region=TN_AREA3,
            predictors={DRNAREA: 1.0},
            latitude=35.1495,
            longitude=-90.0490,
        )
        d = RoiAnalysis._geo_distance(nashville, memphis)
        # Haversine Nashville-Memphis ≈ 320 km
        assert 300 < d < 340

    def test_compute_weights_sums_to_one(self, roi):
        roi.compute_weights()
        active = [s for s in roi.candidate_sites if s.weight > 0]
        assert active, "No active sites after weighting"
        total = sum(s.weight for s in active)
        assert total == pytest.approx(1.0, rel=1e-9)

    def test_compute_weights_nearest_gets_most(self, roi):
        """The site closest to target in char space should get highest weight."""
        roi.compute_weights()
        # Site G: DA=120, slope=5.0 — closest to target DA=100, slope=5.0
        # Site B: DA=90, slope=5.1 — also close
        site_g = next(s for s in roi.candidate_sites if s.site_no == "G")
        site_d = next(s for s in roi.candidate_sites if s.site_no == "D")  # furthest
        assert site_g.weight > site_d.weight

    def test_not_enough_sites_raises(self, target):
        few_sites = [make_roi_site("X", 100.0, 5.0)]  # only 1 site
        roi_small = RoiAnalysis(target=target, candidate_sites=few_sites, min_sites=5)
        with pytest.raises(RuntimeError, match="need at least"):
            roi_small.compute_weights()

    def test_combined_weight_method(self, roi, target):
        """Combined method should produce positive weights summing to 1."""
        roi_comb = RoiAnalysis(
            target=target,
            candidate_sites=roi.candidate_sites,
            weight_method="combined",
            min_sites=3,
        )
        roi_comb.compute_weights()
        total = sum(s.weight for s in roi_comb.candidate_sites)
        assert total == pytest.approx(1.0, rel=1e-9)


# ---------------------------------------------------------------------------
# RoiAnalysis — quantile computation
# ---------------------------------------------------------------------------


class TestRoiQuantiles:
    def test_compute_roi_quantiles_returns_dict(self, roi):
        roi.compute_weights()
        qs = roi.compute_roi_quantiles()
        assert isinstance(qs, dict)
        assert 0.01 in qs
        assert qs[0.01] > 0

    def test_roi_quantile_monotonic_with_return_period(self, roi):
        """Higher return period should yield higher flow."""
        roi.compute_weights()
        qs = roi.compute_roi_quantiles()
        sorted_aeps = sorted(qs.keys(), reverse=True)  # high AEP = low return period
        flows = [qs[a] for a in sorted_aeps]
        assert all(flows[i] <= flows[i + 1] for i in range(len(flows) - 1))

    def test_roi_quantiles_no_weights_raises(self, roi):
        with pytest.raises(RuntimeError, match="call compute_weights"):
            roi.compute_roi_quantiles()

    def test_roi_variance_finite(self, roi):
        roi.compute_weights()
        variances = roi.roi_variance()
        for aep, v in variances.items():
            if not math.isnan(v):
                assert v >= 0.0

    def test_summary_string(self, roi):
        s = roi.summary()
        assert "UNGAGED" in s
        assert "Candidate sites" in s

    # ------------------------------------------------------------------
    # Integration: fetch → run_ema → compute_weights → quantiles
    # ------------------------------------------------------------------

    def test_fetch_and_ema_integration(self, roi):
        """
        Mock NWIS fetch and verify the full pipeline: fetch → EMA → weights → Q.
        Runs using real B17CEngine with synthetic data.
        """
        import numpy as np

        rng = np.random.default_rng(42)

        # Produce a synthetic RDB-style response for each site
        def make_rdb(site_no: str, n: int = 40) -> str:
            flows = rng.lognormal(mean=10, sigma=0.5, size=n)
            rows = "\n".join(
                f"USGS\t{site_no}\t{1970 + i}-01-01\t{f:.0f}\t" for i, f in enumerate(flows)
            )
            header = "agency_cd\tsite_no\tpeak_dt\tpeak_va\tpeak_cd"
            typeline = "5s\t15s\t10d\t10d\t10d"
            return f"# comment\n{header}\n{typeline}\n{rows}\n"

        def fake_get(url, params=None, timeout=30):
            mock = MagicMock()
            mock.text = make_rdb(params.get("site_no", "00000000"))
            mock.raise_for_status = MagicMock()
            return mock

        with patch("hydrolib.regression.roi_analysis.requests.get", side_effect=fake_get):
            counts = roi.fetch_nwis_peaks()

        assert all(n >= 40 for n in counts.values())
        assert all(s.record_source == "nwis" for s in roi.candidate_sites)

        errors = roi.run_ema(min_record_length=10)
        fitted = [s for s in roi.candidate_sites if s.has_lp3]
        assert len(fitted) >= 5

        roi.compute_weights()
        qs = roi.compute_roi_quantiles()
        assert qs[0.01] > qs[0.5]  # Q100 > Q2
