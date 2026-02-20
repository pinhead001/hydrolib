"""
Tests for hydrolib.regression.basin_chars and hydrolib.regression.sir2024_5130.

Uses synthetic (hypothetical) GLS coefficients to test the computational
logic independently of the specific Table 4 values from SIR 2024-5130.
"""

from __future__ import annotations

import csv
import math
import tempfile
from pathlib import Path

import pytest

from hydrolib.regression.basin_chars import BasinCharacteristics, HydrologicArea
from hydrolib.regression.sir2024_5130 import (
    SIR2024_5130,
    STANDARD_AEPS,
    GlsEquation,
    write_table4_template,
)

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def area2_basin() -> BasinCharacteristics:
    """Area 2 basin (DA + S1085 predictors)."""
    return BasinCharacteristics(
        site_no="03606500",
        site_name="Big Sandy River at Bruceton, TN",
        drainage_area_sqmi=779.0,
        hydrologic_area=HydrologicArea.AREA2,
        slope_1085_ftmi=2.4,
    )


@pytest.fixture
def area3_basin() -> BasinCharacteristics:
    """Area 3 basin (DA-only predictor)."""
    return BasinCharacteristics(
        site_no="03430100",
        site_name="Example Creek, TN",
        drainage_area_sqmi=120.5,
        hydrologic_area=HydrologicArea.AREA3,
    )


@pytest.fixture
def area1_basin() -> BasinCharacteristics:
    """Area 1 basin (DA + 2-yr climate factor)."""
    return BasinCharacteristics(
        site_no="07030392",
        site_name="Wolf River near Memphis, TN",
        drainage_area_sqmi=888.0,
        hydrologic_area=HydrologicArea.AREA1,
        climate_factor_2yr=2.15,
    )


@pytest.fixture
def area4_basin() -> BasinCharacteristics:
    """Area 4 basin (DA + % impervious)."""
    return BasinCharacteristics(
        site_no="03431700",
        site_name="Urban Creek, TN",
        drainage_area_sqmi=22.0,
        hydrologic_area=HydrologicArea.AREA4,
        pct_impervious=18.0,
    )


@pytest.fixture
def synthetic_coeff_dict() -> dict:
    """
    Synthetic coefficient table for all 4 areas and 8 standard AEPs.
    These are HYPOTHETICAL values for testing purposes only and must NOT
    be used for engineering design.  Actual coefficients must be sourced
    from SIR 2024-5130 Table 4.
    """
    areas = {
        HydrologicArea.AREA1: dict(b2_label="I2", b2=0.45),
        HydrologicArea.AREA2: dict(b2_label="S1085", b2=0.38),
        HydrologicArea.AREA3: dict(b2_label=None, b2=None),
        HydrologicArea.AREA4: dict(b2_label="Imperv", b2=0.12),
    }
    # Intercept increases with return period; DA exponent ~ 0.75
    aep_b0 = {
        0.5: 1.80,
        0.2: 2.10,
        0.1: 2.30,
        0.04: 2.55,
        0.02: 2.70,
        0.01: 2.90,
        0.005: 3.05,
        0.002: 3.25,
    }
    aep_sep = {
        0.5: 44.0,
        0.2: 40.0,
        0.1: 38.0,
        0.04: 36.0,
        0.02: 35.0,
        0.01: 34.0,
        0.005: 35.0,
        0.002: 38.0,
    }
    result = {}
    for area, meta in areas.items():
        for aep, b0 in aep_b0.items():
            result[(area, aep)] = {
                "b0": b0,
                "b1": 0.75,
                "b2": meta["b2"],
                "sep_pct": aep_sep[aep],
                "pseudo_r2": 0.92,
                "eyr": 18.0,
            }
    return result


@pytest.fixture
def gls_table(synthetic_coeff_dict) -> SIR2024_5130:
    return SIR2024_5130.load_from_dict(synthetic_coeff_dict)


# ---------------------------------------------------------------------------
# HydrologicArea tests
# ---------------------------------------------------------------------------


class TestHydrologicArea:
    def test_from_int_valid(self):
        assert HydrologicArea.from_int(1) == HydrologicArea.AREA1
        assert HydrologicArea.from_int(4) == HydrologicArea.AREA4

    def test_from_int_invalid(self):
        with pytest.raises(ValueError, match="must be 1-4"):
            HydrologicArea.from_int(5)

    def test_label_format(self):
        assert "1" in HydrologicArea.AREA1.label or "A" in HydrologicArea.AREA1.label

    def test_value_string(self):
        assert HydrologicArea.AREA2.value == "area2"


# ---------------------------------------------------------------------------
# BasinCharacteristics tests
# ---------------------------------------------------------------------------


class TestBasinCharacteristics:
    def test_valid_area2(self, area2_basin):
        assert area2_basin.drainage_area_sqmi == 779.0
        assert area2_basin.slope_1085_ftmi == 2.4
        assert area2_basin.hydrologic_area == HydrologicArea.AREA2

    def test_area2_requires_slope(self):
        """Area 2 must have slope; should raise ValueError."""
        with pytest.raises(ValueError, match="slope_1085_ftmi is required"):
            BasinCharacteristics(
                site_no="X",
                site_name="X",
                drainage_area_sqmi=50.0,
                hydrologic_area=HydrologicArea.AREA2,
                slope_1085_ftmi=None,
            )

    def test_negative_da_raises(self):
        with pytest.raises(ValueError, match="must be > 0"):
            BasinCharacteristics(
                site_no="X",
                site_name="X",
                drainage_area_sqmi=-10.0,
                hydrologic_area=HydrologicArea.AREA3,
            )

    def test_zero_da_raises(self):
        with pytest.raises(ValueError):
            BasinCharacteristics(
                site_no="X",
                site_name="X",
                drainage_area_sqmi=0.0,
                hydrologic_area=HydrologicArea.AREA3,
            )

    def test_area1_warns_without_cf(self, caplog):
        """Area 1 without climate_factor_2yr should log a warning."""
        import logging

        with caplog.at_level(logging.WARNING, logger="hydrolib.regression.basin_chars"):
            BasinCharacteristics(
                site_no="X",
                site_name="X",
                drainage_area_sqmi=100.0,
                hydrologic_area=HydrologicArea.AREA1,
            )
        assert any("climate_factor_2yr" in r.message for r in caplog.records)

    def test_from_streamstats_valid(self):
        params = {
            "DRNAREA": 85.3,
            "CSL1085LFP": 4.7,
            "LAT": 35.5,
            "LNG": -86.2,
        }
        basin = BasinCharacteristics.from_streamstats(
            "SS01", "StreamStats Site", HydrologicArea.AREA2, params
        )
        assert basin.drainage_area_sqmi == 85.3
        assert basin.slope_1085_ftmi == 4.7
        assert basin.latitude == 35.5

    def test_from_streamstats_missing_da(self):
        with pytest.raises(KeyError, match="DRNAREA"):
            BasinCharacteristics.from_streamstats(
                "X", "X", HydrologicArea.AREA3, {"CSL1085LFP": 3.0}
            )

    def test_summary_contains_site_no(self, area2_basin):
        s = area2_basin.summary()
        assert "03606500" in s
        assert "779" in s


# ---------------------------------------------------------------------------
# GlsEquation tests
# ---------------------------------------------------------------------------


class TestGlsEquation:
    def test_area2_compute(self, area2_basin):
        """Area 2 equation: Q = 10^(b0 + b1*log(DA) + b2*log(S1085))."""
        eq = GlsEquation(
            hydrologic_area=HydrologicArea.AREA2,
            aep=0.01,
            b0=2.90,
            b1=0.75,
            b2=0.38,
            sep_pct=34.0,
            pseudo_r2=0.92,
        )
        expected_log = 2.90 + 0.75 * math.log10(779.0) + 0.38 * math.log10(2.4)
        assert math.isclose(math.log10(eq.compute(area2_basin)), expected_log, rel_tol=1e-9)

    def test_area3_compute_no_slope(self, area3_basin):
        """Area 3 equation (DA only): b2=None."""
        eq = GlsEquation(
            hydrologic_area=HydrologicArea.AREA3,
            aep=0.01,
            b0=2.90,
            b1=0.75,
        )
        expected_log = 2.90 + 0.75 * math.log10(120.5)
        assert math.isclose(math.log10(eq.compute(area3_basin)), expected_log, rel_tol=1e-9)

    def test_wrong_area_raises(self, area2_basin):
        """Equation for AREA3 applied to AREA2 basin should raise ValueError."""
        eq = GlsEquation(
            hydrologic_area=HydrologicArea.AREA3,
            aep=0.01,
            b0=2.90,
            b1=0.75,
        )
        with pytest.raises(ValueError, match="Hydrologic Area"):
            eq.compute(area2_basin)

    def test_area2_missing_slope_raises(self):
        """Area 2 equation applied to basin without slope raises ValueError."""
        basin = BasinCharacteristics(
            site_no="X",
            site_name="X",
            drainage_area_sqmi=100.0,
            hydrologic_area=HydrologicArea.AREA2,
            slope_1085_ftmi=2.0,  # bypass __post_init__
        )
        # Manually remove slope after construction
        object.__setattr__(basin, "slope_1085_ftmi", None)
        eq = GlsEquation(hydrologic_area=HydrologicArea.AREA2, aep=0.01, b0=2.9, b1=0.75, b2=0.38)
        with pytest.raises(ValueError, match="slope_1085_ftmi is required"):
            eq.compute(basin)

    def test_return_period(self):
        eq = GlsEquation(hydrologic_area=HydrologicArea.AREA3, aep=0.01, b0=2.0, b1=0.7)
        assert math.isclose(eq.return_period, 100.0)

    def test_variance_log10_from_sep(self):
        """Variance should be (log10(1 + SEP/100))^2."""
        eq = GlsEquation(
            hydrologic_area=HydrologicArea.AREA3, aep=0.01, b0=2.0, b1=0.7, sep_pct=40.0
        )
        expected = (math.log10(1.4)) ** 2
        assert math.isclose(eq.variance_log10, expected, rel_tol=1e-9)

    def test_variance_none_without_sep(self):
        eq = GlsEquation(hydrologic_area=HydrologicArea.AREA3, aep=0.01, b0=2.0, b1=0.7)
        assert eq.variance_log10 is None


# ---------------------------------------------------------------------------
# SIR2024_5130 table tests
# ---------------------------------------------------------------------------


class TestSIR2024_5130:
    def test_load_from_dict_count(self, gls_table):
        # 4 areas × 8 AEPs = 32 equations
        assert len(gls_table) == 32

    def test_estimate_area2(self, gls_table, area2_basin):
        q = gls_table.estimate(area2_basin, aep=0.01)
        assert q > 0

    def test_estimate_area3_da_only(self, gls_table, area3_basin):
        q = gls_table.estimate(area3_basin, aep=0.01)
        expected_log = 2.90 + 0.75 * math.log10(120.5)  # b2=None for area3
        assert math.isclose(math.log10(q), expected_log, rel_tol=1e-6)

    def test_estimate_all_aeps_count(self, gls_table, area2_basin):
        results = gls_table.estimate_all_aeps(area2_basin)
        assert len(results) == 8
        assert all(q > 0 for q in results.values())

    def test_missing_aep_raises(self, gls_table, area2_basin):
        with pytest.raises(KeyError):
            gls_table.get_equation(HydrologicArea.AREA2, 0.333)

    def test_available_aeps(self, gls_table):
        aeps = gls_table.available_aeps(HydrologicArea.AREA2)
        assert 0.01 in aeps
        assert 0.5 in aeps
        assert len(aeps) == 8

    def test_available_areas(self, gls_table):
        areas = gls_table.available_areas()
        assert HydrologicArea.AREA2 in areas
        assert len(areas) == 4

    def test_get_variance_area2(self, gls_table):
        v = gls_table.get_variance(HydrologicArea.AREA2, 0.01)
        assert v is not None
        assert v > 0

    def test_summary_table_structure(self, gls_table, area2_basin):
        rows = gls_table.summary_table(area2_basin)
        assert len(rows) == 8
        for row in rows:
            assert "aep" in row
            assert "flow_cfs" in row
            assert row["flow_cfs"] > 0

    def test_repr(self, gls_table):
        r = repr(gls_table)
        assert "SIR2024_5130" in r

    # ------------------------------------------------------------------
    # CSV round-trip
    # ------------------------------------------------------------------

    def test_load_from_csv_roundtrip(self, synthetic_coeff_dict):
        """Write a CSV and reload it; check a few values match."""
        import io

        headers = ["hydrologic_area", "aep", "b0", "b1", "b2", "sep_pct", "pseudo_r2", "eyr"]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, newline="") as fh:
            tmp_path = fh.name
            writer = csv.DictWriter(fh, fieldnames=headers)
            writer.writeheader()
            for (area, aep), kw in synthetic_coeff_dict.items():
                writer.writerow(
                    {
                        "hydrologic_area": area.value,
                        "aep": aep,
                        "b0": kw["b0"],
                        "b1": kw["b1"],
                        "b2": "" if kw.get("b2") is None else kw["b2"],
                        "sep_pct": kw.get("sep_pct", ""),
                        "pseudo_r2": kw.get("pseudo_r2", ""),
                        "eyr": kw.get("eyr", ""),
                    }
                )

        loaded = SIR2024_5130.load_from_csv(tmp_path)
        assert len(loaded) == 32
        eq = loaded.get_equation(HydrologicArea.AREA2, 0.01)
        assert math.isclose(eq.b0, synthetic_coeff_dict[(HydrologicArea.AREA2, 0.01)]["b0"])
        assert math.isclose(eq.b1, 0.75)

        Path(tmp_path).unlink()

    def test_load_from_csv_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            SIR2024_5130.load_from_csv("/nonexistent/path.csv")

    # ------------------------------------------------------------------
    # Template helper
    # ------------------------------------------------------------------

    def test_write_table4_template(self):
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as fh:
            tmp_path = fh.name

        write_table4_template(tmp_path)
        content = Path(tmp_path).read_text()
        assert "hydrologic_area" in content
        assert "aep" in content
        assert "b0" in content
        # Should have 4 areas × 8 AEPs = 32 data rows + header
        lines = [l for l in content.splitlines() if l.strip()]
        assert len(lines) == 33

        Path(tmp_path).unlink()
