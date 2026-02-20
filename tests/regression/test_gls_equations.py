"""
Tests for the generic regression infrastructure:
  hydrolib.regression.region          — HydrologicRegion, RegionRegistry
  hydrolib.regression.basin_chars     — BasinCharacteristics
  hydrolib.regression.regression_table — GlsEquation, RegressionTable
  hydrolib.regression.sir2024_5130    — TN regions, HydrologicArea, SIR2024_5130

Synthetic (hypothetical) coefficients are used throughout so that
computational logic can be verified without the actual SIR 2024-5130 PDF.
Do NOT use these coefficients for engineering design.
"""

from __future__ import annotations

import csv
import math
import tempfile
from pathlib import Path
from typing import Dict, Tuple

import pytest

from hydrolib.regression.basin_chars import (
    CSL1085LFP,
    DRNAREA,
    ELEV,
    IMPERV,
    BasinCharacteristics,
)
from hydrolib.regression.region import HydrologicRegion, RegionRegistry
from hydrolib.regression.regression_table import (
    STANDARD_AEPS,
    GlsEquation,
    RegressionTable,
)
from hydrolib.regression.sir2024_5130 import (
    SIR2024_5130,
    TN_AREA1,
    TN_AREA2,
    TN_AREA3,
    TN_AREA4,
    TN_REGIONS,
    HydrologicArea,
)

# ---------------------------------------------------------------------------
# Shared region fixtures (one TN, one KY, one multi-predictor)
# ---------------------------------------------------------------------------

TN2 = TN_AREA2  # DA + slope
TN3 = TN_AREA3  # DA only

KY_REG1 = HydrologicRegion(
    code="KY_REGION1",
    label="Kentucky Flood Region 1",
    state="KY",
    required_predictors=(DRNAREA, CSL1085LFP),
    publication="WRIR 03-4180",
)

GA_BLUE_RIDGE = HydrologicRegion(
    code="GA_BLUE_RIDGE",
    label="Georgia Blue Ridge",
    state="GA",
    required_predictors=(DRNAREA, ELEV),
    publication="SIR 2017-5038",
)


def _make_tn2_basin(da: float = 100.0, slope: float = 5.0) -> BasinCharacteristics:
    return BasinCharacteristics(
        site_no="SITE01",
        site_name="Test Creek TN",
        region=TN2,
        predictors={DRNAREA: da, CSL1085LFP: slope},
    )


def _make_ky_basin(da: float = 300.0, slope: float = 3.0) -> BasinCharacteristics:
    return BasinCharacteristics(
        site_no="SITE_KY",
        site_name="Test Creek KY",
        region=KY_REG1,
        predictors={DRNAREA: da, CSL1085LFP: slope},
    )


def _make_ga_basin(da: float = 80.0, elev: float = 2500.0) -> BasinCharacteristics:
    return BasinCharacteristics(
        site_no="SITE_GA",
        site_name="Test Creek GA",
        region=GA_BLUE_RIDGE,
        predictors={DRNAREA: da, ELEV: elev},
    )


def _make_coeff_dict(
    regions_and_preds: list,
) -> Dict[Tuple[HydrologicRegion, float], dict]:
    """
    Build a synthetic coefficient dict for multiple regions across all 8 AEPs.

    regions_and_preds: list of (region, predictor_codes)
    """
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
    result = {}
    for region, pred_codes in regions_and_preds:
        for aep, b0 in b0_by_aep.items():
            coeffs = {code: 0.75 if i == 0 else 0.35 for i, code in enumerate(pred_codes)}
            result[(region, aep)] = {
                "intercept": b0,
                "coefficients": coeffs,
                "sep_pct": 36.0,
                "pseudo_r2": 0.92,
                "eyr": 18.0,
            }
    return result


# ---------------------------------------------------------------------------
# HydrologicRegion tests
# ---------------------------------------------------------------------------


class TestHydrologicRegion:
    def test_frozen_and_hashable(self):
        """Regions must be hashable (used as dict keys)."""
        s = {TN2, TN3, KY_REG1}
        assert len(s) == 3

    def test_required_predictors_tuple(self):
        assert isinstance(TN2.required_predictors, tuple)
        assert DRNAREA in TN2.required_predictors
        assert CSL1085LFP in TN2.required_predictors

    def test_validate_predictors_ok(self):
        TN2.validate_predictors({DRNAREA: 100.0, CSL1085LFP: 5.0})  # should not raise

    def test_validate_predictors_missing(self):
        with pytest.raises(ValueError, match="missing predictors"):
            TN2.validate_predictors({DRNAREA: 100.0})  # missing slope

    def test_validate_predictors_nonpositive(self):
        with pytest.raises(ValueError, match="non-positive"):
            TN2.validate_predictors({DRNAREA: 100.0, CSL1085LFP: 0.0})

    def test_has_predictor(self):
        assert TN2.has_predictor(CSL1085LFP)
        assert not TN2.has_predictor(IMPERV)

    def test_from_dict_round_trip(self):
        d = TN2.to_dict()
        r2 = HydrologicRegion.from_dict(d)
        assert r2.code == TN2.code
        assert r2.state == TN2.state
        assert set(r2.required_predictors) == set(TN2.required_predictors)

    def test_str(self):
        assert "TN_AREA2" in str(TN2)

    def test_area3_da_only(self):
        assert TN_AREA3.required_predictors == (DRNAREA,)

    def test_tn_regions_tuple(self):
        assert len(TN_REGIONS) == 4
        assert TN_AREA1 in TN_REGIONS


class TestRegionRegistry:
    def test_register_and_get(self):
        reg = RegionRegistry()
        reg.register(KY_REG1)
        assert reg["KY_REGION1"] is KY_REG1

    def test_contains(self):
        reg = RegionRegistry()
        reg.register(TN2)
        assert "TN_AREA2" in reg
        assert "TN_AREA1" not in reg

    def test_missing_key_raises(self):
        reg = RegionRegistry()
        with pytest.raises(KeyError):
            _ = reg["NONEXISTENT"]

    def test_all_regions(self):
        reg = RegionRegistry()
        for r in [TN2, TN3, KY_REG1]:
            reg.register(r)
        assert len(reg.all_regions()) == 3

    def test_load_from_csv(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, newline="") as fh:
            tmp = fh.name
            writer = csv.DictWriter(
                fh, fieldnames=["code", "label", "state", "required_predictors", "publication"]
            )
            writer.writeheader()
            writer.writerow(
                {
                    "code": "OR_COAST",
                    "label": "Oregon Coast Range",
                    "state": "OR",
                    "required_predictors": "DRNAREA, PRECIP",
                    "publication": "SIR 2014-5048",
                }
            )

        reg = RegionRegistry.load_from_csv(tmp)
        Path(tmp).unlink()
        assert "OR_COAST" in reg
        r = reg["OR_COAST"]
        assert "DRNAREA" in r.required_predictors
        assert "PRECIP" in r.required_predictors


# ---------------------------------------------------------------------------
# BasinCharacteristics tests
# ---------------------------------------------------------------------------


class TestBasinCharacteristics:
    def test_valid_tn2(self):
        b = _make_tn2_basin()
        assert b.drainage_area_sqmi == 100.0
        assert b.slope_1085_ftmi == 5.0
        assert b.region is TN2

    def test_valid_ky(self):
        b = _make_ky_basin()
        assert b.drainage_area_sqmi == 300.0

    def test_valid_ga_with_elev(self):
        b = _make_ga_basin()
        assert b.predictors[ELEV] == 2500.0

    def test_missing_required_raises(self):
        with pytest.raises(ValueError, match="missing predictors"):
            BasinCharacteristics(
                site_no="X",
                site_name="X",
                region=TN2,
                predictors={DRNAREA: 100.0},  # missing CSL1085LFP
            )

    def test_zero_predictor_raises(self):
        with pytest.raises(ValueError, match="non-positive"):
            BasinCharacteristics(
                site_no="X",
                site_name="X",
                region=TN3,
                predictors={DRNAREA: 0.0},
            )

    def test_convenience_properties_none_when_absent(self):
        b = _make_tn2_basin()
        assert b.climate_factor_2yr is None  # I2 not in predictors
        assert b.pct_impervious is None

    def test_state_from_region(self):
        b = _make_tn2_basin()
        assert b.state == "TN"
        b_ky = _make_ky_basin()
        assert b_ky.state == "KY"

    def test_predictor_value_ok(self):
        b = _make_tn2_basin(da=250.0, slope=3.5)
        assert b.predictor_value(DRNAREA) == 250.0

    def test_predictor_value_missing_raises(self):
        b = _make_tn2_basin()
        with pytest.raises(KeyError):
            b.predictor_value(ELEV)

    def test_extra_predictors_stored(self):
        """Extra predictors beyond required_predictors are preserved."""
        b = BasinCharacteristics(
            site_no="X",
            site_name="X",
            region=TN2,
            predictors={DRNAREA: 100.0, CSL1085LFP: 5.0, ELEV: 1200.0},
        )
        assert b.predictors[ELEV] == 1200.0

    def test_from_streamstats_valid(self):
        params = {DRNAREA: 85.3, CSL1085LFP: 4.7, ELEV: 980.0}
        b = BasinCharacteristics.from_streamstats(
            "SS01", "SS Site", TN2, params, latitude=36.0, longitude=-86.5
        )
        assert b.drainage_area_sqmi == 85.3
        assert b.slope_1085_ftmi == 4.7
        assert b.latitude == 36.0

    def test_from_streamstats_skips_nonnumeric(self):
        params = {DRNAREA: 100.0, CSL1085LFP: 5.0, "DESCRIPTION": "some text"}
        b = BasinCharacteristics.from_streamstats("X", "X", TN2, params)
        assert "DESCRIPTION" not in b.predictors

    def test_from_streamstats_missing_required_raises(self):
        with pytest.raises(ValueError):
            BasinCharacteristics.from_streamstats("X", "X", TN2, {DRNAREA: 100.0})

    def test_summary_contains_site_and_region(self):
        s = _make_tn2_basin().summary()
        assert "SITE01" in s
        assert "TN_AREA2" in s


# ---------------------------------------------------------------------------
# GlsEquation tests
# ---------------------------------------------------------------------------


class TestGlsEquation:
    def test_compute_two_predictors(self):
        eq = GlsEquation(
            region=TN2,
            aep=0.01,
            intercept=2.90,
            coefficients={DRNAREA: 0.75, CSL1085LFP: 0.38},
        )
        b = _make_tn2_basin(da=779.0, slope=2.4)
        expected_log = 2.90 + 0.75 * math.log10(779.0) + 0.38 * math.log10(2.4)
        assert math.isclose(math.log10(eq.compute(b)), expected_log, rel_tol=1e-9)

    def test_compute_da_only(self):
        eq = GlsEquation(
            region=TN3,
            aep=0.01,
            intercept=2.85,
            coefficients={DRNAREA: 0.72},
        )
        b = BasinCharacteristics(
            site_no="X",
            site_name="X",
            region=TN3,
            predictors={DRNAREA: 120.5},
        )
        expected_log = 2.85 + 0.72 * math.log10(120.5)
        assert math.isclose(math.log10(eq.compute(b)), expected_log, rel_tol=1e-9)

    def test_compute_three_predictors(self):
        """Arbitrary predictor set (GA Blue Ridge: DA + ELEV)."""
        eq = GlsEquation(
            region=GA_BLUE_RIDGE,
            aep=0.01,
            intercept=1.50,
            coefficients={DRNAREA: 0.80, ELEV: 0.25},
        )
        b = _make_ga_basin(da=80.0, elev=2500.0)
        expected_log = 1.50 + 0.80 * math.log10(80.0) + 0.25 * math.log10(2500.0)
        assert math.isclose(math.log10(eq.compute(b)), expected_log, rel_tol=1e-9)

    def test_wrong_region_raises(self):
        eq = GlsEquation(region=TN3, aep=0.01, intercept=2.9, coefficients={DRNAREA: 0.75})
        with pytest.raises(ValueError, match="region"):
            eq.compute(_make_tn2_basin())  # TN2 basin, TN3 equation

    def test_missing_predictor_raises(self):
        eq = GlsEquation(
            region=TN2,
            aep=0.01,
            intercept=2.90,
            coefficients={DRNAREA: 0.75, CSL1085LFP: 0.38},
        )
        # Build a basin that bypasses validation (empty predictors would fail)
        # Instead, force a missing predictor at compute time
        b = BasinCharacteristics(
            site_no="X",
            site_name="X",
            region=TN2,
            predictors={DRNAREA: 100.0, CSL1085LFP: 5.0},
        )
        # Monkey-patch to remove slope after construction
        b.predictors.pop(CSL1085LFP)
        with pytest.raises(KeyError):
            eq.compute(b)

    def test_return_period(self):
        eq = GlsEquation(region=TN3, aep=0.01, intercept=2.0, coefficients={})
        assert math.isclose(eq.return_period, 100.0)

    def test_variance_log10_from_sep(self):
        eq = GlsEquation(region=TN3, aep=0.01, intercept=2.0, coefficients={}, sep_pct=40.0)
        expected = (math.log10(1.4)) ** 2
        assert math.isclose(eq.variance_log10, expected, rel_tol=1e-9)

    def test_variance_none_without_sep(self):
        eq = GlsEquation(region=TN3, aep=0.01, intercept=2.0, coefficients={})
        assert eq.variance_log10 is None

    def test_predictor_codes_sorted(self):
        eq = GlsEquation(
            region=TN2, aep=0.01, intercept=2.0, coefficients={CSL1085LFP: 0.38, DRNAREA: 0.75}
        )
        assert eq.predictor_codes == sorted([DRNAREA, CSL1085LFP])

    def test_to_dict_has_region_and_aep(self):
        eq = GlsEquation(
            region=TN2, aep=0.01, intercept=2.90, coefficients={DRNAREA: 0.75, CSL1085LFP: 0.38}
        )
        d = eq.to_dict()
        assert d["region_code"] == "TN_AREA2"
        assert d["aep"] == 0.01
        assert d["intercept"] == 2.90
        assert d[DRNAREA] == 0.75


# ---------------------------------------------------------------------------
# RegressionTable — generic tests (multi-state)
# ---------------------------------------------------------------------------


class TestRegressionTableGeneric:
    @pytest.fixture
    def multi_state_table(self) -> RegressionTable:
        """Table spanning TN (Area2, Area3), KY, GA regions."""
        coeff_dict = _make_coeff_dict(
            [
                (TN2, [DRNAREA, CSL1085LFP]),
                (TN3, [DRNAREA]),
                (KY_REG1, [DRNAREA, CSL1085LFP]),
                (GA_BLUE_RIDGE, [DRNAREA, ELEV]),
            ]
        )
        return RegressionTable.load_from_dict(coeff_dict)

    def test_available_states(self, multi_state_table):
        states = multi_state_table.available_states()
        assert "TN" in states
        assert "KY" in states
        assert "GA" in states

    def test_available_regions(self, multi_state_table):
        regions = multi_state_table.available_regions()
        codes = {r.code for r in regions}
        assert "TN_AREA2" in codes
        assert "KY_REGION1" in codes
        assert "GA_BLUE_RIDGE" in codes

    def test_estimate_tn2(self, multi_state_table):
        b = _make_tn2_basin()
        q = multi_state_table.estimate(b, aep=0.01)
        assert q > 0

    def test_estimate_ky(self, multi_state_table):
        b = _make_ky_basin()
        q = multi_state_table.estimate(b, aep=0.01)
        assert q > 0

    def test_estimate_ga_with_elev(self, multi_state_table):
        b = _make_ga_basin()
        q = multi_state_table.estimate(b, aep=0.01)
        assert q > 0

    def test_estimate_all_aeps_count(self, multi_state_table):
        b = _make_tn2_basin()
        results = multi_state_table.estimate_all_aeps(b)
        assert len(results) == 8
        assert all(q > 0 for q in results.values())

    def test_missing_region_raises(self, multi_state_table):
        unknown = HydrologicRegion("ZZ_X", "Unknown", "ZZ", required_predictors=(DRNAREA,))
        b = BasinCharacteristics(
            site_no="X", site_name="X", region=unknown, predictors={DRNAREA: 100.0}
        )
        with pytest.raises(KeyError):
            multi_state_table.estimate(b, aep=0.01)

    def test_missing_aep_raises(self, multi_state_table):
        b = _make_tn2_basin()
        with pytest.raises(KeyError):
            multi_state_table.get_equation(TN2, 0.333)

    def test_repr(self, multi_state_table):
        r = repr(multi_state_table)
        assert "RegressionTable" in r
        assert "TN" in r

    def test_len(self, multi_state_table):
        # 4 regions × 8 AEPs = 32
        assert len(multi_state_table) == 32

    # ------------------------------------------------------------------
    # CSV round-trip
    # ------------------------------------------------------------------

    def test_csv_round_trip(self, multi_state_table):
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as fh:
            tmp = fh.name

        multi_state_table.to_csv(tmp)
        loaded = RegressionTable.load_from_csv(tmp)
        Path(tmp).unlink()

        assert len(loaded) == len(multi_state_table)
        eq_orig = multi_state_table.get_equation(TN2, 0.01)
        eq_load = loaded.get_by_region_code("TN_AREA2", 0.01)
        assert math.isclose(eq_orig.intercept, eq_load.intercept)
        assert math.isclose(eq_orig.coefficients[DRNAREA], eq_load.coefficients[DRNAREA])

    def test_csv_round_trip_multi_predictor(self, multi_state_table):
        """GA equations use ELEV not CSL1085LFP — verify no cross-contamination."""
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as fh:
            tmp = fh.name

        multi_state_table.to_csv(tmp)
        loaded = RegressionTable.load_from_csv(tmp)
        Path(tmp).unlink()

        eq = loaded.get_by_region_code("GA_BLUE_RIDGE", 0.01)
        assert ELEV in eq.coefficients
        # KY equation should not have ELEV
        eq_ky = loaded.get_by_region_code("KY_REGION1", 0.01)
        assert ELEV not in eq_ky.coefficients

    def test_load_from_csv_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            RegressionTable.load_from_csv("/nonexistent/equations.csv")

    # ------------------------------------------------------------------
    # write_template helper
    # ------------------------------------------------------------------

    def test_write_template_row_count(self):
        regions = [TN2, TN3, KY_REG1]
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as fh:
            tmp = fh.name
        RegressionTable.write_template(tmp, regions=regions)
        content = Path(tmp).read_text()
        Path(tmp).unlink()
        lines = [l for l in content.splitlines() if l.strip()]
        # header + 3 regions × 8 AEPs = 25
        assert len(lines) == 25

    def test_write_template_has_predictor_columns(self):
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as fh:
            tmp = fh.name
        RegressionTable.write_template(tmp, regions=[TN2, GA_BLUE_RIDGE])
        header = Path(tmp).read_text().splitlines()[0]
        Path(tmp).unlink()
        assert DRNAREA in header
        assert CSL1085LFP in header
        assert ELEV in header


# ---------------------------------------------------------------------------
# SIR2024_5130 and HydrologicArea (backward-compat) tests
# ---------------------------------------------------------------------------


class TestHydrologicArea:
    def test_class_attrs_are_regions(self):
        assert isinstance(HydrologicArea.AREA1, HydrologicRegion)
        assert isinstance(HydrologicArea.AREA2, HydrologicRegion)

    def test_area2_code(self):
        assert HydrologicArea.AREA2.code == "TN_AREA2"

    def test_from_int_valid(self):
        assert HydrologicArea.from_int(1) is HydrologicArea.AREA1
        assert HydrologicArea.from_int(4) is HydrologicArea.AREA4

    def test_from_int_invalid(self):
        with pytest.raises(ValueError, match="1-4"):
            HydrologicArea.from_int(5)

    def test_all_regions(self):
        regions = HydrologicArea.all_regions()
        assert len(regions) == 4
        assert HydrologicArea.AREA2 in regions


class TestSIR2024_5130:
    @pytest.fixture
    def legacy_coeff_dict(self) -> dict:
        """Legacy b0/b1/b2 format coefficient dict for all four TN areas."""
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
        areas = {
            HydrologicArea.AREA1: dict(b1=0.75, b2=0.45),
            HydrologicArea.AREA2: dict(b1=0.75, b2=0.38),
            HydrologicArea.AREA3: dict(b1=0.72, b2=None),
            HydrologicArea.AREA4: dict(b1=0.70, b2=0.12),
        }
        result = {}
        for area, kw in areas.items():
            for aep, b0 in b0_by_aep.items():
                result[(area, aep)] = {
                    "b0": b0,
                    "b1": kw["b1"],
                    "b2": kw["b2"],
                    "sep_pct": 36.0,
                    "pseudo_r2": 0.92,
                    "eyr": 18.0,
                }
        return result

    @pytest.fixture
    def gls_table(self, legacy_coeff_dict) -> SIR2024_5130:
        return SIR2024_5130.load_from_dict(legacy_coeff_dict)

    def test_is_regression_table(self, gls_table):
        assert isinstance(gls_table, RegressionTable)

    def test_load_from_dict_legacy_count(self, gls_table):
        assert len(gls_table) == 32

    def test_legacy_b2_maps_to_correct_predictor(self, gls_table):
        """b2 for Area2 must map to CSL1085LFP (the 2nd required_predictor)."""
        eq = gls_table.get_equation(HydrologicArea.AREA2, 0.01)
        assert CSL1085LFP in eq.coefficients
        assert math.isclose(eq.coefficients[CSL1085LFP], 0.38)

    def test_legacy_area3_no_b2(self, gls_table):
        """Area 3 (DA-only) must not have a second predictor coefficient."""
        eq = gls_table.get_equation(HydrologicArea.AREA3, 0.01)
        assert CSL1085LFP not in eq.coefficients
        assert len(eq.coefficients) == 1

    def test_estimate_area2(self, gls_table):
        basin = BasinCharacteristics(
            site_no="03606500",
            site_name="Big Sandy River",
            region=HydrologicArea.AREA2,
            predictors={DRNAREA: 779.0, CSL1085LFP: 2.4},
        )
        q = gls_table.estimate(basin, aep=0.01)
        expected_log = 2.90 + 0.75 * math.log10(779.0) + 0.38 * math.log10(2.4)
        assert math.isclose(math.log10(q), expected_log, rel_tol=1e-6)

    def test_estimate_area3_da_only(self, gls_table):
        basin = BasinCharacteristics(
            site_no="X",
            site_name="X",
            region=HydrologicArea.AREA3,
            predictors={DRNAREA: 120.5},
        )
        q = gls_table.estimate(basin, aep=0.01)
        expected_log = 2.90 + 0.72 * math.log10(120.5)
        assert math.isclose(math.log10(q), expected_log, rel_tol=1e-6)

    def test_get_variance_area2(self, gls_table):
        v = gls_table.get_variance(HydrologicArea.AREA2, 0.01)
        assert v is not None and v > 0

    def test_estimate_all_aeps_count(self, gls_table):
        basin = BasinCharacteristics(
            site_no="X",
            site_name="X",
            region=HydrologicArea.AREA2,
            predictors={DRNAREA: 100.0, CSL1085LFP: 5.0},
        )
        results = gls_table.estimate_all_aeps(basin)
        assert len(results) == 8

    def test_publication_set(self, gls_table):
        assert "2024-5130" in gls_table.publication or "sir20245130" in gls_table.publication

    def test_write_table4_template(self):
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as fh:
            tmp = fh.name
        SIR2024_5130.write_table4_template(tmp)
        content = Path(tmp).read_text()
        Path(tmp).unlink()
        lines = [l for l in content.splitlines() if l.strip()]
        # header + 4 areas × 8 AEPs = 33
        assert len(lines) == 33
        assert DRNAREA in content

    def test_csv_round_trip_via_sir_subclass(self, gls_table):
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as fh:
            tmp = fh.name
        gls_table.to_csv(tmp)
        loaded = SIR2024_5130.load_from_csv(tmp)
        Path(tmp).unlink()

        assert len(loaded) == 32
        eq = loaded.get_equation(HydrologicArea.AREA2, 0.01)
        assert math.isclose(eq.intercept, 2.90)
        assert math.isclose(eq.coefficients[DRNAREA], 0.75)
