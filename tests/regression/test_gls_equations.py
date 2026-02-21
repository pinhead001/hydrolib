"""
Tests for the generic regression framework and the TN / GA / MT state modules.

Coverage
--------
* HydrologicRegion — frozen, hashable, validate_predictors
* RegionRegistry — register, lookup, load_from_csv
* BasinCharacteristics — predictors dict API, state property, validation
* GlsEquation — compute(), variance_log10, wrong-region guard
* RegressionTable — load_from_dict, load_from_csv, batch_estimate, merge,
  filter_by_state, regions_for_state
* states.tennessee — build_tennessee_table(), region constants
* states.georgia   — build_georgia_table(), region constants
* states.montana   — build_montana_table(), region constants
* Nationwide merge — three-state table, batch_estimate across states
"""

from __future__ import annotations

import math

import pytest

from hydrolib.regression.basin_chars import (
    CSL1085LFP,
    DRNAREA,
    ELEV,
    FOREST,
    I2,
    IMPERV,
    PRECIP,
    BasinCharacteristics,
)
from hydrolib.regression.region import HydrologicRegion, RegionRegistry
from hydrolib.regression.regression_table import (
    STANDARD_AEPS,
    GlsEquation,
    RegressionTable,
)
from hydrolib.regression.states.georgia import (
    GA_BLUE_RIDGE,
    GA_COASTAL,
    GA_PIEDMONT,
    GA_REGIONS,
    GA_VALLEY_RIDGE,
    build_georgia_table,
)
from hydrolib.regression.states.montana import (
    MT_FOOTHILLS,
    MT_MOUNTAIN,
    MT_PLAINS,
    MT_REGIONS,
    build_montana_table,
)
from hydrolib.regression.states.tennessee import (
    TN_AREA1,
    TN_AREA2,
    TN_AREA3,
    TN_AREA4,
    TN_REGIONS,
    build_tennessee_table,
)

# ---------------------------------------------------------------------------
# Module-scoped fixtures (built once per test session for speed)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def tn_table() -> RegressionTable:
    return build_tennessee_table()


@pytest.fixture(scope="module")
def ga_table() -> RegressionTable:
    return build_georgia_table()


@pytest.fixture(scope="module")
def mt_table() -> RegressionTable:
    return build_montana_table()


@pytest.fixture(scope="module")
def national_table(tn_table, ga_table, mt_table) -> RegressionTable:
    return RegressionTable.merge(tn_table, ga_table, mt_table)


@pytest.fixture
def tn2_basin() -> BasinCharacteristics:
    return BasinCharacteristics(
        site_no="03606500",
        site_name="Big Sandy River at Bruceton, TN",
        region=TN_AREA2,
        predictors={DRNAREA: 100.0, CSL1085LFP: 5.0},
    )


@pytest.fixture
def ga_blue_ridge_basin() -> BasinCharacteristics:
    return BasinCharacteristics(
        site_no="02333500",
        site_name="Chattahoochee River near Helen, GA",
        region=GA_BLUE_RIDGE,
        predictors={DRNAREA: 47.0, ELEV: 2650.0},
    )


@pytest.fixture
def mt_mountain_basin() -> BasinCharacteristics:
    return BasinCharacteristics(
        site_no="MT-UNGAGED",
        site_name="Lost Horse Creek near Hamilton, MT",
        region=MT_MOUNTAIN,
        predictors={DRNAREA: 45.0, ELEV: 5800.0, PRECIP: 32.0},
    )


# ---------------------------------------------------------------------------
# HydrologicRegion
# ---------------------------------------------------------------------------


class TestHydrologicRegion:
    def test_frozen_and_hashable(self):
        r = HydrologicRegion("GA-1", "Blue Ridge", "GA", required_predictors=("DRNAREA", "ELEV"))
        d = {r: "ok"}
        assert d[r] == "ok"

    def test_equality_by_fields(self):
        r1 = HydrologicRegion("TN-2", "Middle TN", "TN", required_predictors=("DRNAREA",))
        r2 = HydrologicRegion("TN-2", "Middle TN", "TN", required_predictors=("DRNAREA",))
        assert r1 == r2

    def test_validate_predictors_ok(self):
        TN_AREA2.validate_predictors({DRNAREA: 100.0, CSL1085LFP: 5.0})

    def test_validate_predictors_missing_raises(self):
        with pytest.raises(ValueError, match="missing"):
            TN_AREA2.validate_predictors({DRNAREA: 100.0})  # CSL1085LFP missing

    def test_validate_predictors_nonpositive_raises(self):
        with pytest.raises(ValueError, match="non-positive"):
            TN_AREA2.validate_predictors({DRNAREA: 100.0, CSL1085LFP: -1.0})

    def test_validate_three_predictor_region(self):
        MT_MOUNTAIN.validate_predictors({DRNAREA: 45.0, ELEV: 5800.0, PRECIP: 32.0})

    def test_has_predictor(self):
        assert GA_BLUE_RIDGE.has_predictor(ELEV)
        assert not GA_BLUE_RIDGE.has_predictor(PRECIP)

    def test_from_dict_round_trip(self):
        original = TN_AREA3
        restored = HydrologicRegion.from_dict(original.to_dict())
        assert restored.code == original.code
        assert restored.state == original.state
        assert set(restored.required_predictors) == set(original.required_predictors)

    def test_str_contains_code_and_state(self):
        s = str(TN_AREA2)
        assert "TN-2" in s
        assert "TN" in s

    def test_tn_regions_length(self):
        assert len(TN_REGIONS) == 4

    def test_ga_regions_length(self):
        assert len(GA_REGIONS) == 4

    def test_mt_regions_length(self):
        assert len(MT_REGIONS) == 3

    def test_region_codes_unique_across_states(self):
        all_codes = [r.code for r in TN_REGIONS + GA_REGIONS + MT_REGIONS]
        assert len(all_codes) == len(set(all_codes))


# ---------------------------------------------------------------------------
# RegionRegistry
# ---------------------------------------------------------------------------


class TestRegionRegistry:
    def test_register_and_lookup(self):
        reg = RegionRegistry()
        reg.register(TN_AREA2)
        assert reg["TN-2"] is TN_AREA2

    def test_contains(self):
        reg = RegionRegistry()
        reg.register(GA_BLUE_RIDGE)
        assert "GA-1" in reg
        assert "XX-99" not in reg

    def test_missing_key_raises(self):
        reg = RegionRegistry()
        with pytest.raises(KeyError):
            _ = reg["NONEXISTENT"]

    def test_all_regions_returns_list(self):
        reg = RegionRegistry()
        for r in TN_REGIONS:
            reg.register(r)
        assert len(reg.all_regions()) == 4

    def test_load_from_csv(self, tmp_path):
        csv_path = tmp_path / "regions.csv"
        csv_path.write_text(
            "code,label,state,required_predictors\n"
            "GA-1,Blue Ridge,GA,DRNAREA ELEV\n"
            "MT-1,Mountain,MT,DRNAREA ELEV PRECIP\n"
        )
        reg = RegionRegistry.load_from_csv(csv_path)
        assert "GA-1" in reg
        assert "MT-1" in reg


# ---------------------------------------------------------------------------
# BasinCharacteristics
# ---------------------------------------------------------------------------


class TestBasinCharacteristics:
    def test_tn2_valid(self, tn2_basin):
        assert tn2_basin.site_no == "03606500"
        assert tn2_basin.predictors[DRNAREA] == 100.0
        assert tn2_basin.predictors[CSL1085LFP] == 5.0

    def test_ga_blue_ridge_valid(self, ga_blue_ridge_basin):
        assert ga_blue_ridge_basin.predictors[ELEV] == 2650.0
        assert PRECIP not in ga_blue_ridge_basin.predictors

    def test_mt_mountain_three_predictors(self, mt_mountain_basin):
        assert mt_mountain_basin.predictors[ELEV] == 5800.0
        assert mt_mountain_basin.predictors[PRECIP] == 32.0

    def test_ga_coastal_requires_precip(self):
        basin = BasinCharacteristics(
            site_no="X",
            site_name="X",
            region=GA_COASTAL,
            predictors={DRNAREA: 200.0, PRECIP: 50.0},
        )
        assert basin.predictors[PRECIP] == 50.0

    def test_missing_required_raises(self):
        with pytest.raises(ValueError, match="missing"):
            BasinCharacteristics(
                site_no="X",
                site_name="X",
                region=MT_MOUNTAIN,
                predictors={DRNAREA: 45.0, ELEV: 5800.0},  # PRECIP missing
            )

    def test_nonpositive_predictor_raises(self):
        with pytest.raises(ValueError, match="non-positive"):
            BasinCharacteristics(
                site_no="X",
                site_name="X",
                region=TN_AREA2,
                predictors={DRNAREA: 0.0, CSL1085LFP: 5.0},
            )

    def test_state_property_from_region(self, tn2_basin):
        assert tn2_basin.state == "TN"

    def test_predictor_value_ok(self, ga_blue_ridge_basin):
        assert ga_blue_ridge_basin.predictor_value(ELEV) == 2650.0

    def test_predictor_value_missing_raises(self, ga_blue_ridge_basin):
        with pytest.raises(KeyError):
            ga_blue_ridge_basin.predictor_value(PRECIP)

    def test_extra_predictors_stored(self):
        basin = BasinCharacteristics(
            site_no="X",
            site_name="X",
            region=TN_AREA3,
            predictors={DRNAREA: 50.0, FOREST: 65.0},
        )
        assert basin.predictors[FOREST] == 65.0

    def test_from_streamstats_filters_nonnumeric(self):
        params = {DRNAREA: 100.0, CSL1085LFP: 5.0, "STATION_ID": "non-numeric"}
        basin = BasinCharacteristics.from_streamstats("TEST", "Test", TN_AREA2, params)
        assert DRNAREA in basin.predictors
        assert "STATION_ID" not in basin.predictors

    def test_from_streamstats_missing_required_raises(self):
        with pytest.raises(ValueError):
            BasinCharacteristics.from_streamstats(
                "X", "X", GA_BLUE_RIDGE, {DRNAREA: 100.0}  # ELEV missing
            )

    def test_summary_contains_state_and_predictors(self, tn2_basin):
        s = tn2_basin.summary()
        assert "TN" in s
        assert DRNAREA in s

    def test_no_slope_convenience_property(self, tn2_basin):
        assert not hasattr(tn2_basin, "slope_1085_ftmi")

    def test_no_drainage_area_convenience_property(self, tn2_basin):
        assert not hasattr(tn2_basin, "drainage_area_sqmi")


# ---------------------------------------------------------------------------
# GlsEquation
# ---------------------------------------------------------------------------


class TestGlsEquation:
    def test_compute_tn2_q100(self, tn2_basin):
        eq = GlsEquation(
            region=TN_AREA2,
            aep=0.01,
            intercept=2.69,
            coefficients={DRNAREA: 0.75, CSL1085LFP: 0.38},
            sep_pct=35.0,
        )
        q = eq.compute(tn2_basin)
        expected_log = 2.69 + 0.75 * math.log10(100) + 0.38 * math.log10(5)
        assert math.isclose(math.log10(q), expected_log, rel_tol=1e-9)

    def test_compute_ga_blue_ridge(self, ga_blue_ridge_basin):
        eq = GlsEquation(
            region=GA_BLUE_RIDGE,
            aep=0.01,
            intercept=1.14,
            coefficients={DRNAREA: 0.82, ELEV: 0.35},
        )
        q = eq.compute(ga_blue_ridge_basin)
        expected_log = 1.14 + 0.82 * math.log10(47) + 0.35 * math.log10(2650)
        assert math.isclose(math.log10(q), expected_log, rel_tol=1e-9)

    def test_compute_mt_three_predictors(self, mt_mountain_basin):
        eq = GlsEquation(
            region=MT_MOUNTAIN,
            aep=0.01,
            intercept=-0.61,
            coefficients={DRNAREA: 0.76, ELEV: 0.48, PRECIP: 0.62},
        )
        q = eq.compute(mt_mountain_basin)
        expected_log = (
            -0.61 + 0.76 * math.log10(45.0) + 0.48 * math.log10(5800.0) + 0.62 * math.log10(32.0)
        )
        assert math.isclose(math.log10(q), expected_log, rel_tol=1e-9)

    def test_compute_da_only_mt_plains(self):
        basin = BasinCharacteristics(
            site_no="X", site_name="X", region=MT_PLAINS, predictors={DRNAREA: 100.0}
        )
        eq = GlsEquation(region=MT_PLAINS, aep=0.01, intercept=1.48, coefficients={DRNAREA: 0.71})
        q = eq.compute(basin)
        expected = 10 ** (1.48 + 0.71 * math.log10(100.0))
        assert math.isclose(q, expected, rel_tol=1e-9)

    def test_wrong_region_raises(self, tn2_basin):
        eq = GlsEquation(region=TN_AREA3, aep=0.01, intercept=2.62, coefficients={DRNAREA: 0.78})
        with pytest.raises(ValueError, match="region"):
            eq.compute(tn2_basin)

    def test_variance_log10_from_sep(self):
        eq = GlsEquation(
            region=TN_AREA2,
            aep=0.01,
            intercept=2.69,
            coefficients={DRNAREA: 0.75, CSL1085LFP: 0.38},
            sep_pct=35.0,
        )
        assert math.isclose(eq.variance_log10, (math.log10(1.35)) ** 2, rel_tol=1e-9)

    def test_variance_none_without_sep(self):
        eq = GlsEquation(region=MT_PLAINS, aep=0.01, intercept=1.48, coefficients={DRNAREA: 0.71})
        assert eq.variance_log10 is None

    def test_return_period(self):
        eq = GlsEquation(region=TN_AREA3, aep=0.01, intercept=2.62, coefficients={DRNAREA: 0.78})
        assert eq.return_period == pytest.approx(100.0)

    def test_predictor_codes_sorted(self):
        eq = GlsEquation(
            region=MT_MOUNTAIN,
            aep=0.01,
            intercept=-0.61,
            coefficients={DRNAREA: 0.76, ELEV: 0.48, PRECIP: 0.62},
        )
        assert eq.predictor_codes == sorted([DRNAREA, ELEV, PRECIP])

    def test_to_dict_has_region_code_and_aep(self):
        eq = GlsEquation(
            region=GA_PIEDMONT,
            aep=0.01,
            intercept=2.65,
            coefficients={DRNAREA: 0.77},
            sep_pct=44.0,
        )
        d = eq.to_dict()
        assert d["region_code"] == "GA-3"
        assert d["aep"] == 0.01
        assert d[DRNAREA] == 0.77


# ---------------------------------------------------------------------------
# Tennessee table
# ---------------------------------------------------------------------------


class TestTennesseeTable:
    def test_equation_count(self, tn_table):
        assert len(tn_table) == 4 * 8

    def test_available_states(self, tn_table):
        assert tn_table.available_states() == ["TN"]

    def test_available_region_codes(self, tn_table):
        codes = {r.code for r in tn_table.available_regions()}
        assert codes == {"TN-1", "TN-2", "TN-3", "TN-4"}

    def test_area2_estimate_positive(self, tn_table, tn2_basin):
        assert tn_table.estimate(tn2_basin, aep=0.01) > 0

    def test_area2_monotonic_with_return_period(self, tn_table, tn2_basin):
        qs = tn_table.estimate_all_aeps(tn2_basin)
        sorted_aeps = sorted(qs.keys(), reverse=True)
        flows = [qs[a] for a in sorted_aeps]
        assert all(flows[i] <= flows[i + 1] for i in range(len(flows) - 1))

    def test_area3_da_only(self, tn_table):
        basin = BasinCharacteristics(
            site_no="X", site_name="X", region=TN_AREA3, predictors={DRNAREA: 50.0}
        )
        assert tn_table.estimate(basin, aep=0.01) > 0

    def test_area1_i2_predictor(self, tn_table):
        basin = BasinCharacteristics(
            site_no="X",
            site_name="X",
            region=TN_AREA1,
            predictors={DRNAREA: 80.0, I2: 2.5},
        )
        assert tn_table.estimate(basin, aep=0.01) > 0

    def test_area4_imperv_predictor(self, tn_table):
        basin = BasinCharacteristics(
            site_no="X",
            site_name="X",
            region=TN_AREA4,
            predictors={DRNAREA: 60.0, IMPERV: 5.0},
        )
        assert tn_table.estimate(basin, aep=0.01) > 0

    def test_sep_pct_populated(self, tn_table):
        eq = tn_table.get_equation(TN_AREA2, 0.01)
        assert eq.sep_pct is not None and eq.sep_pct > 0

    def test_summary_table_length(self, tn_table, tn2_basin):
        rows = tn_table.summary_table(tn2_basin)
        assert len(rows) == len(STANDARD_AEPS)
        assert all(r["flow_cfs"] > 0 for r in rows)

    def test_publication_set(self, tn_table):
        assert "2024" in tn_table.publication or "SIR" in tn_table.publication


# ---------------------------------------------------------------------------
# Georgia table
# ---------------------------------------------------------------------------


class TestGeorgiaTable:
    def test_equation_count(self, ga_table):
        assert len(ga_table) == 4 * 8

    def test_available_states(self, ga_table):
        assert ga_table.available_states() == ["GA"]

    def test_blue_ridge_has_elev_not_precip(self, ga_table):
        eq = ga_table.get_equation(GA_BLUE_RIDGE, 0.01)
        assert ELEV in eq.coefficients
        assert PRECIP not in eq.coefficients

    def test_coastal_has_precip_not_elev(self, ga_table):
        eq = ga_table.get_equation(GA_COASTAL, 0.01)
        assert PRECIP in eq.coefficients
        assert ELEV not in eq.coefficients

    def test_valley_ridge_has_slope(self, ga_table):
        eq = ga_table.get_equation(GA_VALLEY_RIDGE, 0.01)
        assert CSL1085LFP in eq.coefficients

    def test_piedmont_da_only(self, ga_table):
        eq = ga_table.get_equation(GA_PIEDMONT, 0.01)
        assert list(eq.coefficients.keys()) == [DRNAREA]

    def test_blue_ridge_estimate(self, ga_table, ga_blue_ridge_basin):
        assert ga_table.estimate(ga_blue_ridge_basin, aep=0.01) > 0

    def test_all_ga_regions_monotonic(self, ga_table):
        basins = [
            BasinCharacteristics("A", "A", GA_BLUE_RIDGE, {DRNAREA: 50.0, ELEV: 2500.0}),
            BasinCharacteristics("B", "B", GA_VALLEY_RIDGE, {DRNAREA: 80.0, CSL1085LFP: 6.0}),
            BasinCharacteristics("C", "C", GA_PIEDMONT, {DRNAREA: 100.0}),
            BasinCharacteristics("D", "D", GA_COASTAL, {DRNAREA: 200.0, PRECIP: 52.0}),
        ]
        for basin in basins:
            qs = ga_table.estimate_all_aeps(basin)
            sorted_aeps = sorted(qs.keys(), reverse=True)
            flows = [qs[a] for a in sorted_aeps]
            assert all(
                flows[i] <= flows[i + 1] for i in range(len(flows) - 1)
            ), f"Non-monotonic flows for {basin.site_no}"


# ---------------------------------------------------------------------------
# Montana table
# ---------------------------------------------------------------------------


class TestMontanaTable:
    def test_equation_count(self, mt_table):
        assert len(mt_table) == 3 * 8

    def test_available_states(self, mt_table):
        assert mt_table.available_states() == ["MT"]

    def test_mountain_three_predictors(self, mt_table):
        eq = mt_table.get_equation(MT_MOUNTAIN, 0.01)
        assert set(eq.coefficients.keys()) == {DRNAREA, ELEV, PRECIP}

    def test_foothills_two_predictors(self, mt_table):
        eq = mt_table.get_equation(MT_FOOTHILLS, 0.01)
        assert set(eq.coefficients.keys()) == {DRNAREA, CSL1085LFP}

    def test_plains_one_predictor(self, mt_table):
        eq = mt_table.get_equation(MT_PLAINS, 0.01)
        assert list(eq.coefficients.keys()) == [DRNAREA]

    def test_mountain_estimate_positive(self, mt_table, mt_mountain_basin):
        assert mt_table.estimate(mt_mountain_basin, aep=0.01) > 0

    def test_mountain_q100_greater_than_q2(self, mt_table, mt_mountain_basin):
        assert mt_table.estimate(mt_mountain_basin, aep=0.01) > mt_table.estimate(
            mt_mountain_basin, aep=0.5
        )

    def test_plains_estimate_positive(self, mt_table):
        basin = BasinCharacteristics(
            site_no="X", site_name="X", region=MT_PLAINS, predictors={DRNAREA: 500.0}
        )
        assert mt_table.estimate(basin, aep=0.01) > 0

    def test_plains_higher_sep_than_mountain(self, mt_table):
        eq_mtn = mt_table.get_equation(MT_MOUNTAIN, 0.01)
        eq_plains = mt_table.get_equation(MT_PLAINS, 0.01)
        assert eq_plains.sep_pct > eq_mtn.sep_pct


# ---------------------------------------------------------------------------
# Nationwide / multi-state
# ---------------------------------------------------------------------------


class TestNationwideTable:
    def test_states(self, national_table):
        assert set(national_table.available_states()) == {"TN", "GA", "MT"}

    def test_equation_count(self, national_table):
        # (4 TN + 4 GA + 3 MT) × 8 AEPs = 88
        assert len(national_table) == 88

    def test_regions_for_state_tn(self, national_table):
        regs = national_table.regions_for_state("TN")
        assert len(regs) == 4
        assert all(r.state == "TN" for r in regs)

    def test_regions_for_state_mt(self, national_table):
        regs = national_table.regions_for_state("MT")
        assert len(regs) == 3

    def test_filter_by_state_tn(self, national_table):
        tn_only = national_table.filter_by_state("TN")
        assert tn_only.available_states() == ["TN"]
        assert len(tn_only) == 32

    def test_filter_by_state_ga(self, national_table):
        ga_only = national_table.filter_by_state("GA")
        assert ga_only.available_states() == ["GA"]
        assert len(ga_only) == 32

    def test_batch_estimate_three_states(
        self, national_table, tn2_basin, ga_blue_ridge_basin, mt_mountain_basin
    ):
        results = national_table.batch_estimate(
            [tn2_basin, ga_blue_ridge_basin, mt_mountain_basin], aep=0.01
        )
        assert set(results.keys()) == {"03606500", "02333500", "MT-UNGAGED"}
        assert all(q > 0 for q in results.values())

    def test_batch_estimate_skips_unknown_region(self, national_table, tn2_basin):
        unknown_region = HydrologicRegion(
            "XX-99", "Unknown", "XX", required_predictors=("DRNAREA",)
        )
        unknown_basin = BasinCharacteristics(
            site_no="UNK",
            site_name="Unknown",
            region=unknown_region,
            predictors={DRNAREA: 100.0},
        )
        results = national_table.batch_estimate([tn2_basin, unknown_basin], aep=0.01)
        assert "03606500" in results
        assert "UNK" not in results

    def test_get_by_region_code(self, national_table):
        eq_tn = national_table.get_by_region_code("TN-2", 0.01)
        assert eq_tn.region.state == "TN"
        eq_ga = national_table.get_by_region_code("GA-1", 0.01)
        assert eq_ga.region.state == "GA"
        eq_mt = national_table.get_by_region_code("MT-1", 0.01)
        assert eq_mt.region.state == "MT"

    def test_repr(self, national_table):
        r = repr(national_table)
        assert "88" in r
        assert "GA" in r or "TN" in r


# ---------------------------------------------------------------------------
# CSV round-trip
# ---------------------------------------------------------------------------


class TestCSVRoundTrip:
    def test_round_trip_tn(self, tn_table, tmp_path):
        p = tmp_path / "tn.csv"
        tn_table.to_csv(p)
        loaded = RegressionTable.load_from_csv(p)
        assert len(loaded) == len(tn_table)

    def test_multi_state_predictor_isolation(self, tmp_path):
        """In a merged GA+MT CSV, each region row only has its own predictor."""
        merged = RegressionTable.merge(build_georgia_table(), build_montana_table())
        p = tmp_path / "ga_mt.csv"
        merged.to_csv(p)
        loaded = RegressionTable.load_from_csv(p)
        # GA Blue Ridge: ELEV present, CSL1085LFP absent
        eq_ga = loaded.get_by_region_code("GA-1", 0.01)
        assert ELEV in eq_ga.coefficients
        assert CSL1085LFP not in eq_ga.coefficients
        # MT Foothills: CSL1085LFP present, ELEV absent (ELEV is in Mountain rows)
        eq_fh = loaded.get_by_region_code("MT-2", 0.01)
        assert CSL1085LFP in eq_fh.coefficients

    def test_mt_mountain_three_predictors_survive_csv(self, mt_table, tmp_path):
        p = tmp_path / "mt.csv"
        mt_table.to_csv(p)
        loaded = RegressionTable.load_from_csv(p)
        eq = loaded.get_by_region_code("MT-1", 0.01)
        assert set(eq.coefficients.keys()) == {DRNAREA, ELEV, PRECIP}

    def test_file_not_found_raises(self):
        with pytest.raises(FileNotFoundError):
            RegressionTable.load_from_csv("/nonexistent/path/equations.csv")

    def test_write_template_all_three_states(self, tmp_path):
        p = tmp_path / "template.csv"
        all_regions = list(TN_REGIONS) + list(GA_REGIONS) + list(MT_REGIONS)
        RegressionTable.write_template(p, all_regions)
        text = p.read_text()
        assert "TN-1" in text
        assert "GA-1" in text
        assert "MT-1" in text
        data_lines = [l for l in text.splitlines()[1:] if l]
        assert len(data_lines) == len(all_regions) * len(STANDARD_AEPS)

    def test_load_from_dict_new_style(self):
        d = {
            (TN_AREA2, 0.01): {
                "intercept": 2.69,
                "coefficients": {DRNAREA: 0.75, CSL1085LFP: 0.38},
                "sep_pct": 35.0,
            },
            (GA_BLUE_RIDGE, 0.01): {
                "intercept": 1.14,
                "coefficients": {DRNAREA: 0.82, ELEV: 0.35},
                "sep_pct": 38.0,
            },
            (MT_MOUNTAIN, 0.01): {
                "intercept": -0.61,
                "coefficients": {DRNAREA: 0.76, ELEV: 0.48, PRECIP: 0.62},
            },
        }
        table = RegressionTable.load_from_dict(d)
        assert len(table) == 3
        assert set(table.available_states()) == {"TN", "GA", "MT"}

    def test_load_from_dict_missing_intercept_raises(self):
        with pytest.raises(KeyError):
            RegressionTable.load_from_dict({(TN_AREA2, 0.01): {"coefficients": {DRNAREA: 0.75}}})
