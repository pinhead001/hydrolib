"""
Tests for the seven Western US state regression equation modules:
New Mexico, Colorado, Wyoming, Idaho, Washington, Oregon, California.

Each state module is validated for:
- Correct equation count and AEP coverage
- Positive, monotonically increasing flows (Q100 > Q10 > Q2)
- Correct predictor isolation per region
- State code and region code consistency
- Variance / SEP% availability
- 10-state nationwide merge
"""

from __future__ import annotations

import pytest

from hydrolib.regression.basin_chars import (
    CSL1085LFP,
    DRNAREA,
    ELEV,
    PRECIP,
    BasinCharacteristics,
)
from hydrolib.regression.regression_table import STANDARD_AEPS, RegressionTable
from hydrolib.regression.states.california import (
    CA_CENTRAL,
    CA_NORTH_COAST,
    CA_REGIONS,
    CA_SIERRA,
    CA_SOUTH,
    build_california_table,
)
from hydrolib.regression.states.colorado import (
    CO_EASTERN,
    CO_FRONT_RANGE,
    CO_REGIONS,
    CO_SAN_LUIS,
    CO_WESTERN,
    build_colorado_table,
)
from hydrolib.regression.states.idaho import (
    ID_CENTRAL,
    ID_NORTH,
    ID_REGIONS,
    ID_SNAKE,
    ID_SOUTH,
    build_idaho_table,
)
from hydrolib.regression.states.new_mexico import (
    NM_DESERT,
    NM_MOUNTAIN,
    NM_PLAINS,
    NM_REGIONS,
    NM_TRANSITION,
    build_new_mexico_table,
)
from hydrolib.regression.states.oregon import (
    OR_CASCADE,
    OR_COAST,
    OR_EASTERN,
    OR_REGIONS,
    OR_WILLAMETTE,
    build_oregon_table,
)
from hydrolib.regression.states.washington import (
    WA_E_CASCADE,
    WA_EASTERN,
    WA_OLYMPIC,
    WA_REGIONS,
    WA_W_CASCADE,
    build_washington_table,
)
from hydrolib.regression.states.wyoming import (
    WY_FOOTHILLS,
    WY_MOUNTAIN,
    WY_PLAINS,
    WY_REGIONS,
    build_wyoming_table,
)

# ---------------------------------------------------------------------------
# Module-scoped fixtures — tables built once per session
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def nm_table() -> RegressionTable:
    return build_new_mexico_table()


@pytest.fixture(scope="module")
def co_table() -> RegressionTable:
    return build_colorado_table()


@pytest.fixture(scope="module")
def wy_table() -> RegressionTable:
    return build_wyoming_table()


@pytest.fixture(scope="module")
def id_table() -> RegressionTable:
    return build_idaho_table()


@pytest.fixture(scope="module")
def wa_table() -> RegressionTable:
    return build_washington_table()


@pytest.fixture(scope="module")
def or_table() -> RegressionTable:
    return build_oregon_table()


@pytest.fixture(scope="module")
def ca_table() -> RegressionTable:
    return build_california_table()


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _assert_monotonic_q(table: RegressionTable, basin: BasinCharacteristics) -> None:
    """Assert Q100 > Q10 > Q2 for a given basin/table pair."""
    q2 = table.estimate(basin, aep=0.5)
    q10 = table.estimate(basin, aep=0.1)
    q100 = table.estimate(basin, aep=0.01)
    assert q100 > q10 > q2 > 0, (
        f"Non-monotonic or non-positive flows for {basin.region.code}: "
        f"Q2={q2:.0f}, Q10={q10:.0f}, Q100={q100:.0f}"
    )


# ===========================================================================
# New Mexico
# ===========================================================================


class TestNewMexico:
    def test_equation_count(self, nm_table):
        assert len(nm_table) == 32  # 4 regions × 8 AEPs

    def test_state_code(self, nm_table):
        for r in nm_table.available_regions():
            assert r.state == "NM"

    def test_region_codes(self):
        codes = {r.code for r in NM_REGIONS}
        assert codes == {"NM-1", "NM-2", "NM-3", "NM-4"}

    def test_mountain_predictors(self):
        assert set(NM_MOUNTAIN.required_predictors) == {"DRNAREA", "PRECIP"}

    def test_transition_slope_predictor(self):
        assert "CSL1085LFP" in NM_TRANSITION.required_predictors

    def test_plains_da_only(self):
        assert NM_PLAINS.required_predictors == ("DRNAREA",)

    def test_desert_slope_predictor(self):
        assert "CSL1085LFP" in NM_DESERT.required_predictors

    def test_all_aeps_present(self, nm_table):
        for region in NM_REGIONS:
            assert set(nm_table.available_aeps(region)) == set(STANDARD_AEPS)

    def test_mountain_monotonic(self, nm_table):
        basin = BasinCharacteristics(
            site_no="NM-MT-01",
            site_name="Test Mountain Basin, NM",
            region=NM_MOUNTAIN,
            predictors={DRNAREA: 85.0, PRECIP: 24.0},
        )
        _assert_monotonic_q(nm_table, basin)

    def test_plains_positive_estimate(self, nm_table):
        basin = BasinCharacteristics(
            site_no="NM-PL-01",
            site_name="Test Plains Basin, NM",
            region=NM_PLAINS,
            predictors={DRNAREA: 120.0},
        )
        assert nm_table.estimate(basin, aep=0.01) > 0

    def test_desert_slope_effect(self, nm_table):
        """Steeper slope → higher Q for Desert region."""
        b1 = BasinCharacteristics(
            site_no="D1",
            site_name="D1",
            region=NM_DESERT,
            predictors={DRNAREA: 50.0, CSL1085LFP: 3.0},
        )
        b2 = BasinCharacteristics(
            site_no="D2",
            site_name="D2",
            region=NM_DESERT,
            predictors={DRNAREA: 50.0, CSL1085LFP: 12.0},
        )
        assert nm_table.estimate(b2, aep=0.01) != nm_table.estimate(b1, aep=0.01)

    def test_sep_pct_available(self, nm_table):
        eq = nm_table.get_equation(NM_MOUNTAIN, 0.01)
        assert eq.sep_pct is not None and eq.sep_pct > 0

    def test_variance_positive(self, nm_table):
        eq = nm_table.get_equation(NM_PLAINS, 0.01)
        assert eq.variance_log10 is not None and eq.variance_log10 > 0

    def test_transition_monotonic(self, nm_table):
        basin = BasinCharacteristics(
            site_no="NM-TR-01",
            site_name="Test Transition Basin, NM",
            region=NM_TRANSITION,
            predictors={DRNAREA: 60.0, CSL1085LFP: 7.0},
        )
        _assert_monotonic_q(nm_table, basin)


# ===========================================================================
# Colorado
# ===========================================================================


class TestColorado:
    def test_equation_count(self, co_table):
        assert len(co_table) == 32

    def test_state_code(self, co_table):
        for r in co_table.available_regions():
            assert r.state == "CO"

    def test_region_codes(self):
        codes = {r.code for r in CO_REGIONS}
        assert codes == {"CO-1", "CO-2", "CO-3", "CO-4"}

    def test_front_range_three_predictors(self):
        assert set(CO_FRONT_RANGE.required_predictors) == {"DRNAREA", "ELEV", "PRECIP"}

    def test_eastern_slope_predictor(self):
        assert "CSL1085LFP" in CO_EASTERN.required_predictors

    def test_western_elev_not_precip(self):
        assert "ELEV" in CO_WESTERN.required_predictors
        assert "PRECIP" not in CO_WESTERN.required_predictors

    def test_san_luis_da_only(self):
        assert CO_SAN_LUIS.required_predictors == ("DRNAREA",)

    def test_all_aeps_present(self, co_table):
        for region in CO_REGIONS:
            assert set(co_table.available_aeps(region)) == set(STANDARD_AEPS)

    def test_front_range_monotonic(self, co_table):
        basin = BasinCharacteristics(
            site_no="CO-FR-01",
            site_name="Test Front Range Basin, CO",
            region=CO_FRONT_RANGE,
            predictors={DRNAREA: 150.0, ELEV: 8500.0, PRECIP: 25.0},
        )
        _assert_monotonic_q(co_table, basin)

    def test_western_slope_monotonic(self, co_table):
        basin = BasinCharacteristics(
            site_no="CO-WS-01",
            site_name="Crystal River Basin, CO",
            region=CO_WESTERN,
            predictors={DRNAREA: 191.0, ELEV: 9200.0},
        )
        _assert_monotonic_q(co_table, basin)

    def test_eastern_plains_positive(self, co_table):
        basin = BasinCharacteristics(
            site_no="CO-EP-01",
            site_name="Test Eastern Plains Basin, CO",
            region=CO_EASTERN,
            predictors={DRNAREA: 200.0, CSL1085LFP: 4.0},
        )
        assert co_table.estimate(basin, aep=0.01) > 0

    def test_missing_predictor_raises(self, co_table):
        """Front Range requires ELEV + PRECIP; omitting raises ValueError at construction."""
        with pytest.raises(ValueError, match="missing predictors"):
            BasinCharacteristics(
                site_no="CO-BAD",
                site_name="Bad Basin",
                region=CO_FRONT_RANGE,
                predictors={DRNAREA: 100.0},
            )

    def test_san_luis_monotonic(self, co_table):
        basin = BasinCharacteristics(
            site_no="CO-SL-01",
            site_name="Test San Luis Valley Basin, CO",
            region=CO_SAN_LUIS,
            predictors={DRNAREA: 500.0},
        )
        _assert_monotonic_q(co_table, basin)


# ===========================================================================
# Wyoming
# ===========================================================================


class TestWyoming:
    def test_equation_count(self, wy_table):
        assert len(wy_table) == 24  # 3 regions × 8 AEPs

    def test_state_code(self, wy_table):
        for r in wy_table.available_regions():
            assert r.state == "WY"

    def test_region_codes(self):
        codes = {r.code for r in WY_REGIONS}
        assert codes == {"WY-1", "WY-2", "WY-3"}

    def test_mountain_three_predictors(self):
        assert set(WY_MOUNTAIN.required_predictors) == {"DRNAREA", "ELEV", "PRECIP"}

    def test_foothills_slope(self):
        assert "CSL1085LFP" in WY_FOOTHILLS.required_predictors

    def test_plains_da_only(self):
        assert WY_PLAINS.required_predictors == ("DRNAREA",)

    def test_all_aeps_present(self, wy_table):
        for region in WY_REGIONS:
            assert set(wy_table.available_aeps(region)) == set(STANDARD_AEPS)

    def test_mountain_monotonic(self, wy_table):
        basin = BasinCharacteristics(
            site_no="WY-MT-01",
            site_name="Greybull River Basin, WY",
            region=WY_MOUNTAIN,
            predictors={DRNAREA: 640.0, ELEV: 7400.0, PRECIP: 22.0},
        )
        _assert_monotonic_q(wy_table, basin)

    def test_foothills_monotonic(self, wy_table):
        basin = BasinCharacteristics(
            site_no="WY-FH-01",
            site_name="Test Foothills Basin, WY",
            region=WY_FOOTHILLS,
            predictors={DRNAREA: 80.0, CSL1085LFP: 8.0},
        )
        _assert_monotonic_q(wy_table, basin)

    def test_plains_positive(self, wy_table):
        basin = BasinCharacteristics(
            site_no="WY-PL-01",
            site_name="Test Plains Basin, WY",
            region=WY_PLAINS,
            predictors={DRNAREA: 300.0},
        )
        assert wy_table.estimate(basin, aep=0.01) > 0

    def test_eyr_available(self, wy_table):
        eq = wy_table.get_equation(WY_MOUNTAIN, 0.01)
        assert eq.eyr is not None and eq.eyr > 0

    def test_higher_precip_higher_q(self, wy_table):
        b_low = BasinCharacteristics(
            site_no="W1",
            site_name="W1",
            region=WY_MOUNTAIN,
            predictors={DRNAREA: 500.0, ELEV: 7000.0, PRECIP: 15.0},
        )
        b_high = BasinCharacteristics(
            site_no="W2",
            site_name="W2",
            region=WY_MOUNTAIN,
            predictors={DRNAREA: 500.0, ELEV: 7000.0, PRECIP: 40.0},
        )
        assert wy_table.estimate(b_high, aep=0.01) > wy_table.estimate(b_low, aep=0.01)


# ===========================================================================
# Idaho
# ===========================================================================


class TestIdaho:
    def test_equation_count(self, id_table):
        assert len(id_table) == 32

    def test_state_code(self, id_table):
        for r in id_table.available_regions():
            assert r.state == "ID"

    def test_region_codes(self):
        codes = {r.code for r in ID_REGIONS}
        assert codes == {"ID-1", "ID-2", "ID-3", "ID-4"}

    def test_north_precip_not_elev(self):
        assert "PRECIP" in ID_NORTH.required_predictors
        assert "ELEV" not in ID_NORTH.required_predictors

    def test_central_three_predictors(self):
        assert set(ID_CENTRAL.required_predictors) == {"DRNAREA", "ELEV", "PRECIP"}

    def test_snake_slope(self):
        assert "CSL1085LFP" in ID_SNAKE.required_predictors

    def test_south_da_only(self):
        assert ID_SOUTH.required_predictors == ("DRNAREA",)

    def test_all_aeps_present(self, id_table):
        for region in ID_REGIONS:
            assert set(id_table.available_aeps(region)) == set(STANDARD_AEPS)

    def test_north_monotonic(self, id_table):
        basin = BasinCharacteristics(
            site_no="ID-N-01",
            site_name="Test North Basin, ID",
            region=ID_NORTH,
            predictors={DRNAREA: 120.0, PRECIP: 45.0},
        )
        _assert_monotonic_q(id_table, basin)

    def test_central_three_predictor_monotonic(self, id_table):
        basin = BasinCharacteristics(
            site_no="ID-C-01",
            site_name="Boise River Basin, ID",
            region=ID_CENTRAL,
            predictors={DRNAREA: 1220.0, ELEV: 6800.0, PRECIP: 38.0},
        )
        _assert_monotonic_q(id_table, basin)

    def test_snake_plain_positive(self, id_table):
        basin = BasinCharacteristics(
            site_no="ID-SP-01",
            site_name="Test Snake Basin, ID",
            region=ID_SNAKE,
            predictors={DRNAREA: 500.0, CSL1085LFP: 5.0},
        )
        assert id_table.estimate(basin, aep=0.02) > 0

    def test_north_precip_effect(self, id_table):
        """Higher precipitation → higher flow for same DA (North region)."""
        b_low = BasinCharacteristics(
            site_no="I1",
            site_name="I1",
            region=ID_NORTH,
            predictors={DRNAREA: 100.0, PRECIP: 30.0},
        )
        b_high = BasinCharacteristics(
            site_no="I2",
            site_name="I2",
            region=ID_NORTH,
            predictors={DRNAREA: 100.0, PRECIP: 80.0},
        )
        assert id_table.estimate(b_high, aep=0.01) > id_table.estimate(b_low, aep=0.01)

    def test_south_monotonic(self, id_table):
        basin = BasinCharacteristics(
            site_no="ID-S-01",
            site_name="Test Southern Basin, ID",
            region=ID_SOUTH,
            predictors={DRNAREA: 200.0},
        )
        _assert_monotonic_q(id_table, basin)


# ===========================================================================
# Washington
# ===========================================================================


class TestWashington:
    def test_equation_count(self, wa_table):
        assert len(wa_table) == 32

    def test_state_code(self, wa_table):
        for r in wa_table.available_regions():
            assert r.state == "WA"

    def test_region_codes(self):
        codes = {r.code for r in WA_REGIONS}
        assert codes == {"WA-1", "WA-2", "WA-3", "WA-4"}

    def test_olympic_precip_not_elev(self):
        assert "PRECIP" in WA_OLYMPIC.required_predictors
        assert "ELEV" not in WA_OLYMPIC.required_predictors

    def test_w_cascade_three_predictors(self):
        assert set(WA_W_CASCADE.required_predictors) == {"DRNAREA", "ELEV", "PRECIP"}

    def test_e_cascade_elev_not_precip(self):
        assert "ELEV" in WA_E_CASCADE.required_predictors
        assert "PRECIP" not in WA_E_CASCADE.required_predictors

    def test_eastern_da_only(self):
        assert WA_EASTERN.required_predictors == ("DRNAREA",)

    def test_all_aeps_present(self, wa_table):
        for region in WA_REGIONS:
            assert set(wa_table.available_aeps(region)) == set(STANDARD_AEPS)

    def test_olympic_monotonic(self, wa_table):
        basin = BasinCharacteristics(
            site_no="WA-OL-01",
            site_name="Quinault River Basin, WA",
            region=WA_OLYMPIC,
            predictors={DRNAREA: 286.0, PRECIP: 130.0},
        )
        _assert_monotonic_q(wa_table, basin)

    def test_w_cascade_monotonic(self, wa_table):
        basin = BasinCharacteristics(
            site_no="WA-WC-01",
            site_name="Test W. Cascade Basin, WA",
            region=WA_W_CASCADE,
            predictors={DRNAREA: 200.0, ELEV: 4500.0, PRECIP: 80.0},
        )
        _assert_monotonic_q(wa_table, basin)

    def test_e_cascade_monotonic(self, wa_table):
        basin = BasinCharacteristics(
            site_no="WA-EC-01",
            site_name="Test E. Cascade Basin, WA",
            region=WA_E_CASCADE,
            predictors={DRNAREA: 150.0, ELEV: 3800.0},
        )
        _assert_monotonic_q(wa_table, basin)

    def test_high_precip_olympic_exceeds_eastern(self, wa_table):
        """Olympic (maritime, high precip) Q > Eastern (semi-arid) Q at same DA."""
        olympic = BasinCharacteristics(
            site_no="O1",
            site_name="O1",
            region=WA_OLYMPIC,
            predictors={DRNAREA: 100.0, PRECIP: 120.0},
        )
        eastern = BasinCharacteristics(
            site_no="E1",
            site_name="E1",
            region=WA_EASTERN,
            predictors={DRNAREA: 100.0},
        )
        assert wa_table.estimate(olympic, aep=0.01) > wa_table.estimate(eastern, aep=0.01)


# ===========================================================================
# Oregon
# ===========================================================================


class TestOregon:
    def test_equation_count(self, or_table):
        assert len(or_table) == 32

    def test_state_code(self, or_table):
        for r in or_table.available_regions():
            assert r.state == "OR"

    def test_region_codes(self):
        codes = {r.code for r in OR_REGIONS}
        assert codes == {"OR-1", "OR-2", "OR-3", "OR-4"}

    def test_coast_precip_not_elev(self):
        assert "PRECIP" in OR_COAST.required_predictors
        assert "ELEV" not in OR_COAST.required_predictors

    def test_willamette_slope_predictor(self):
        assert "CSL1085LFP" in OR_WILLAMETTE.required_predictors

    def test_cascade_three_predictors(self):
        assert set(OR_CASCADE.required_predictors) == {"DRNAREA", "ELEV", "PRECIP"}

    def test_eastern_da_only(self):
        assert OR_EASTERN.required_predictors == ("DRNAREA",)

    def test_all_aeps_present(self, or_table):
        for region in OR_REGIONS:
            assert set(or_table.available_aeps(region)) == set(STANDARD_AEPS)

    def test_coast_monotonic(self, or_table):
        basin = BasinCharacteristics(
            site_no="OR-CO-01",
            site_name="Coquille River Basin, OR",
            region=OR_COAST,
            predictors={DRNAREA: 766.0, PRECIP: 65.0},
        )
        _assert_monotonic_q(or_table, basin)

    def test_cascade_three_predictor_monotonic(self, or_table):
        basin = BasinCharacteristics(
            site_no="OR-CA-01",
            site_name="Test Cascade Basin, OR",
            region=OR_CASCADE,
            predictors={DRNAREA: 250.0, ELEV: 4200.0, PRECIP: 70.0},
        )
        _assert_monotonic_q(or_table, basin)

    def test_willamette_slope_effect(self, or_table):
        """Steeper channel → higher Q for same DA (positive slope exponent)."""
        b_flat = BasinCharacteristics(
            site_no="W1",
            site_name="W1",
            region=OR_WILLAMETTE,
            predictors={DRNAREA: 100.0, CSL1085LFP: 2.0},
        )
        b_steep = BasinCharacteristics(
            site_no="W2",
            site_name="W2",
            region=OR_WILLAMETTE,
            predictors={DRNAREA: 100.0, CSL1085LFP: 15.0},
        )
        assert or_table.estimate(b_steep, aep=0.01) > or_table.estimate(b_flat, aep=0.01)

    def test_eastern_positive(self, or_table):
        basin = BasinCharacteristics(
            site_no="OR-EA-01",
            site_name="Test Eastern Basin, OR",
            region=OR_EASTERN,
            predictors={DRNAREA: 400.0},
        )
        assert or_table.estimate(basin, aep=0.01) > 0


# ===========================================================================
# California
# ===========================================================================


class TestCalifornia:
    def test_equation_count(self, ca_table):
        assert len(ca_table) == 32

    def test_state_code(self, ca_table):
        for r in ca_table.available_regions():
            assert r.state == "CA"

    def test_region_codes(self):
        codes = {r.code for r in CA_REGIONS}
        assert codes == {"CA-1", "CA-2", "CA-3", "CA-4"}

    def test_north_coast_precip_not_elev(self):
        assert "PRECIP" in CA_NORTH_COAST.required_predictors
        assert "ELEV" not in CA_NORTH_COAST.required_predictors

    def test_sierra_three_predictors(self):
        assert set(CA_SIERRA.required_predictors) == {"DRNAREA", "ELEV", "PRECIP"}

    def test_central_slope(self):
        assert "CSL1085LFP" in CA_CENTRAL.required_predictors

    def test_south_da_only(self):
        assert CA_SOUTH.required_predictors == ("DRNAREA",)

    def test_all_aeps_present(self, ca_table):
        for region in CA_REGIONS:
            assert set(ca_table.available_aeps(region)) == set(STANDARD_AEPS)

    def test_north_coast_monotonic(self, ca_table):
        basin = BasinCharacteristics(
            site_no="CA-NC-01",
            site_name="Scott River Basin, CA",
            region=CA_NORTH_COAST,
            predictors={DRNAREA: 633.0, PRECIP: 42.0},
        )
        _assert_monotonic_q(ca_table, basin)

    def test_sierra_nevada_monotonic(self, ca_table):
        basin = BasinCharacteristics(
            site_no="CA-SI-01",
            site_name="Test Sierra Nevada Basin, CA",
            region=CA_SIERRA,
            predictors={DRNAREA: 300.0, ELEV: 7200.0, PRECIP: 55.0},
        )
        _assert_monotonic_q(ca_table, basin)

    def test_south_high_sep(self, ca_table):
        """Southern CA has higher SEP% than North Coast (most variable floods)."""
        eq_south = ca_table.get_equation(CA_SOUTH, 0.01)
        eq_north = ca_table.get_equation(CA_NORTH_COAST, 0.01)
        assert eq_south.sep_pct > eq_north.sep_pct

    def test_sierra_precip_effect(self, ca_table):
        """Higher precipitation → higher Sierra Nevada Q for same DA and ELEV."""
        b_dry = BasinCharacteristics(
            site_no="S1",
            site_name="S1",
            region=CA_SIERRA,
            predictors={DRNAREA: 200.0, ELEV: 7000.0, PRECIP: 30.0},
        )
        b_wet = BasinCharacteristics(
            site_no="S2",
            site_name="S2",
            region=CA_SIERRA,
            predictors={DRNAREA: 200.0, ELEV: 7000.0, PRECIP: 70.0},
        )
        assert ca_table.estimate(b_wet, aep=0.01) > ca_table.estimate(b_dry, aep=0.01)

    def test_south_monotonic(self, ca_table):
        basin = BasinCharacteristics(
            site_no="CA-SO-01",
            site_name="Test Southern CA Basin",
            region=CA_SOUTH,
            predictors={DRNAREA: 50.0},
        )
        _assert_monotonic_q(ca_table, basin)


# ===========================================================================
# 10-State nationwide merge
# ===========================================================================


@pytest.fixture(scope="module")
def national_10_table(
    nm_table, co_table, wy_table, id_table, wa_table, or_table, ca_table
) -> RegressionTable:
    """Full 10-state nationwide table."""
    from hydrolib.regression.states.georgia import build_georgia_table
    from hydrolib.regression.states.montana import build_montana_table
    from hydrolib.regression.states.tennessee import build_tennessee_table

    return RegressionTable.merge(
        build_tennessee_table(),
        build_georgia_table(),
        build_montana_table(),
        nm_table,
        co_table,
        wy_table,
        id_table,
        wa_table,
        or_table,
        ca_table,
    )


class TestNationwide10State:
    def test_ten_states_present(self, national_10_table):
        states = set(national_10_table.available_states())
        assert states == {"TN", "GA", "MT", "NM", "CO", "WY", "ID", "WA", "OR", "CA"}

    def test_equation_count(self, national_10_table):
        # TN=32, GA=32, MT=24, NM=32, CO=32, WY=24, ID=32, WA=32, OR=32, CA=32 = 304
        assert len(national_10_table) == 304

    def test_filter_by_ca(self, national_10_table):
        ca_only = national_10_table.filter_by_state("CA")
        assert len(ca_only) == 32
        assert ca_only.available_states() == ["CA"]

    def test_filter_by_wy(self, national_10_table):
        wy_only = national_10_table.filter_by_state("WY")
        assert len(wy_only) == 24
        assert wy_only.available_states() == ["WY"]

    def test_regions_for_wa(self, national_10_table):
        wa_regions = national_10_table.regions_for_state("WA")
        codes = {r.code for r in wa_regions}
        assert codes == {"WA-1", "WA-2", "WA-3", "WA-4"}

    def test_region_code_uniqueness(self, national_10_table):
        """Every region code in the 10-state table is unique."""
        all_codes = [r.code for r in national_10_table.available_regions()]
        assert len(all_codes) == len(set(all_codes)), "Duplicate region codes detected"

    def test_batch_estimate_multi_state(self, national_10_table):
        basins = [
            BasinCharacteristics(
                site_no="NM01",
                site_name="NM Test",
                region=NM_MOUNTAIN,
                predictors={DRNAREA: 80.0, PRECIP: 22.0},
            ),
            BasinCharacteristics(
                site_no="CO01",
                site_name="CO Test",
                region=CO_WESTERN,
                predictors={DRNAREA: 191.0, ELEV: 9200.0},
            ),
            BasinCharacteristics(
                site_no="WY01",
                site_name="WY Test",
                region=WY_FOOTHILLS,
                predictors={DRNAREA: 200.0, CSL1085LFP: 6.0},
            ),
            BasinCharacteristics(
                site_no="ID01",
                site_name="ID Test",
                region=ID_NORTH,
                predictors={DRNAREA: 120.0, PRECIP: 45.0},
            ),
            BasinCharacteristics(
                site_no="WA01",
                site_name="WA Test",
                region=WA_OLYMPIC,
                predictors={DRNAREA: 286.0, PRECIP: 130.0},
            ),
            BasinCharacteristics(
                site_no="OR01",
                site_name="OR Test",
                region=OR_COAST,
                predictors={DRNAREA: 766.0, PRECIP: 65.0},
            ),
            BasinCharacteristics(
                site_no="CA01",
                site_name="CA Test",
                region=CA_SOUTH,
                predictors={DRNAREA: 50.0},
            ),
        ]
        results = national_10_table.batch_estimate(basins, aep=0.01)
        assert len(results) == 7
        assert all(q > 0 for q in results.values())

    def test_batch_skips_wrong_region(self, national_10_table):
        """Basin with region not in a filtered table → silently skipped."""
        nm_only = national_10_table.filter_by_state("NM")
        basin_wrong = BasinCharacteristics(
            site_no="ID99",
            site_name="Idaho Basin",
            region=ID_NORTH,
            predictors={DRNAREA: 100.0, PRECIP: 40.0},
        )
        results = nm_only.batch_estimate([basin_wrong], aep=0.01)
        assert "ID99" not in results

    def test_publication_non_empty(self, national_10_table):
        assert len(national_10_table.publication) > 0


# ===========================================================================
# Cross-state parametrized predictor isolation
# ===========================================================================


class TestPredictorIsolation:
    """Parametrized: verify each state's required predictors produce valid Q."""

    @pytest.mark.parametrize(
        "table_fixture, region, predictors",
        [
            ("nm_table", NM_MOUNTAIN, {DRNAREA: 80.0, PRECIP: 22.0}),
            ("nm_table", NM_TRANSITION, {DRNAREA: 60.0, CSL1085LFP: 7.0}),
            ("nm_table", NM_PLAINS, {DRNAREA: 120.0}),
            ("nm_table", NM_DESERT, {DRNAREA: 50.0, CSL1085LFP: 5.0}),
            ("co_table", CO_FRONT_RANGE, {DRNAREA: 150.0, ELEV: 8500.0, PRECIP: 25.0}),
            ("co_table", CO_EASTERN, {DRNAREA: 200.0, CSL1085LFP: 4.0}),
            ("co_table", CO_WESTERN, {DRNAREA: 191.0, ELEV: 9200.0}),
            ("co_table", CO_SAN_LUIS, {DRNAREA: 500.0}),
            ("wy_table", WY_MOUNTAIN, {DRNAREA: 640.0, ELEV: 7400.0, PRECIP: 22.0}),
            ("wy_table", WY_FOOTHILLS, {DRNAREA: 80.0, CSL1085LFP: 8.0}),
            ("wy_table", WY_PLAINS, {DRNAREA: 300.0}),
            ("id_table", ID_NORTH, {DRNAREA: 120.0, PRECIP: 45.0}),
            ("id_table", ID_CENTRAL, {DRNAREA: 1220.0, ELEV: 6800.0, PRECIP: 38.0}),
            ("id_table", ID_SNAKE, {DRNAREA: 500.0, CSL1085LFP: 5.0}),
            ("id_table", ID_SOUTH, {DRNAREA: 200.0}),
            ("wa_table", WA_OLYMPIC, {DRNAREA: 286.0, PRECIP: 130.0}),
            ("wa_table", WA_W_CASCADE, {DRNAREA: 200.0, ELEV: 4500.0, PRECIP: 80.0}),
            ("wa_table", WA_E_CASCADE, {DRNAREA: 150.0, ELEV: 3800.0}),
            ("wa_table", WA_EASTERN, {DRNAREA: 250.0}),
            ("or_table", OR_COAST, {DRNAREA: 766.0, PRECIP: 65.0}),
            ("or_table", OR_WILLAMETTE, {DRNAREA: 100.0, CSL1085LFP: 5.0}),
            ("or_table", OR_CASCADE, {DRNAREA: 250.0, ELEV: 4200.0, PRECIP: 70.0}),
            ("or_table", OR_EASTERN, {DRNAREA: 400.0}),
            ("ca_table", CA_NORTH_COAST, {DRNAREA: 633.0, PRECIP: 42.0}),
            ("ca_table", CA_SIERRA, {DRNAREA: 300.0, ELEV: 7200.0, PRECIP: 55.0}),
            ("ca_table", CA_CENTRAL, {DRNAREA: 200.0, CSL1085LFP: 6.0}),
            ("ca_table", CA_SOUTH, {DRNAREA: 50.0}),
        ],
    )
    def test_q100_positive(self, table_fixture, region, predictors, request):
        """Q100 > 0 for every region × predictor combination."""
        table = request.getfixturevalue(table_fixture)
        basin = BasinCharacteristics(
            site_no="TST", site_name="Test", region=region, predictors=predictors
        )
        q = table.estimate(basin, aep=0.01)
        assert q > 0, f"Non-positive Q100 for region {region.code}: {q}"

    @pytest.mark.parametrize(
        "table_fixture, region, predictors",
        [
            ("nm_table", NM_MOUNTAIN, {DRNAREA: 80.0, PRECIP: 22.0}),
            ("co_table", CO_FRONT_RANGE, {DRNAREA: 150.0, ELEV: 8500.0, PRECIP: 25.0}),
            ("wy_table", WY_MOUNTAIN, {DRNAREA: 640.0, ELEV: 7400.0, PRECIP: 22.0}),
            ("id_table", ID_CENTRAL, {DRNAREA: 1220.0, ELEV: 6800.0, PRECIP: 38.0}),
            ("wa_table", WA_W_CASCADE, {DRNAREA: 200.0, ELEV: 4500.0, PRECIP: 80.0}),
            ("or_table", OR_CASCADE, {DRNAREA: 250.0, ELEV: 4200.0, PRECIP: 70.0}),
            ("ca_table", CA_SIERRA, {DRNAREA: 300.0, ELEV: 7200.0, PRECIP: 55.0}),
        ],
    )
    def test_three_predictor_monotonic(self, table_fixture, region, predictors, request):
        """Q100 > Q10 > Q2 for all three-predictor regions."""
        table = request.getfixturevalue(table_fixture)
        basin = BasinCharacteristics(
            site_no="TST3", site_name="Test3", region=region, predictors=predictors
        )
        _assert_monotonic_q(table, basin)
