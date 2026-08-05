"""Tests for hydrolib.regional_skew module."""

import pytest

from hydrolib.regional_skew import (
    NATIONWIDE_FALLBACK,
    NATIONWIDE_SKEW,
    NATIONWIDE_SKEW_SE,
    STATE_SKEW,
    get_regional_skew,
    list_available_states,
    normalize_state,
)


class TestNormalizeState:
    def test_accepts_abbreviation(self):
        assert normalize_state("VT") == "VT"
        assert normalize_state("vt") == "VT"

    def test_accepts_full_name(self):
        assert normalize_state("Vermont") == "VT"
        assert normalize_state("vermont") == "VT"
        assert normalize_state("NORTH CAROLINA") == "NC"

    def test_strips_whitespace(self):
        assert normalize_state("  Vermont  ") == "VT"

    def test_unknown_state_raises(self):
        with pytest.raises(ValueError):
            normalize_state("Atlantis")

    def test_unknown_abbreviation_raises(self):
        with pytest.raises(ValueError):
            normalize_state("ZZ")


class TestGetRegionalSkew:
    def test_known_state_returns_study(self):
        estimate = get_regional_skew("Vermont")
        assert estimate.value == pytest.approx(0.44)
        assert estimate.mse == pytest.approx(0.078)
        assert "VT" in estimate.states

    def test_multi_state_study_shared_across_states(self):
        ga = get_regional_skew("GA")
        nc = get_regional_skew("NC")
        sc = get_regional_skew("SC")
        assert ga is nc is sc
        assert ga.value == pytest.approx(0.048)
        assert ga.mse == pytest.approx(0.092)

    def test_state_deferring_to_nationwide_value_still_documents_source(self):
        sd = get_regional_skew("SD")
        assert sd.value == NATIONWIDE_SKEW
        assert sd.se == NATIONWIDE_SKEW_SE
        assert sd is not NATIONWIDE_FALLBACK
        assert "South Dakota" in sd.source

    def test_unlisted_state_falls_back_to_nationwide(self):
        estimate = get_regional_skew("California")
        assert estimate is NATIONWIDE_FALLBACK
        assert estimate.value == NATIONWIDE_SKEW
        assert estimate.se == NATIONWIDE_SKEW_SE

    def test_montana_uses_nationwide_value_with_recalibrated_se(self):
        mt = get_regional_skew("Montana")
        assert mt.value == NATIONWIDE_SKEW
        assert mt.se == pytest.approx(0.64)
        assert mt.se != NATIONWIDE_SKEW_SE
        assert "Parrett" in mt.source

    def test_washington_defers_to_nationwide_value_and_se(self):
        wa = get_regional_skew("WA")
        assert wa.value == NATIONWIDE_SKEW
        assert wa.se == NATIONWIDE_SKEW_SE
        assert wa is not NATIONWIDE_FALLBACK
        assert "Washington" in wa.source

    def test_idaho_and_oregon_not_yet_covered(self):
        # Both use spatially-varying skew (Idaho: 3 maps by flood type;
        # Oregon: 3 flood regions) that don't fit a single (value, se) --
        # they should fall back to the nationwide default, not a guessed
        # constant.
        idaho = get_regional_skew("Idaho")
        oregon = get_regional_skew("Oregon")
        assert idaho is NATIONWIDE_FALLBACK
        assert oregon is NATIONWIDE_FALLBACK

    def test_arizona_uses_statewide_constant(self):
        az = get_regional_skew("Arizona")
        assert az.value == pytest.approx(-0.09)
        assert az.mse == pytest.approx(0.08)
        assert "Paretti" in az.source

    def test_other_southwest_states_not_yet_covered(self):
        # New Mexico, Nevada, Utah, Colorado, Wyoming, and California were
        # all investigated: no single citable statewide constant could be
        # confirmed, and several (Wyoming: 6 regions; California: skew as a
        # continuous function of elevation) are confirmed to use spatially
        # varying models that don't fit a single (value, se). They should
        # fall back to the nationwide default rather than a guessed number.
        for state in ["New Mexico", "Nevada", "Utah", "Colorado", "Wyoming", "California"]:
            assert get_regional_skew(state) is NATIONWIDE_FALLBACK, state

    def test_missouri_uses_statewide_constant(self):
        mo = get_regional_skew("Missouri")
        assert mo.value == pytest.approx(-0.30)
        assert mo.mse == pytest.approx(0.14)

    def test_maine_uses_statewide_constant(self):
        me = get_regional_skew("Maine")
        assert me.value == pytest.approx(0.029)
        assert me.se == pytest.approx(0.297)

    def test_new_england_multi_state_study_shared_across_states(self):
        ct = get_regional_skew("CT")
        ri = get_regional_skew("RI")
        ma = get_regional_skew("MA")
        assert ct is ri is ma
        assert ct.value == pytest.approx(0.37)
        assert ct.mse == pytest.approx(0.14)

    def test_maine_and_new_hampshire_excluded_from_new_england_constant(self):
        # Maine has its own separate constant (tested above); New Hampshire
        # uses an older spatially-varying map, not the New England SIR
        # 2017-5037 constant -- neither should be swept into the CT/RI/MA
        # study.
        assert get_regional_skew("Maine") is not get_regional_skew("CT")
        assert get_regional_skew("New Hampshire") is NATIONWIDE_FALLBACK

    def test_virginia_west_virginia_multi_state_study_shared(self):
        va = get_regional_skew("VA")
        wv = get_regional_skew("WV")
        assert va is wv
        assert va.value == pytest.approx(0.50)
        assert va.se == pytest.approx(0.574)

    def test_kentucky_and_tennessee_not_swept_into_va_wv_constant(self):
        # Grouped with VA/WV in some earlier USGS work but not confirmed to
        # share this specific constant.
        assert get_regional_skew("Kentucky") is NATIONWIDE_FALLBACK
        assert get_regional_skew("Tennessee") is NATIONWIDE_FALLBACK

    def test_north_dakota_defers_to_nationwide_value_and_se(self):
        nd = get_regional_skew("ND")
        assert nd.value == NATIONWIDE_SKEW
        assert nd.se == NATIONWIDE_SKEW_SE
        assert nd is not NATIONWIDE_FALLBACK
        assert "North Dakota" in nd.source

    def test_arkansas_and_hawaii_not_added_without_confirmed_se(self):
        # Both have a confirmed skew value (-0.17 and -0.14 respectively)
        # but no confirmed standard error/MSE was found, so they are
        # deliberately left out of the code dataset rather than paired
        # with a guessed uncertainty.
        assert get_regional_skew("Arkansas") is NATIONWIDE_FALLBACK
        assert get_regional_skew("Hawaii") is NATIONWIDE_FALLBACK

    def test_invalid_state_raises(self):
        with pytest.raises(ValueError):
            get_regional_skew("Not A State")

    def test_mse_is_se_squared(self):
        estimate = get_regional_skew("Vermont")
        assert estimate.mse == pytest.approx(estimate.se**2)


class TestListAvailableStates:
    def test_returns_sorted_codes(self):
        states = list_available_states()
        assert states == sorted(states)
        assert "VT" in states
        assert set(states) == set(STATE_SKEW.keys())

    def test_matches_nationwide_default_constant(self):
        assert NATIONWIDE_SKEW == -0.302
        assert NATIONWIDE_SKEW_SE == 0.55
