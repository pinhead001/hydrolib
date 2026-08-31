"""Native EMA against peakfq 8.1.0 on two Wyoming/Montana sites.

Big Sandy was the only parity site for a long time, and it cannot answer two
questions on its own.

Does the fit agree when nothing is censored?
    Powder River is 85 contiguous water years with no historic peaks, no zero
    flows and no PILFs -- a plain systematic record fitted with a weighted
    regional skew. It answers yes, to machine precision, which is what says
    the in-loop regional-skew weighting is right rather than merely closer.

What does censoring cost?
    Cains Coulee is 32 contiguous years whose MGBT finds 11 PILFs, so the
    censoring is produced by the fit rather than supplied as input. Its
    reference ``Wd`` is 0.184 -- the first case in this repository where the
    Halloween determinant ratio is far from 1, since ``emafit.f:763`` only
    reaches ``detrat`` when the at-site skew clears 0.04 and Big Sandy's is
    0.0066. hydrolib has no ``detrat``, so it uses 1.

Read together they localise the open P3 defect precisely: with no censoring
the native fit is exact, and everything that remains is censored-path work.
That contrast has already paid for itself once. Cains Coulee's at-site skew
was 0.122 off and is now 2e-5 off, because the pair of cases isolated the
cause to the moment bias-correction convention -- ``moms_p3`` builds c2 and c3
from the exact-peak count (``bcf = 1997``, ``emafit.f:3898``) and hydrolib used
the total interval count, which is invisible on an uncensored record. What is
left is the skew *weighting*: the ADJE bias adjustment on the at-site skew MSE,
and ``detrat`` itself.

The peakfq 7.4 columns in the fixture CSVs are a sanity cross-check only.
Parity is against the committed goldens, generated from the vendored 8.1.0
Fortran by ``tools/gen_fortran_golden.py``.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.fixtures.paths import SKIP_REASON, TESTDATA_AVAILABLE
from tests.fortran_parity.conftest import load_golden

# The marker alone does not skip anything -- it only labels. Without the
# skipif the module-scope fixtures raise FileNotFoundError and every test in
# here reports as an ERROR rather than a skip when the reference tree is
# absent. Same pairing as tests/test_r_fixtures.py.
pytestmark = [
    pytest.mark.requires_peakfqr_testdata,
    pytest.mark.skipif(not TESTDATA_AVAILABLE, reason=SKIP_REASON),
]


def _native(site_no: str):
    from hydrolib.bulletin17c import Bulletin17C
    from tests.fixtures.wymt_peaks import load_site

    site = load_site(site_no)
    years = sorted(site.peaks)
    b17c = Bulletin17C(
        peak_flows=[site.peaks[y] for y in years],
        water_years=years,
        regional_skew=site.regional_skew,
        regional_skew_mse=site.regional_skew_mse,
    )
    b17c.run_analysis(method="ema")
    return b17c.results


def _reference(name: str):
    golden = load_golden(name)
    if golden is None:
        pytest.skip(f"golden file missing for {name} (run tools/gen_fortran_golden.py)")
    cmoms = np.asarray(golden["outputs"]["cmoms"], dtype=float)
    return {
        "mean_log": cmoms[0][0],
        "std_log": float(np.sqrt(cmoms[1][0])),
        "skew_weighted": cmoms[2][0],
        "skew_at_site": cmoms[2][1],
        "aeps": np.asarray(golden["inputs"]["aeps"], dtype=float),
        "quantiles": 10.0 ** np.asarray(golden["outputs"]["quantiles"]["yp"], dtype=float),
        "mgbt": golden["outputs"]["mgbt"],
        "wd": golden["outputs"]["skew"]["Wdout"],
    }


def _quantile_errors(results, ref) -> dict:
    """Percent error per AEP, for the AEPs the native fit reports."""
    q = results.quantiles
    errors = {}
    for aep, ref_q in zip(ref["aeps"], ref["quantiles"]):
        row = q[np.isclose(q["aep"], aep)]
        if not row.empty:
            errors[float(aep)] = abs(float(row["flow_cfs"].iloc[0]) / ref_q - 1) * 100
    return errors


@pytest.fixture(scope="module")
def powder_river():
    """(native results, reference) for USGS 06326500."""
    return _native("06326500.00"), _reference("powder_river_06326500")


@pytest.fixture(scope="module")
def cains_coulee():
    """(native results, reference) for USGS 06327450."""
    return _native("06327450.00"), _reference("cains_coulee_06327450")


class TestPowderRiverUncensored:
    """USGS 06326500, 1938-2022. No censoring anywhere, and it agrees exactly.

    This is the test that makes the weighting fix a claim rather than a
    direction: on a systematic record with a weighted regional skew, the native
    EMA reproduces peakfq 8.1.0 to machine precision. Measured absolute
    differences are 0.0 on the mean, 3.7e-14 on the standard deviation,
    4.5e-12 on the at-site skew and 7.5e-11 on the weighted skew. The bound
    below is 1e-6 -- four orders of magnitude looser than that, and four
    tighter than anything the open defect could hide in.
    """

    @pytest.mark.parametrize("field", ["mean_log", "std_log", "skew_weighted", "skew_at_site"])
    def test_moment_matches_to_machine_precision(self, powder_river, field):
        results, ref = powder_river
        native = {
            "mean_log": results.mean_log,
            "std_log": results.std_log,
            "skew_weighted": results.skew_weighted,
            "skew_at_site": results.skew_station,
        }[field]
        assert (
            abs(native - ref[field]) < 1e-6
        ), f"{field}: native {native!r} vs peakfq {ref[field]!r}"

    def test_no_pilfs_either_side(self, powder_river):
        """MGBT must agree that there is nothing to censor."""
        results, ref = powder_river
        assert results.n_low_outliers == ref["mgbt"]["gbnlow"] == 0

    def test_determinant_ratio_is_one_without_censoring(self, powder_river):
        """Wd = 1 here, so the missing detrat cannot be what is being measured."""
        _, ref = powder_river
        assert ref["wd"] == pytest.approx(1.0)

    def test_quantiles_agree_within_half_a_percent(self, powder_river):
        results, ref = powder_river
        errors = _quantile_errors(results, ref)
        worst = max(errors, key=errors.get)
        assert errors[worst] < 0.5, f"worst at AEP {worst}: {errors[worst]:.3f}%"


class TestCainsCouleeCensored:
    """USGS 06327450, 1991-2022. MGBT censors 11 peaks, and the skew weighting diverges.

    Everything discrete matches: the same 11 PILFs at the same 332 cfs cut.
    So does the at-site fit now -- its skew agrees to 2e-5, since the moment
    bias-correction factors were corrected to the exact-peak count.

    What does not match is the *weighted* skew, and the reference ``Wd`` of
    0.184 is most of why: hydrolib has no ``detrat`` and uses 1, so it weights
    the regional skew more than fivefold too heavily here. The at-site skew
    MSE's ADJE adjustment is the rest.
    """

    def test_mgbt_finds_the_same_pilfs(self, cains_coulee):
        """MGBT is the part verified line-by-line against the Fortran; it holds here."""
        results, ref = cains_coulee
        assert results.n_low_outliers == ref["mgbt"]["gbnlow"] == 11
        assert results.low_outlier_threshold == pytest.approx(332.0, rel=1e-3)

    def test_reference_exercises_the_determinant_ratio(self, cains_coulee):
        """The point of this case: Wd is nowhere near 1, unlike every other case."""
        _, ref = cains_coulee
        assert ref["wd"] == pytest.approx(0.184, abs=5e-4)

    def test_mean_still_agrees_closely(self, cains_coulee):
        """Censoring moves the mean least; measured 5.9e-3 in log10, ~1.4% in flow."""
        results, ref = cains_coulee
        assert abs(results.mean_log - ref["mean_log"]) < 0.01

    def test_at_site_skew_matches(self, cains_coulee):
        """Was xfail at a gap of 0.122; now agrees to 2e-5.

        The gap was the moment bias-correction convention, not the censored
        moment code: ``moms_p3`` builds c2 and c3 from the exact-peak count
        (``bcf = 1997``, ``emafit.f:3898``) and hydrolib used the total
        interval count. Since ``n_e == n`` without censoring, this 11-PILF
        record is the only parity case that could ever have shown it.

        The *weighted* skew is a separate defect and is still xfail below.
        """
        results, ref = cains_coulee
        assert abs(results.skew_station - ref["skew_at_site"]) < 0.02

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "Skew-weighting defect, and the last one on this case. peakfq's "
            "at-site skew MSE uses the ADJE censoring bias adjustment (from "
            "var_mom, not ported) and its Wd comes from detrat (also not "
            "ported, and 0.184 here against the 1 hydrolib assumes). Measured "
            "gap 0.124 on the weighted skew, against an at-site skew that now "
            "matches to 2e-5. See docs/VAR_MOM_PORT_PLAN.md phases 4 and 5."
        ),
    )
    def test_weighted_skew_matches(self, cains_coulee):
        results, ref = cains_coulee
        assert abs(results.skew_weighted - ref["skew_weighted"]) < 0.02

    def test_quantile_error_is_bounded_and_worst_in_the_lower_tail(self, cains_coulee):
        """Recorded, not asserted away: 21.07% at worst (AEP 0.995), 4.69% at Q100.

        The lower tail is where the censoring bites, which is the signature of
        the open defect rather than of a broken fit.

        Both moved the wrong way with the bias-correction fix -- 18.7% and
        2.7% before it -- and both are downstream of the weighted skew, which
        is still 0.124 off while the at-site skew is now 2e-5 off. Correcting
        the at-site skew moves ``nG`` through ``_b17b_skew_mse``, so a fit that
        is right about the at-site moments is *further* from peakfq's quantiles
        until ``detrat`` and ADJE land. The bounds are unchanged rather than
        widened; at 1.19x and 1.07x of headroom they will alarm on anything
        else that shifts, which is what they are for.
        """
        results, ref = cains_coulee
        errors = _quantile_errors(results, ref)
        assert max(errors.values()) < 25.0
        assert errors[0.01] < 5.0
        worst_aep = max(errors, key=errors.get)
        assert worst_aep > 0.5, f"worst error at AEP {worst_aep}, expected the lower tail"
