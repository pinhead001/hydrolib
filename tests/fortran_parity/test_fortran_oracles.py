"""Per-routine oracles from the vendored Fortran.

``emafitpr`` alone is a poor way to check a ported routine. The EMA fit is a
fixed point whose condition number is around 1e13 -- one ulp of input moves the
converged at-site skew by 3e-3 -- so a correct routine can look wrong through
it and a wrong one can look right. ``build_fortran/_emafort.pyf`` therefore
also exposes ``mseg_all_sub``, ``detratsub``, ``var_mom`` and ``moms_p3``, all
of which peakfqr's own R code already calls through ``.Fortran()``.

These tests pin what those oracles say, so the ``var_mom`` port (TODO.md P3)
can be written against them one routine at a time instead of end to end.

They also record two things the oracles found immediately:

* ``_b17b_skew_mse`` matches the Fortran's ``mseg()`` exactly up to n = 150 and
  diverges above it, because ADJE's bias adjustment partially undoes the cap.
* ``_ema_iteration`` diverged from ``moms_p3`` on censored rows only, which was
  first read as a defect in the truncated-moment code and turned out to be the
  bias-correction convention -- ``bcf = 1997`` builds c2 and c3 from the
  exact-peak count. Fixed; both tests below are now equality tests, and the
  contrast between them is what localised it.
"""

from __future__ import annotations

import collections

import numpy as np
import pytest

from tests.fortran_parity.conftest import load_golden

pytest.importorskip(
    "hydrolib.peakfqr",
    reason="Fortran extension not built; run python build_fortran/build.py",
)

pytestmark = pytest.mark.requires_fortran

CASES = ("big_sandy_03606500", "powder_river_06326500", "cains_coulee_06327450")

QMIN, QMAX = 1e-20, 1e20


def _threshold_groups(tl, tu):
    """Distinct (tl, tu) perception-threshold pairs and their counts.

    These routines take threshold *groups*, not per-observation arrays:
    ``nobs(i)`` counts the observations sharing the pair ``(tl(i), tu(i))``.

    Grouped on exact values, deliberately. Rounding to 12 decimals to be
    "robust" perturbs log10(18000) enough to move Big Sandy's at-site skew MSE
    by 2.2e-4 -- the fixed point amplifies a 1e-12 input change that far. The
    values here are exact duplicates already, so there is nothing to round.
    """
    counter = collections.Counter(zip(tl, tu))
    return (
        np.array(list(counter.values()), dtype=float),
        np.array([k[0] for k in counter], dtype=float),
        np.array([k[1] for k in counter], dtype=float),
    )


def _golden(name):
    golden = load_golden(name)
    if golden is None:
        pytest.skip(f"golden file missing for {name}")
    return golden


def _at_site_moments(golden):
    """Column 2 of cmoms: the at-site fit, which is what these routines take."""
    cmoms = np.asarray(golden["outputs"]["cmoms"], dtype=float)
    return np.array([cmoms[0][1], cmoms[1][1], cmoms[2][1]])


class TestSkewMseOracle:
    """mseg_all_sub: the at-site skew MSE, through the default ADJE option."""

    @pytest.mark.parametrize("name", CASES)
    def test_reproduces_emafitpr_as_g_mse(self, name):
        """Called directly it must give exactly what the full fit reports.

        This is the oracle the ADJE port is written against: it already
        includes the censoring bias adjustment hydrolib does not implement --
        0.0944 on Big Sandy against an uncensored 0.0636.
        """
        from hydrolib.peakfqr._emafort import mseg_all_sub

        golden = _golden(name)
        nobs, tl, tu = _threshold_groups(golden["inputs"]["tl"], golden["inputs"]["tu"])
        direct = float(mseg_all_sub(nobs, tl, tu, _at_site_moments(golden)))
        assert direct == pytest.approx(golden["outputs"]["skew"]["as_G_mse_o"], rel=1e-12)

    @pytest.mark.parametrize("n, skew", [(84, 0.0066), (85, -0.181), (32, -0.708), (50, 0.5)])
    def test_matches_b17b_formula_when_nothing_is_censored(self, n, skew):
        """With no censoring ADJE's bias adjustment is 1, so it reduces to mseg().

        That makes this a direct check of hydrolib's ``_b17b_skew_mse``.
        """
        from hydrolib.bulletin17c import _b17b_skew_mse
        from hydrolib.peakfqr._emafort import mseg_all_sub

        fortran = float(
            mseg_all_sub(
                np.array([float(n)]),
                np.array([-99.0]),
                np.array([99.0]),
                np.array([0.0, 0.09, skew]),
            )
        )
        assert fortran == pytest.approx(_b17b_skew_mse(n, skew), rel=1e-6)

    def test_b17b_formula_diverges_above_150_observations(self):
        """A real limitation, recorded rather than hidden.

        ``mseg_all`` evaluates ``mseg()`` at ``min(n, 150)`` and then multiplies
        by a bias adjustment that partially undoes the cap.
        ``_b17b_skew_mse`` applies the cap without the adjustment, so above 150
        it overestimates -- 0.0479 against the Fortran's 0.0365 at n = 200, 24%
        high, which over-weights the regional skew on a long record. Fixing it
        needs ``mse_ema``, hence ``var_mom``; see TODO.md P3. No parity case
        reaches n = 150, so nothing else detects this.
        """
        from hydrolib.bulletin17c import _b17b_skew_mse
        from hydrolib.peakfqr._emafort import mseg_all_sub

        mc = np.array([0.0, 0.09, 0.3])
        fortran = float(mseg_all_sub(np.array([200.0]), np.array([-99.0]), np.array([99.0]), mc))
        assert fortran == pytest.approx(0.03648, abs=1e-5)
        assert _b17b_skew_mse(200, 0.3) == pytest.approx(0.04789, abs=1e-5)
        assert _b17b_skew_mse(200, 0.3) / fortran == pytest.approx(1.31, abs=0.02)


class TestDeterminantRatioOracle:
    """detratsub: the Halloween Wd, which hydrolib does not implement."""

    def test_needs_the_post_mgbt_thresholds(self):
        """The usage note that matters, and it is not obvious.

        Cains Coulee's *input* thresholds are uncensored throughout, so calling
        detrat with them returns 1.0. Its censoring is created by MGBT inside
        the fit, which raises the lower threshold to log10(332) = 2.521, and
        only with those does detrat give the 0.184 that emafitpr reports.
        emafitpr returns them as tlema/tuema.
        """
        from hydrolib.peakfqr._emafort import detratsub, emafitpr
        from tests.fortran_parity.cases import CASES as CASE_FACTORIES
        from tests.fortran_parity.cases import build_emafit_inputs

        args = build_emafit_inputs(CASE_FACTORIES["cains_coulee_06327450"]())
        out = emafitpr(
            args["ql"],
            args["qu"],
            args["tl"],
            args["tu"],
            args["dtype"],
            args["reg_m"],
            args["reg_m_mse"],
            args["reg_sd"],
            args["reg_sd_mse"],
            args["r_g"],
            args["r_g_mse"],
            args["gbthrsh0"],
            args["pq"],
            args["eps"],
            args["wght_opt_n"],
        )
        cmoms, wdout = np.asarray(out[9], dtype=float), float(out[14])
        n = len(args["ql"])
        tlema, tuema = np.asarray(out[17], dtype=float), np.asarray(out[18], dtype=float)
        mc = np.array([cmoms[0][1], cmoms[1][1], cmoms[2][1]])

        as_supplied = float(detratsub(mc, n, *_threshold_groups(args["tl"], args["tu"])))
        assert as_supplied == pytest.approx(1.0), "input thresholds carry no censoring"

        after_mgbt = float(detratsub(mc, n, *_threshold_groups(tlema[:n], tuema[:n])))
        assert after_mgbt == pytest.approx(wdout, rel=1e-9)
        assert after_mgbt == pytest.approx(0.184, abs=5e-4)


class TestMomentIterationOracle:
    """moms_p3 against _ema_iteration, the transcription it was written from."""

    @staticmethod
    def _fit(name):
        from hydrolib.bulletin17c import Bulletin17C
        from tests.fortran_parity.cases import CASES as CASE_FACTORIES

        case = CASE_FACTORIES[name]()
        b17c = Bulletin17C(
            peak_flows=list(case.systematic.values()) + list(case.historical.values()),
            water_years=list(case.systematic) + list(case.historical),
            regional_skew=case.regional_skew,
            regional_skew_mse=case.regional_skew_mse,
        )
        b17c.run_analysis(method="ema")
        return b17c

    @staticmethod
    def _same_rows(analyzer):
        """The exact intervals _ema_iteration is working on, as ql/qu.

        Both sides must be handed the same rows. MGBT censors peaks *inside*
        run_analysis, so the case's original ql/qu are not what the native
        iteration sees.
        """
        ql, qu = [], []
        for interval in analyzer.intervals:
            low = interval.lower if interval.lower > 0 else QMIN
            high = interval.upper if np.isfinite(interval.upper) else QMAX
            ql.append(np.log10(low))
            qu.append(np.log10(low if not interval.is_censored else high))
        return np.array(ql), np.array(qu)

    @staticmethod
    def _compare(b17c):
        from hydrolib.peakfqr._emafort import moms_p3

        analyzer = b17c._analyzer
        ql, qu = TestMomentIterationOracle._same_rows(analyzer)
        results = b17c.results
        start = np.array([results.mean_log, results.std_log**2, results.skew_station])
        fortran = np.asarray(moms_p3(ql, qu, 0.0, 0.0, start), dtype=float)
        mean, std, skew = analyzer._ema_iteration(
            start[0], np.sqrt(start[1]), start[2], n_regional=0.0
        )
        return fortran, np.array([mean, std**2, skew]), int((ql != qu).sum())

    @pytest.mark.parametrize("name", ["powder_river_06326500", "big_sandy_03606500"])
    def test_exact_on_uncensored_rows(self, name):
        """Where no interval is censored the transcription is exact.

        Measured differences: 0.0 on the mean, ~1e-14 on the variance and
        ~1e-12 on the skew. That is the formulas in moms_p3 being right.
        """
        fortran, mine, n_censored = self._compare(self._fit(name))
        assert n_censored == 0, "this case should have no censored intervals here"
        assert np.allclose(fortran, mine, rtol=0.0, atol=1e-9)

    def test_exact_on_censored_rows_too(self):
        """Cains Coulee censors 11 of 32, and the two agree there as well.

        This test used to assert the opposite -- a 0.70% variance gap and a
        4.94% skew gap -- and concluded from it that hydrolib's truncated-P3
        moment code was wrong on censored intervals. It was not. ``moms_p3``
        reads its bias-correction convention from ``common /tac002/``, whose
        blockdata default is ``data bcf/1997/`` (``emafit.f:3898``): c2 and c3
        are built from the *exact-peak* count, not the total. ``n_e == n``
        whenever nothing is censored, which is why the uncensored case above
        passed throughout and why nothing here could see the difference.

        With :func:`_bias_correction_factors` fixed, the measured differences
        are ~7e-11 on the mean, 1.7e-10 on the variance and ~2e-10 on the
        skew. The bound below is 1e-8.

        What that leaves for the ``var_mom`` port is the skew *weighting*
        path -- the ADJE bias adjustment and the determinant ratio -- not the
        moment code. See ``docs/VAR_MOM_PORT_PLAN.md``.
        """
        fortran, mine, n_censored = self._compare(self._fit("cains_coulee_06327450"))
        assert n_censored == 11, "this case should censor 11 of its 32 rows"
        assert np.allclose(fortran, mine, rtol=0.0, atol=1e-8), f"{fortran!r} vs {mine!r}"


class TestVarianceOfMomentsOracle:
    """var_mom: the root of the dependency tree the whole port hangs off."""

    @pytest.mark.parametrize("name", CASES)
    def test_returns_a_symmetric_positive_definite_matrix(self, name):
        """No independent reference for this one, so check what must hold.

        It is a covariance matrix of (mean, variance, skew) estimators, so it
        must be 3x3, symmetric and positive definite. A port that produces
        anything else is wrong regardless of the numbers.
        """
        from hydrolib.peakfqr._emafort import var_mom

        golden = _golden(name)
        nobs, tl, tu = _threshold_groups(golden["inputs"]["tl"], golden["inputs"]["tu"])
        varm = np.asarray(var_mom(nobs, tl, tu, _at_site_moments(golden)), dtype=float)
        assert varm.shape == (3, 3)
        assert np.allclose(varm, varm.T, rtol=1e-9), "covariance must be symmetric"
        assert np.all(np.linalg.eigvalsh(varm) > 0), "covariance must be positive definite"
        assert np.all(np.diag(varm) > 0)
