"""Tests for multi-site batch processing.

``hydrolib.batch`` had no tests at all, and that is how ``analyze_sites``
shipped broken: ``fetch_nwis_batch`` returns plain dicts, ``B17CEngine.fit``
reads ``record.flow``, and ``run_multi_site``'s ``except Exception`` turned the
resulting ``AttributeError`` into ``{"error": ...}`` for every site. The
function reported failure for every input and never said why.

Network access is mocked throughout; no test in this module contacts NWIS.
"""

from __future__ import annotations

from unittest.mock import patch

import pandas as pd

from hydrolib.batch import analyze_sites, batch_summary_table, run_multi_site
from hydrolib.core import PeakRecord

# Enough spread to give a finite skew and a well-conditioned LP3 fit.
_FLOWS = [12000.0, 8500.0, 15200.0, 9800.0, 22000.0, 7400.0, 11300.0, 18600.0]


def _records(flows=_FLOWS, start_year=1990):
    return [PeakRecord(year=start_year + i, flow=q) for i, q in enumerate(flows)]


def _dicts(flows=_FLOWS, start_year=1990):
    """The shape ``fetch_nwis_batch`` actually returns."""
    return [{"year": start_year + i, "flow": q, "source": "USGS"} for i, q in enumerate(flows)]


class TestRunMultiSite:
    def test_fits_each_site(self):
        results = run_multi_site({"01010000": _records(), "02020000": _records()})

        assert set(results) == {"01010000", "02020000"}
        for site in results:
            assert "error" not in results[site]
            mu, sigma, skew = results[site]["params"]
            assert 3.0 < mu < 4.5
            assert sigma > 0
            assert results[site]["n"] == len(_FLOWS)

    def test_quantiles_increase_with_return_period(self):
        results = run_multi_site({"01010000": _records()}, return_periods=(2, 10, 100))
        quantiles = results["01010000"]["quantiles"]

        assert quantiles[2] < quantiles[10] < quantiles[100]

    def test_confidence_interval_brackets_the_estimate(self):
        results = run_multi_site({"01010000": _records()}, return_periods=(100,))
        ci = results["01010000"]["ci"][100]

        assert ci["lower"] < ci["estimate"] < ci["upper"]

    def test_one_bad_site_does_not_stop_the_others(self):
        # Two records is below B17CEngine.fit's minimum of three.
        data = {"good": _records(), "bad": _records(flows=[100.0, 200.0])}
        results = run_multi_site(data)

        assert "error" not in results["good"]
        assert "error" in results["bad"]
        assert "3" in results["bad"]["error"]


class TestAnalyzeSites:
    def test_converts_fetched_dicts_into_peak_records(self):
        """The regression test for the defect described in this module's docstring."""
        fetched = ({"01010000": _dicts()}, {})

        with patch("hydrolib.batch.fetch_nwis_batch", return_value=fetched):
            results, fetch_errors = analyze_sites(["01010000"])

        assert fetch_errors == {}
        assert "error" not in results["01010000"], results["01010000"]
        assert results["01010000"]["n"] == len(_FLOWS)

    def test_matches_a_direct_run_multi_site_call(self):
        """Going through the fetch path must not change the numbers."""
        fetched = ({"01010000": _dicts()}, {})

        with patch("hydrolib.batch.fetch_nwis_batch", return_value=fetched):
            via_fetch, _ = analyze_sites(["01010000"], return_periods=(10, 100))
        direct = run_multi_site({"01010000": _records()}, return_periods=(10, 100))

        assert via_fetch["01010000"]["params"] == direct["01010000"]["params"]
        assert via_fetch["01010000"]["quantiles"] == direct["01010000"]["quantiles"]

    def test_censored_flow_survives_the_conversion(self):
        rows = _dicts() + [{"year": 1998, "flow": None, "source": "USGS"}]
        fetched = ({"01010000": rows}, {})

        with patch("hydrolib.batch.fetch_nwis_batch", return_value=fetched):
            results, _ = analyze_sites(["01010000"])

        # fit() drops the null, so n is the count of real observations.
        assert results["01010000"]["n"] == len(_FLOWS)

    def test_fetch_errors_pass_through_untouched(self):
        fetched = ({"01010000": _dicts()}, {"09999999": "HTTP 404"})

        with patch("hydrolib.batch.fetch_nwis_batch", return_value=fetched):
            results, fetch_errors = analyze_sites(["01010000", "09999999"])

        assert fetch_errors == {"09999999": "HTTP 404"}
        assert "09999999" not in results

    def test_workers_reaches_the_fetcher(self):
        with patch("hydrolib.batch.fetch_nwis_batch", return_value=({}, {})) as fetch:
            analyze_sites(["01010000"], workers=3)

        assert fetch.call_args.kwargs["workers"] == 3


class TestBatchSummaryTable:
    def test_one_row_per_site(self):
        results = run_multi_site({"01010000": _records(), "02020000": _records()})
        table = batch_summary_table(results, return_periods=(10, 100))

        assert len(table) == 2
        assert set(table["Site"]) == {"01010000", "02020000"}
        assert list(table.columns[:5]) == ["Site", "n", "Mean (log)", "Std (log)", "Skew"]

    def test_carries_a_column_per_requested_return_period(self):
        results = run_multi_site({"01010000": _records()})
        table = batch_summary_table(results, return_periods=(10, 100))

        assert "Q10" in table.columns
        assert "Q100" in table.columns
        assert "Q50" not in table.columns
        assert table.loc[0, "Q10"] < table.loc[0, "Q100"]

    def test_a_return_period_the_fit_did_not_compute_is_skipped(self):
        # run_multi_site fitted the standard set, which has no T = 1000.
        results = run_multi_site({"01010000": _records()})
        table = batch_summary_table(results, return_periods=(100, 1000))

        assert "Q100" in table.columns
        assert "Q1000" not in table.columns

    def test_failed_sites_appear_with_an_error_instead_of_quantiles(self):
        data = {"good": _records(), "bad": _records(flows=[100.0, 200.0])}
        table = batch_summary_table(run_multi_site(data))

        assert set(table["Site"]) == {"good", "bad"}
        bad = table.set_index("Site").loc["bad"]
        assert isinstance(bad["Error"], str) and bad["Error"]
        assert pd.isna(bad["Q100"])

    def test_empty_results_give_an_empty_table(self):
        table = batch_summary_table({})

        assert len(table) == 0
