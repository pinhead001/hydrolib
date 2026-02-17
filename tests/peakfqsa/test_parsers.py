"""Tests for PeakfqSA output parser."""

from __future__ import annotations

import pytest

from hydrolib.peakfqsa.parsers import PeakfqSAResult, parse_peakfqsa_output


class TestPeakfqSAResult:
    """Tests for the PeakfqSAResult dataclass."""

    def test_default_construction(self) -> None:
        """Default result has empty/zero fields."""
        result = PeakfqSAResult()
        assert result.station_name == ""
        assert result.n_peaks == 0
        assert result.quantiles == {}
        assert result.raw_output == ""

    def test_raw_output_stored(self) -> None:
        """raw_output field preserves full text."""
        text = "some output text"
        result = PeakfqSAResult(raw_output=text)
        assert result.raw_output == text


class TestParseOutput:
    """Tests for parse_peakfqsa_output function."""

    # TODO: Add parser tests with fixture output text
    pass
