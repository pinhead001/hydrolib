"""Tests for PeakfqSA subprocess wrapper."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from hydrolib.peakfqsa.config import PeakfqSAConfig
from hydrolib.peakfqsa.wrapper import (
    PeakfqSAExecutionError,
    PeakfqSATimeoutError,
    PeakfqSAWrapper,
)

# Mark for tests requiring the real PeakfqSA binary
requires_peakfqsa = pytest.mark.requires_peakfqsa


class TestPeakfqSAWrapper:
    """Tests for PeakfqSAWrapper (mock-based, no real binary needed)."""

    def test_is_available_without_binary(self) -> None:
        """is_available() returns False when binary not found."""
        config = PeakfqSAConfig(executable_path=Path("/nonexistent/peakfqsa"))
        wrapper = PeakfqSAWrapper(config)
        assert wrapper.is_available() is False

    def test_is_available_with_binary(self, tmp_path: Path) -> None:
        """is_available() returns True when binary exists."""
        exe = tmp_path / "peakfqsa"
        exe.touch()
        config = PeakfqSAConfig(executable_path=exe)
        wrapper = PeakfqSAWrapper(config)
        assert wrapper.is_available() is True

    # TODO: Add more mock-based tests as wrapper is implemented


@requires_peakfqsa
class TestPeakfqSAWrapperIntegration:
    """Integration tests requiring the real PeakfqSA binary.

    These tests are skipped when PeakfqSA is not installed.
    """

    # TODO: Add integration tests when wrapper.run() is implemented
    pass
