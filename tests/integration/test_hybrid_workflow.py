"""Integration tests for hybrid Bulletin 17C workflow."""

from __future__ import annotations

import pytest

# Mark for tests requiring the real PeakfqSA binary
requires_peakfqsa = pytest.mark.requires_peakfqsa


class TestBigSandyNativeOnly:
    """Test Big Sandy analysis using native EMA only (no PeakfqSA needed)."""

    # TODO: Implement native-only quantile validation test
    pass


@requires_peakfqsa
class TestBigSandyEndToEnd:
    """Full integration: native + PeakfqSA comparison.

    Skipped when PeakfqSA is not installed.
    """

    # TODO: Implement end-to-end comparison test
    pass
