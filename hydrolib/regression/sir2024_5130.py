"""
hydrolib.regression.sir2024_5130 - Tennessee-specific helpers for SIR 2024-5130.

Reference
---------
Ladd, D.E., and Ensminger, P.A., 2025, Estimating the magnitude and
frequency of floods at ungaged locations on streams in Tennessee through
the 2013 water year: U.S. Geological Survey Scientific Investigations
Report 2024-5130, 19 p., https://doi.org/10.3133/sir20245130.

Contents
--------
* ``TN_AREA1`` … ``TN_AREA4`` — pre-built :class:`HydrologicRegion` instances
  for the four Tennessee hydrologic areas defined in SIR 2024-5130.
* ``HydrologicArea`` — backward-compatible class that exposes the four
  regions as class attributes (``HydrologicArea.AREA1``, etc.) and
  provides ``from_int()``.
* ``SIR2024_5130`` — thin :class:`RegressionTable` subclass pre-configured
  with the correct publication string and a ``write_table4_template()``
  helper.

Populating Coefficients
-----------------------
Table 4 values are **not hardcoded** because the full report must be
consulted directly.  Populate via::

    table = SIR2024_5130.load_from_dict(coeff_dict)
    table = SIR2024_5130.load_from_csv("sir2024_5130_table4.csv")
    SIR2024_5130.write_table4_template("table4_template.csv")

Legacy dict format (b0/b1/b2) is still accepted::

    coeff_dict = {
        (HydrologicArea.AREA2, 0.01): {
            "b0": 2.90, "b1": 0.75, "b2": 0.38,
            "sep_pct": 34.0, "pseudo_r2": 0.93,
        },
    }
    table = SIR2024_5130.load_from_dict(coeff_dict)

New-style dict (preferred)::

    coeff_dict = {
        (HydrologicArea.AREA2, 0.01): {
            "intercept": 2.90,
            "coefficients": {"DRNAREA": 0.75, "CSL1085LFP": 0.38},
            "sep_pct": 34.0, "pseudo_r2": 0.93,
        },
    }
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from hydrolib.regression.region import HydrologicRegion, RegionRegistry
from hydrolib.regression.regression_table import (
    STANDARD_AEPS,
    STANDARD_RETURN_PERIODS,
    GlsEquation,
    RegressionTable,
    _dict_to_equation,
)

log = logging.getLogger(__name__)

_PUB = "SIR 2024-5130; https://doi.org/10.3133/sir20245130"

# ---------------------------------------------------------------------------
# Pre-defined Tennessee hydrologic regions (SIR 2024-5130)
# ---------------------------------------------------------------------------

TN_AREA1 = HydrologicRegion(
    code="TN_AREA1",
    label="Tennessee Hydrologic Area 1 — West TN Lowlands",
    state="TN",
    required_predictors=("DRNAREA", "I2"),
    publication=_PUB,
    notes="Predictor I2 = 2-year recurrence-interval precipitation intensity (in/hr). "
    "Applicable to non-urban, unregulated streams in west Tennessee.",
)

TN_AREA2 = HydrologicRegion(
    code="TN_AREA2",
    label="Tennessee Hydrologic Area 2 — Middle TN",
    state="TN",
    required_predictors=("DRNAREA", "CSL1085LFP"),
    publication=_PUB,
    notes="Predictor CSL1085LFP = 10-85 channel slope (ft/mi). "
    "Applicable to non-urban, unregulated streams in middle Tennessee.",
)

TN_AREA3 = HydrologicRegion(
    code="TN_AREA3",
    label="Tennessee Hydrologic Area 3 — Upper Cumberland / East TN Valley",
    state="TN",
    required_predictors=("DRNAREA",),
    publication=_PUB,
    notes="Drainage area (sq mi) is the sole predictor. "
    "Applicable to non-urban, unregulated streams in upper Cumberland and "
    "east Tennessee valley.",
)

TN_AREA4 = HydrologicRegion(
    code="TN_AREA4",
    label="Tennessee Hydrologic Area 4 — Ridge and Valley / Blue Ridge",
    state="TN",
    required_predictors=("DRNAREA", "IMPERV"),
    publication=_PUB,
    notes="Predictor IMPERV = percent impervious cover (%). "
    "Applicable to non-urban, unregulated streams in Ridge and Valley "
    "and Blue Ridge physiographic provinces of Tennessee.",
)

# Ordered tuple for iteration
TN_REGIONS: Tuple[HydrologicRegion, ...] = (TN_AREA1, TN_AREA2, TN_AREA3, TN_AREA4)


# ---------------------------------------------------------------------------
# HydrologicArea — backward-compatible class attribute accessor
# ---------------------------------------------------------------------------


class HydrologicArea:
    """
    Backward-compatible accessor for the four SIR 2024-5130 hydrologic areas.

    Attributes
    ----------
    AREA1, AREA2, AREA3, AREA4 : HydrologicRegion
        Pre-built region objects for Tennessee areas 1–4.

    Examples
    --------
    >>> from hydrolib.regression.sir2024_5130 import HydrologicArea
    >>> region = HydrologicArea.AREA2
    >>> region.code
    'TN_AREA2'
    >>> HydrologicArea.from_int(2) is HydrologicArea.AREA2
    True
    """

    AREA1: HydrologicRegion = TN_AREA1
    AREA2: HydrologicRegion = TN_AREA2
    AREA3: HydrologicRegion = TN_AREA3
    AREA4: HydrologicRegion = TN_AREA4

    @classmethod
    def from_int(cls, n: int) -> HydrologicRegion:
        """Return the area region for integer 1–4."""
        mapping = {1: cls.AREA1, 2: cls.AREA2, 3: cls.AREA3, 4: cls.AREA4}
        if n not in mapping:
            raise ValueError(f"Hydrologic area must be 1-4, got {n}")
        return mapping[n]

    @classmethod
    def all_regions(cls) -> List[HydrologicRegion]:
        """Return all four TN regions as a list."""
        return list(TN_REGIONS)


# ---------------------------------------------------------------------------
# SIR2024_5130 — thin RegressionTable subclass
# ---------------------------------------------------------------------------


class SIR2024_5130(RegressionTable):
    """
    :class:`RegressionTable` pre-configured for SIR 2024-5130 (Tennessee).

    The subclass adds:

    * ``publication`` string set automatically.
    * :meth:`write_table4_template` — writes a blank CSV with all four TN
      areas and the eight standard AEPs, ready to fill from the report PDF.
    * :meth:`load_from_dict` — still accepts the old ``b0``/``b1``/``b2``
      style; the predictor-code mapping is derived from each TN region's
      ``required_predictors`` tuple.

    Usage
    -----
    >>> table = SIR2024_5130.load_from_csv("sir2024_5130_table4.csv")
    >>> q100 = table.estimate(basin, aep=0.01)
    """

    def __init__(
        self,
        region_registry: Optional[RegionRegistry] = None,
    ) -> None:
        registry = region_registry or RegionRegistry()
        for r in TN_REGIONS:
            registry.register(r)
        super().__init__(publication=_PUB, region_registry=registry)

    # ------------------------------------------------------------------
    # Constructors (override to pre-register TN regions)
    # ------------------------------------------------------------------

    @classmethod
    def load_from_dict(  # type: ignore[override]
        cls,
        coeff_dict: Dict[Tuple[HydrologicRegion, float], dict],
        *,
        publication: str = "",
    ) -> "SIR2024_5130":
        """
        Build from a dict of coefficients.

        Accepts both new-style (``"intercept"`` / ``"coefficients"``) and
        legacy-style (``"b0"`` / ``"b1"`` / ``"b2"``) dicts.

        See module docstring for examples.
        """
        table = cls()
        if publication:
            table.publication = publication
        for (region, aep), kw in coeff_dict.items():
            eq = _dict_to_equation(region, aep, kw)
            table._add(eq)
            table._registry.register(region)
        return table

    @classmethod
    def load_from_csv(  # type: ignore[override]
        cls,
        path: str | Path,
        *,
        publication: str = "",
    ) -> "SIR2024_5130":
        """Load equation table from a CSV file (pre-registers TN regions)."""
        table = cls()
        if publication:
            table.publication = publication
        # Delegate to the generic loader, passing our pre-seeded registry
        generic = RegressionTable.load_from_csv(
            path, publication=publication, region_registry=table._registry
        )
        table._eqs = generic._eqs
        if not publication and generic.publication:
            table.publication = generic.publication
        return table

    # ------------------------------------------------------------------
    # Template helper
    # ------------------------------------------------------------------

    @classmethod
    def write_table4_template(cls, path: str | Path) -> None:
        """
        Write a blank Table 4 CSV template for all four TN areas.

        Open SIR 2024-5130 (https://doi.org/10.3133/sir20245130) and fill
        in the coefficients from Table 4 for each row.

        Parameters
        ----------
        path : str or Path
        """
        RegressionTable.write_template(
            path=path,
            regions=list(TN_REGIONS),
            aeps=STANDARD_AEPS,
        )
        log.info("SIR 2024-5130 Table 4 template written to %s", path)

    # ------------------------------------------------------------------
    # Legacy attribute helpers (kept for backward compatibility)
    # ------------------------------------------------------------------

    def get_equation(  # type: ignore[override]
        self, region: HydrologicRegion, aep: float
    ) -> GlsEquation:
        """Retrieve the equation for a given TN area and AEP."""
        return super().get_equation(region, aep)

    def get_variance(self, region: HydrologicRegion, aep: float) -> Optional[float]:
        """Return the prediction variance (log10 space) for an equation."""
        return super().get_variance(region, aep)
