"""
hydrolib.regression.basin_chars - Generic basin characteristics container.

:class:`BasinCharacteristics` stores the hydraulic and physical attributes
of a drainage basin required by regional regression equations.  All predictor
values are stored in a ``predictors`` dict keyed by **USGS StreamStats
parameter codes** so the same class works for any state and any predictor
combination.

Common StreamStats codes (use the constants below to avoid magic strings)
--------------------------------------------------------------------------
``DRNAREA``    — Drainage area (sq mi)
``CSL1085LFP`` — 10-85 channel slope (ft/mi)
``ELEV``       — Mean basin elevation (ft)
``PRECIP``     — Mean annual precipitation (in)
``I2``         — 2-year precipitation intensity (in/hr or in)
``IMPERV``     — Impervious cover (%)
``FOREST``     — Forest cover (%)
``BFIPCT``     — Base-flow index (%)
``SLOPE``      — Mean basin slope (%)
``LC06IMP``    — 2006 NLCD impervious cover (%)

Full list: https://streamstats.usgs.gov/ss/
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, Optional

from hydrolib.regression.region import HydrologicRegion

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# StreamStats code constants (use these to avoid magic strings in call sites)
# ---------------------------------------------------------------------------

DRNAREA = "DRNAREA"  # Drainage area (sq mi) — universal
CSL1085LFP = "CSL1085LFP"  # 10-85 channel slope (ft/mi)
ELEV = "ELEV"  # Mean basin elevation (ft)
PRECIP = "PRECIP"  # Mean annual precipitation (in)
I2 = "I2"  # 2-yr precipitation intensity (in/hr or in)
IMPERV = "IMPERV"  # Impervious cover (%)
FOREST = "FOREST"  # Forest cover (%)
BFIPCT = "BFIPCT"  # Base-flow index (%)
SLOPE = "SLOPE"  # Mean basin slope (%)
LC06IMP = "LC06IMP"  # 2006 NLCD impervious cover (%)


# ---------------------------------------------------------------------------
# BasinCharacteristics
# ---------------------------------------------------------------------------


@dataclass
class BasinCharacteristics:
    """
    Physical and climatic characteristics of a drainage basin.

    All predictor values live in ``predictors``, a dict keyed by USGS
    StreamStats parameter codes.  This design is state-agnostic: Tennessee
    equations need DRNAREA and CSL1085LFP; Montana mountain equations need
    DRNAREA, ELEV, and PRECIP; both use the same class.

    Parameters
    ----------
    site_no : str
        USGS 8-digit site number, or a descriptive label for ungaged sites.
    site_name : str
        Station or site name.
    region : HydrologicRegion
        Hydrologic region to which this site is assigned.
    predictors : dict
        Mapping of StreamStats parameter code → value.  Must include all
        codes listed in ``region.required_predictors`` with positive values.
    latitude : float, optional
        Outlet latitude (decimal degrees, WGS84).
    longitude : float, optional
        Outlet longitude (decimal degrees, WGS84).
    notes : str, optional
        Free-text field (karst observations, regulation flags, etc.).

    Examples
    --------
    Tennessee Middle TN (Area 2) — DA + channel slope::

        from hydrolib.regression.states.tennessee import TN_AREA2
        from hydrolib.regression.basin_chars import DRNAREA, CSL1085LFP

        basin = BasinCharacteristics(
            site_no="03606500",
            site_name="Big Sandy River at Bruceton, TN",
            region=TN_AREA2,
            predictors={DRNAREA: 779.0, CSL1085LFP: 2.4},
        )

    Montana mountain — DA + elevation + precipitation::

        from hydrolib.regression.states.montana import MT_MOUNTAIN
        from hydrolib.regression.basin_chars import DRNAREA, ELEV, PRECIP

        mt_basin = BasinCharacteristics(
            site_no="MT-UNGAGED",
            site_name="Lost Horse Creek near Hamilton, MT",
            region=MT_MOUNTAIN,
            predictors={DRNAREA: 45.0, ELEV: 5800.0, PRECIP: 32.0},
        )

    Georgia Blue Ridge — DA + elevation::

        from hydrolib.regression.states.georgia import GA_BLUE_RIDGE
        from hydrolib.regression.basin_chars import DRNAREA, ELEV

        ga_basin = BasinCharacteristics(
            site_no="02333500",
            site_name="Chattahoochee River near Helen, GA",
            region=GA_BLUE_RIDGE,
            predictors={DRNAREA: 47.0, ELEV: 2650.0},
        )
    """

    site_no: str
    site_name: str
    region: HydrologicRegion
    predictors: Dict[str, float] = field(default_factory=dict)
    latitude: Optional[float] = None
    longitude: Optional[float] = None
    notes: str = ""

    def __post_init__(self) -> None:
        """Validate that all required predictors are present and positive."""
        try:
            self.region.validate_predictors(self.predictors)
        except ValueError as exc:
            raise ValueError(
                f"BasinCharacteristics validation failed for site {self.site_no!r}: {exc}"
            ) from exc

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def state(self) -> str:
        """State abbreviation, inherited from the assigned region."""
        return self.region.state

    # ------------------------------------------------------------------
    # Factory helpers
    # ------------------------------------------------------------------

    @classmethod
    def from_streamstats(
        cls,
        site_no: str,
        site_name: str,
        region: HydrologicRegion,
        streamstats_params: Dict[str, float],
        *,
        latitude: Optional[float] = None,
        longitude: Optional[float] = None,
        notes: str = "",
    ) -> "BasinCharacteristics":
        """
        Construct directly from a StreamStats parameter dictionary.

        All numeric entries from *streamstats_params* are forwarded to
        ``predictors``; non-numeric values are silently dropped.  Required
        predictor validation still runs via ``__post_init__``.

        Parameters
        ----------
        site_no, site_name : str
        region : HydrologicRegion
        streamstats_params : dict
            Full StreamStats response mapping code → value, e.g.::

                {"DRNAREA": 85.3, "CSL1085LFP": 4.7, "ELEV": 1240.0}

        latitude, longitude : float, optional
        notes : str, optional

        Returns
        -------
        BasinCharacteristics
        """
        predictors: Dict[str, float] = {}
        for k, v in streamstats_params.items():
            try:
                predictors[k] = float(v)
            except (TypeError, ValueError):
                log.debug("Skipping non-numeric StreamStats param %s=%r", k, v)

        return cls(
            site_no=site_no,
            site_name=site_name,
            region=region,
            predictors=predictors,
            latitude=latitude,
            longitude=longitude,
            notes=notes,
        )

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------

    def predictor_value(self, code: str) -> float:
        """
        Return the value of a predictor by StreamStats code.

        Parameters
        ----------
        code : str
            StreamStats parameter code.

        Returns
        -------
        float

        Raises
        ------
        KeyError
            If *code* is not in ``predictors``.
        """
        if code not in self.predictors:
            raise KeyError(
                f"Predictor {code!r} not found for site {self.site_no!r}. "
                f"Available: {sorted(self.predictors)}"
            )
        return self.predictors[code]

    def summary(self) -> str:
        """Return a formatted one-paragraph summary."""
        lines = [
            f"Site  : {self.site_no} — {self.site_name}",
            f"State : {self.region.state}",
            f"Region: {self.region}",
        ]
        for code in sorted(self.predictors):
            lines.append(f"  {code:<16} : {self.predictors[code]:.4g}")
        if self.latitude is not None:
            lines.append(f"  {'Lat/Lon':<16} : {self.latitude:.4f}, {self.longitude:.4f}")
        if self.notes:
            lines.append(f"  {'Notes':<16} : {self.notes}")
        return "\n".join(lines)
