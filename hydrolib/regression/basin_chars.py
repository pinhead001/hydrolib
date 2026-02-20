"""
hydrolib.regression.basin_chars - Generic basin characteristics container.

:class:`BasinCharacteristics` stores the hydraulic and physical attributes
of a drainage basin required by regional regression equations.  Predictor
values are stored in a flexible ``predictors`` dict keyed by **USGS
StreamStats parameter codes** (e.g. ``"DRNAREA"``, ``"CSL1085LFP"``), so
the same class works for equations in any state or publication.

Common StreamStats codes
------------------------
``DRNAREA``   — Drainage area (sq mi)
``CSL1085LFP``— 10-85 channel slope (ft/mi)
``I2``        — 2-year precipitation intensity (in/hr)
``IMPERV``    — Impervious cover (%)
``ELEV``      — Mean basin elevation (ft)
``FOREST``    — Forest cover (%)
``PRECIP``    — Mean annual precipitation (in)
``BFIPCT``    — Base-flow index (%)

Convenience properties expose the most common codes as named attributes
so that existing code using ``basin.drainage_area_sqmi`` continues to work.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, Optional

from hydrolib.regression.region import HydrologicRegion

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# BasinCharacteristics
# ---------------------------------------------------------------------------

# StreamStats code constants (avoid magic strings in calling code)
DRNAREA = "DRNAREA"
CSL1085LFP = "CSL1085LFP"
I2 = "I2"
IMPERV = "IMPERV"
ELEV = "ELEV"
FOREST = "FOREST"
PRECIP = "PRECIP"
BFIPCT = "BFIPCT"


@dataclass
class BasinCharacteristics:
    """
    Physical and climatic characteristics of a drainage basin.

    Parameters
    ----------
    site_no : str
        USGS 8-digit site number, or a descriptive label for ungaged sites.
    site_name : str
        Station or site name.
    region : HydrologicRegion
        Hydrologic region to which this site is assigned.
    predictors : dict
        Mapping of StreamStats parameter code -> value.  Must include all
        codes listed in ``region.required_predictors``.
    latitude : float, optional
        Outlet latitude (decimal degrees, WGS84).  Used for ROI geographic
        distance.
    longitude : float, optional
        Outlet longitude (decimal degrees, WGS84).
    notes : str, optional
        Free-text field (karst observations, regulation flags, etc.).

    Examples
    --------
    >>> from hydrolib.regression.region import HydrologicRegion
    >>> tn_area2 = HydrologicRegion(
    ...     code="TN_AREA2", label="TN Area 2", state="TN",
    ...     required_predictors=("DRNAREA", "CSL1085LFP"),
    ... )
    >>> basin = BasinCharacteristics(
    ...     site_no="03606500",
    ...     site_name="Big Sandy River at Bruceton, TN",
    ...     region=tn_area2,
    ...     predictors={"DRNAREA": 779.0, "CSL1085LFP": 2.4},
    ... )
    >>> basin.drainage_area_sqmi
    779.0

    Multi-state / arbitrary-predictor example (Kentucky)::

        ky_reg4 = HydrologicRegion(
            code="KY_REGION4", label="KY Region 4", state="KY",
            required_predictors=("DRNAREA", "CSL1085LFP"),
            publication="WRIR 03-4180",
        )
        ky_basin = BasinCharacteristics(
            site_no="03287500",
            site_name="Licking River at Farmers, KY",
            region=ky_reg4,
            predictors={"DRNAREA": 1140.0, "CSL1085LFP": 1.8},
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
    # Convenience properties (most-common StreamStats codes)
    # ------------------------------------------------------------------

    @property
    def drainage_area_sqmi(self) -> Optional[float]:
        """Drainage area in square miles (``DRNAREA``)."""
        return self.predictors.get(DRNAREA)

    @property
    def slope_1085_ftmi(self) -> Optional[float]:
        """10-85 channel slope in ft/mi (``CSL1085LFP``)."""
        return self.predictors.get(CSL1085LFP)

    @property
    def climate_factor_2yr(self) -> Optional[float]:
        """2-year precipitation intensity (``I2``)."""
        return self.predictors.get(I2)

    @property
    def pct_impervious(self) -> Optional[float]:
        """Percent impervious cover (``IMPERV``)."""
        return self.predictors.get(IMPERV)

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

        All keys from *streamstats_params* are forwarded to ``predictors``,
        so the full StreamStats response can be passed without pre-filtering.
        Required predictor validation still runs via ``__post_init__``.

        Parameters
        ----------
        site_no : str
        site_name : str
        region : HydrologicRegion
        streamstats_params : dict
            Mapping StreamStats parameter code -> numeric value.
            Typically obtained from the StreamStats REST API response::

                {
                    "DRNAREA": 85.3,
                    "CSL1085LFP": 4.7,
                    "ELEV": 1240.0,
                    ...
                }

        latitude, longitude : float, optional
        notes : str, optional

        Returns
        -------
        BasinCharacteristics

        Raises
        ------
        ValueError
            If any required predictor for *region* is missing.
        """
        # Convert all values to float, silently skip non-numeric entries
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
            f"Region: {self.region}",
        ]
        for code in sorted(self.predictors):
            lines.append(f"  {code:<16} : {self.predictors[code]:.4g}")
        if self.latitude is not None:
            lines.append(f"  {'Lat/Lon':<16} : {self.latitude:.4f}, {self.longitude:.4f}")
        if self.notes:
            lines.append(f"  {'Notes':<16} : {self.notes}")
        return "\n".join(lines)
