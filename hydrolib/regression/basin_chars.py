"""
hydrolib.regression.basin_chars - Basin characteristics for regional regression.

Stores StreamStats-derived basin attributes needed by SIR 2024-5130
regression equations (drainage area, 10-85 channel slope, hydrologic area,
and auxiliary predictors).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Optional

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Hydrologic area enumeration (SIR 2024-5130)
# ---------------------------------------------------------------------------


class HydrologicArea(Enum):
    """
    Hydrologic area designation from SIR 2024-5130 (Tennessee).

    Attributes
    ----------
    AREA1 : str
        West Tennessee lowlands; predictor variables: DA + 2-yr climate factor.
    AREA2 : str
        Middle Tennessee; predictor variables: DA + 10-85 channel slope (S1085).
    AREA3 : str
        Upper Cumberland / East Tennessee Valley; predictor variable: DA only.
    AREA4 : str
        Ridge and Valley / Blue Ridge; predictor variables: DA + % impervious.
    """

    AREA1 = "area1"
    AREA2 = "area2"
    AREA3 = "area3"
    AREA4 = "area4"

    @classmethod
    def from_int(cls, n: int) -> "HydrologicArea":
        """Construct from integer (1–4)."""
        mapping = {1: cls.AREA1, 2: cls.AREA2, 3: cls.AREA3, 4: cls.AREA4}
        if n not in mapping:
            raise ValueError(f"Hydrologic area must be 1-4, got {n}")
        return mapping[n]

    @property
    def label(self) -> str:
        """Human-readable label."""
        return f"Hydrologic Area {self.value[-1].upper()}"


# ---------------------------------------------------------------------------
# Basin characteristics dataclass
# ---------------------------------------------------------------------------


@dataclass
class BasinCharacteristics:
    """
    Basin characteristics for a site, derived from StreamStats.

    Parameters
    ----------
    site_no : str
        USGS 8-digit site number (or descriptive label for ungaged sites).
    site_name : str
        Descriptive site name.
    drainage_area_sqmi : float
        Drainage area in square miles (from StreamStats watershed delineation).
    hydrologic_area : HydrologicArea
        SIR 2024-5130 hydrologic area assignment.
    slope_1085_ftmi : float, optional
        10-85 channel slope in feet per mile (required for AREA2).
        Computed by StreamStats from 10% and 85% of channel length.
    climate_factor_2yr : float, optional
        2-year recurrence-interval precipitation intensity (in/hr or inches),
        required for AREA1.  Obtain from NOAA Atlas 14 or StreamStats.
    pct_impervious : float, optional
        Percentage of drainage basin covered by impervious surfaces (0-100),
        required for AREA4.  Obtain from NLCD via StreamStats.
    latitude : float, optional
        Outlet latitude (decimal degrees, WGS84).  Used for ROI distance.
    longitude : float, optional
        Outlet longitude (decimal degrees, WGS84).  Used for ROI distance.
    state : str
        State abbreviation (default "TN").
    notes : str
        Free-text field for site-specific notes (karst, regulation, etc.).

    Examples
    --------
    >>> from hydrolib.regression.basin_chars import BasinCharacteristics, HydrologicArea
    >>> site = BasinCharacteristics(
    ...     site_no="03606500",
    ...     site_name="Big Sandy River at Bruceton, TN",
    ...     drainage_area_sqmi=779.0,
    ...     hydrologic_area=HydrologicArea.AREA2,
    ...     slope_1085_ftmi=2.4,
    ... )
    """

    site_no: str
    site_name: str
    drainage_area_sqmi: float
    hydrologic_area: HydrologicArea
    slope_1085_ftmi: Optional[float] = None
    climate_factor_2yr: Optional[float] = None
    pct_impervious: Optional[float] = None
    latitude: Optional[float] = None
    longitude: Optional[float] = None
    state: str = "TN"
    notes: str = ""

    def __post_init__(self) -> None:
        """Validate required predictors for the assigned hydrologic area."""
        if self.drainage_area_sqmi <= 0:
            raise ValueError(f"drainage_area_sqmi must be > 0, got {self.drainage_area_sqmi}")
        if self.hydrologic_area == HydrologicArea.AREA2 and self.slope_1085_ftmi is None:
            raise ValueError("slope_1085_ftmi is required for HydrologicArea.AREA2 (SIR 2024-5130)")
        if self.hydrologic_area == HydrologicArea.AREA1 and self.climate_factor_2yr is None:
            log.warning(
                "climate_factor_2yr not set for AREA1 site %s; GLS computation will fail.",
                self.site_no,
            )
        if self.hydrologic_area == HydrologicArea.AREA4 and self.pct_impervious is None:
            log.warning(
                "pct_impervious not set for AREA4 site %s; GLS computation will fail.",
                self.site_no,
            )

    @classmethod
    def from_streamstats(
        cls,
        site_no: str,
        site_name: str,
        hydrologic_area: HydrologicArea,
        streamstats_params: dict,
    ) -> "BasinCharacteristics":
        """
        Construct from a StreamStats parameter dictionary.

        StreamStats returns basin characteristics as a list of parameter
        objects; this helper extracts the relevant fields by their
        standard StreamStats parameter codes.

        Parameters
        ----------
        site_no : str
            USGS site number or descriptive label.
        site_name : str
            Station name.
        hydrologic_area : HydrologicArea
            SIR 2024-5130 area designation.
        streamstats_params : dict
            Dict mapping StreamStats parameter code -> value.
            Expected keys (all optional except ``DRNAREA``):

            * ``DRNAREA``  – drainage area (sq mi)
            * ``CSL1085LFP`` – 10-85 channel slope (ft/mi)
            * ``I2`` – 2-year precipitation intensity (in/hr)
            * ``IMPERV`` – % impervious
            * ``LAT`` – outlet latitude
            * ``LNG`` – outlet longitude

        Returns
        -------
        BasinCharacteristics
        """
        da = streamstats_params.get("DRNAREA")
        if da is None:
            raise KeyError("'DRNAREA' not found in streamstats_params")

        return cls(
            site_no=site_no,
            site_name=site_name,
            drainage_area_sqmi=float(da),
            hydrologic_area=hydrologic_area,
            slope_1085_ftmi=streamstats_params.get("CSL1085LFP"),
            climate_factor_2yr=streamstats_params.get("I2"),
            pct_impervious=streamstats_params.get("IMPERV"),
            latitude=streamstats_params.get("LAT"),
            longitude=streamstats_params.get("LNG"),
        )

    def summary(self) -> str:
        """Return a formatted one-paragraph summary."""
        lines = [
            f"Site: {self.site_no} — {self.site_name}",
            f"  Hydrologic area : {self.hydrologic_area.label}",
            f"  Drainage area   : {self.drainage_area_sqmi:.2f} sq mi",
        ]
        if self.slope_1085_ftmi is not None:
            lines.append(f"  S1085 slope     : {self.slope_1085_ftmi:.3f} ft/mi")
        if self.climate_factor_2yr is not None:
            lines.append(f"  2-yr CF         : {self.climate_factor_2yr:.3f}")
        if self.pct_impervious is not None:
            lines.append(f"  % impervious    : {self.pct_impervious:.1f}%")
        if self.notes:
            lines.append(f"  Notes           : {self.notes}")
        return "\n".join(lines)
