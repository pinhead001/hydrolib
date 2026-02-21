"""
hydrolib.regression.states.wyoming - Wyoming flood-frequency regression equations.

Three hydrologic regions reflecting Wyoming's mountain-to-plains gradient.
Equations follow the log10-linear GLS form::

    log10(Q_T) = intercept + b1*log10(DRNAREA) [+ b2*log10(predictor2) + ...]

Region definitions
------------------
``WY_MOUNTAIN``  — Mountain (Wind River, Tetons, Bighorn, etc.) — DRNAREA + ELEV + PRECIP
``WY_FOOTHILLS`` — Foothills and Intermontane Basins — DRNAREA + CSL1085LFP
``WY_PLAINS``    — Eastern Plains and High Desert — DRNAREA only

Usage
-----
::

    from hydrolib.regression.states.wyoming import (
        WY_MOUNTAIN, WY_FOOTHILLS, WY_PLAINS,
        WY_REGIONS, build_wyoming_table,
    )
    from hydrolib.regression.basin_chars import BasinCharacteristics, DRNAREA, ELEV, PRECIP

    table = build_wyoming_table()
    basin = BasinCharacteristics(
        site_no="06279940",
        site_name="Greybull River above Meeteetse, WY",
        region=WY_MOUNTAIN,
        predictors={DRNAREA: 640.0, ELEV: 7400.0, PRECIP: 22.0},
    )
    q100 = table.estimate(basin, aep=0.01)

Note
----
Coefficient values embedded here are representative for illustration and
testing.  Use authoritative published values for engineering design.
Reference: Lowham, H.W. (1988), USGS WRIR 87-4219; and subsequent updates.
"""

from __future__ import annotations

from typing import Dict, Tuple

from hydrolib.regression.region import HydrologicRegion
from hydrolib.regression.regression_table import GlsEquation, RegressionTable

_PUB = "Lowham (1988), USGS WRIR 87-4219; https://doi.org/10.3133/wri874219"

# ---------------------------------------------------------------------------
# Region definitions
# ---------------------------------------------------------------------------

WY_MOUNTAIN = HydrologicRegion(
    code="WY-1",
    label="Wyoming Mountain",
    state="WY",
    required_predictors=("DRNAREA", "ELEV", "PRECIP"),
    publication=_PUB,
    notes="ELEV = mean basin elevation (ft); PRECIP = mean annual precipitation (in). "
    "Wind River Range, Bighorn Mountains, Tetons, and Absaroka Range; "
    "snowmelt-dominated; perennial streams.",
)

WY_FOOTHILLS = HydrologicRegion(
    code="WY-2",
    label="Wyoming Foothills and Intermontane Basins",
    state="WY",
    required_predictors=("DRNAREA", "CSL1085LFP"),
    publication=_PUB,
    notes="CSL1085LFP = 10-85 channel slope (ft/mi). Intermontane basins and "
    "foothills of central and western Wyoming; mixed runoff regime.",
)

WY_PLAINS = HydrologicRegion(
    code="WY-3",
    label="Wyoming Eastern Plains and High Desert",
    state="WY",
    required_predictors=("DRNAREA",),
    publication=_PUB,
    notes="Drainage area is the sole predictor. Semi-arid High Plains and "
    "Wyoming Basin; sparse, episodic runoff; ephemeral and intermittent streams.",
)

WY_REGIONS: Tuple[HydrologicRegion, ...] = (WY_MOUNTAIN, WY_FOOTHILLS, WY_PLAINS)

# ---------------------------------------------------------------------------
# Regression coefficients (illustrative — use published values for design)
# ---------------------------------------------------------------------------
# Format: {aep: (intercept, {code: exponent}, sep_pct, pseudo_r2, eyr)}

_WY1_COEFFS: Dict[float, tuple] = {
    0.500: (-1.82, {"DRNAREA": 0.76, "ELEV": 0.48, "PRECIP": 0.62}, 52.0, 0.85, 9.0),
    0.200: (-1.50, {"DRNAREA": 0.77, "ELEV": 0.49, "PRECIP": 0.64}, 48.0, 0.86, 10.0),
    0.100: (-1.28, {"DRNAREA": 0.77, "ELEV": 0.50, "PRECIP": 0.65}, 45.0, 0.87, 11.0),
    0.040: (-1.03, {"DRNAREA": 0.78, "ELEV": 0.51, "PRECIP": 0.66}, 43.0, 0.88, 12.0),
    0.020: (-0.87, {"DRNAREA": 0.78, "ELEV": 0.51, "PRECIP": 0.67}, 42.0, 0.88, 12.0),
    0.010: (-0.72, {"DRNAREA": 0.78, "ELEV": 0.52, "PRECIP": 0.67}, 42.0, 0.88, 12.0),
    0.005: (-0.57, {"DRNAREA": 0.79, "ELEV": 0.52, "PRECIP": 0.68}, 44.0, 0.87, 11.0),
    0.002: (-0.37, {"DRNAREA": 0.79, "ELEV": 0.53, "PRECIP": 0.68}, 48.0, 0.86, 10.0),
}

_WY2_COEFFS: Dict[float, tuple] = {
    0.500: (0.38, {"DRNAREA": 0.72, "CSL1085LFP": 0.24}, 62.0, 0.79, 7.0),
    0.200: (0.68, {"DRNAREA": 0.73, "CSL1085LFP": 0.25}, 58.0, 0.81, 8.0),
    0.100: (0.89, {"DRNAREA": 0.73, "CSL1085LFP": 0.26}, 55.0, 0.82, 9.0),
    0.040: (1.12, {"DRNAREA": 0.74, "CSL1085LFP": 0.27}, 53.0, 0.83, 9.0),
    0.020: (1.27, {"DRNAREA": 0.74, "CSL1085LFP": 0.27}, 52.0, 0.84, 10.0),
    0.010: (1.40, {"DRNAREA": 0.75, "CSL1085LFP": 0.28}, 52.0, 0.84, 10.0),
    0.005: (1.53, {"DRNAREA": 0.75, "CSL1085LFP": 0.28}, 54.0, 0.83, 9.0),
    0.002: (1.70, {"DRNAREA": 0.76, "CSL1085LFP": 0.29}, 58.0, 0.82, 8.0),
}

_WY3_COEFFS: Dict[float, tuple] = {
    0.500: (0.52, {"DRNAREA": 0.67}, 76.0, 0.66, 4.0),
    0.200: (0.85, {"DRNAREA": 0.68}, 72.0, 0.68, 5.0),
    0.100: (1.06, {"DRNAREA": 0.69}, 69.0, 0.70, 5.0),
    0.040: (1.31, {"DRNAREA": 0.69}, 67.0, 0.71, 6.0),
    0.020: (1.47, {"DRNAREA": 0.70}, 66.0, 0.72, 6.0),
    0.010: (1.61, {"DRNAREA": 0.70}, 66.0, 0.72, 6.0),
    0.005: (1.75, {"DRNAREA": 0.71}, 68.0, 0.71, 6.0),
    0.002: (1.93, {"DRNAREA": 0.71}, 72.0, 0.69, 5.0),
}

_REGION_COEFFS = {
    WY_MOUNTAIN: _WY1_COEFFS,
    WY_FOOTHILLS: _WY2_COEFFS,
    WY_PLAINS: _WY3_COEFFS,
}


def build_wyoming_table() -> RegressionTable:
    """
    Build a :class:`RegressionTable` for all three Wyoming hydrologic regions.

    Returns a table with 24 equations (3 regions × 8 standard AEPs).

    Returns
    -------
    RegressionTable

    Notes
    -----
    Coefficient values are representative for illustration.  Use authoritative
    published values from the applicable USGS report for engineering design.
    """
    table = RegressionTable(publication=_PUB)
    for region in WY_REGIONS:
        for aep, (intercept, coeff_map, sep_pct, pseudo_r2, eyr) in _REGION_COEFFS[region].items():
            table.add_equation(
                GlsEquation(
                    region=region,
                    aep=aep,
                    intercept=intercept,
                    coefficients=dict(coeff_map),
                    sep_pct=sep_pct,
                    pseudo_r2=pseudo_r2,
                    eyr=eyr,
                )
            )
    return table
