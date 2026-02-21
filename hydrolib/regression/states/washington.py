"""
hydrolib.regression.states.washington - Washington flood-frequency regression equations.

Four hydrologic regions reflecting Washington's strong orographic gradient.
Equations follow the log10-linear GLS form::

    log10(Q_T) = intercept + b1*log10(DRNAREA) [+ b2*log10(predictor2) + ...]

Region definitions
------------------
``WA_OLYMPIC``   — Olympic Peninsula and Coast Ranges — DRNAREA + PRECIP
``WA_W_CASCADE`` — Western Cascades — DRNAREA + ELEV + PRECIP
``WA_E_CASCADE`` — Eastern Cascades — DRNAREA + ELEV
``WA_EASTERN``   — Eastern Washington (Columbia Plateau and Basin) — DRNAREA only

Usage
-----
::

    from hydrolib.regression.states.washington import (
        WA_OLYMPIC, WA_W_CASCADE, WA_E_CASCADE, WA_EASTERN,
        WA_REGIONS, build_washington_table,
    )
    from hydrolib.regression.basin_chars import BasinCharacteristics, DRNAREA, PRECIP

    table = build_washington_table()
    basin = BasinCharacteristics(
        site_no="12039300",
        site_name="Quinault River at Quinault, WA",
        region=WA_OLYMPIC,
        predictors={DRNAREA: 286.0, PRECIP: 130.0},
    )
    q100 = table.estimate(basin, aep=0.01)

Note
----
Coefficient values embedded here are representative for illustration and
testing.  Use authoritative published values for engineering design.
Reference: Curran, C.A., and others (2019), USGS SIR 2019-5073.
"""

from __future__ import annotations

from typing import Dict, Tuple

from hydrolib.regression.region import HydrologicRegion
from hydrolib.regression.regression_table import GlsEquation, RegressionTable

_PUB = "Curran and others (2019), USGS SIR 2019-5073; https://doi.org/10.3133/sir20195073"

# ---------------------------------------------------------------------------
# Region definitions
# ---------------------------------------------------------------------------

WA_OLYMPIC = HydrologicRegion(
    code="WA-1",
    label="Washington Olympic Peninsula and Coast",
    state="WA",
    required_predictors=("DRNAREA", "PRECIP"),
    publication=_PUB,
    notes="PRECIP = mean annual precipitation (in). Olympic Peninsula and "
    "Coast Ranges; maritime climate; highest precipitation in CONUS; "
    "rainfall-dominated; perennial streams.",
)

WA_W_CASCADE = HydrologicRegion(
    code="WA-2",
    label="Washington Western Cascades",
    state="WA",
    required_predictors=("DRNAREA", "ELEV", "PRECIP"),
    publication=_PUB,
    notes="ELEV = mean basin elevation (ft); PRECIP = mean annual precipitation (in). "
    "Western slope of Cascade Range; mixed rain-snow; "
    "includes rain-on-snow flood-generating events.",
)

WA_E_CASCADE = HydrologicRegion(
    code="WA-3",
    label="Washington Eastern Cascades",
    state="WA",
    required_predictors=("DRNAREA", "ELEV"),
    publication=_PUB,
    notes="ELEV = mean basin elevation (ft). Eastern slope of the Cascades "
    "in transition from maritime to continental climate; snowmelt-dominated.",
)

WA_EASTERN = HydrologicRegion(
    code="WA-4",
    label="Washington Eastern (Columbia Plateau and Basin)",
    state="WA",
    required_predictors=("DRNAREA",),
    publication=_PUB,
    notes="Drainage area is the sole predictor. Columbia Plateau and "
    "Okanagan Highlands; semi-arid; intermittent streams; "
    "episodic runoff from spring snowmelt and convective storms.",
)

WA_REGIONS: Tuple[HydrologicRegion, ...] = (WA_OLYMPIC, WA_W_CASCADE, WA_E_CASCADE, WA_EASTERN)

# ---------------------------------------------------------------------------
# Regression coefficients (illustrative — use published values for design)
# ---------------------------------------------------------------------------
# Format: {aep: (intercept, {code: exponent}, sep_pct, pseudo_r2, eyr)}

_WA1_COEFFS: Dict[float, tuple] = {
    0.500: (0.32, {"DRNAREA": 0.82, "PRECIP": 0.60}, 42.0, 0.91, 15.0),
    0.200: (0.62, {"DRNAREA": 0.83, "PRECIP": 0.62}, 38.0, 0.92, 17.0),
    0.100: (0.82, {"DRNAREA": 0.83, "PRECIP": 0.63}, 36.0, 0.93, 18.0),
    0.040: (1.04, {"DRNAREA": 0.84, "PRECIP": 0.64}, 35.0, 0.93, 19.0),
    0.020: (1.18, {"DRNAREA": 0.84, "PRECIP": 0.65}, 35.0, 0.93, 19.0),
    0.010: (1.31, {"DRNAREA": 0.85, "PRECIP": 0.65}, 35.0, 0.93, 19.0),
    0.005: (1.44, {"DRNAREA": 0.85, "PRECIP": 0.66}, 37.0, 0.93, 18.0),
    0.002: (1.61, {"DRNAREA": 0.86, "PRECIP": 0.66}, 40.0, 0.92, 17.0),
}

_WA2_COEFFS: Dict[float, tuple] = {
    0.500: (-0.55, {"DRNAREA": 0.80, "ELEV": 0.42, "PRECIP": 0.55}, 45.0, 0.90, 14.0),
    0.200: (-0.24, {"DRNAREA": 0.81, "ELEV": 0.43, "PRECIP": 0.57}, 41.0, 0.91, 15.0),
    0.100: (-0.03, {"DRNAREA": 0.81, "ELEV": 0.44, "PRECIP": 0.58}, 39.0, 0.92, 16.0),
    0.040: (0.21, {"DRNAREA": 0.82, "ELEV": 0.45, "PRECIP": 0.59}, 37.0, 0.92, 17.0),
    0.020: (0.37, {"DRNAREA": 0.82, "ELEV": 0.45, "PRECIP": 0.60}, 37.0, 0.92, 17.0),
    0.010: (0.51, {"DRNAREA": 0.83, "ELEV": 0.46, "PRECIP": 0.60}, 37.0, 0.92, 17.0),
    0.005: (0.65, {"DRNAREA": 0.83, "ELEV": 0.46, "PRECIP": 0.61}, 39.0, 0.92, 16.0),
    0.002: (0.83, {"DRNAREA": 0.84, "ELEV": 0.47, "PRECIP": 0.61}, 43.0, 0.91, 15.0),
}

_WA3_COEFFS: Dict[float, tuple] = {
    0.500: (-0.88, {"DRNAREA": 0.81, "ELEV": 0.55}, 50.0, 0.86, 10.0),
    0.200: (-0.57, {"DRNAREA": 0.82, "ELEV": 0.56}, 46.0, 0.87, 11.0),
    0.100: (-0.36, {"DRNAREA": 0.82, "ELEV": 0.57}, 43.0, 0.88, 12.0),
    0.040: (-0.11, {"DRNAREA": 0.83, "ELEV": 0.58}, 41.0, 0.89, 13.0),
    0.020: (0.05, {"DRNAREA": 0.83, "ELEV": 0.58}, 40.0, 0.90, 14.0),
    0.010: (0.19, {"DRNAREA": 0.84, "ELEV": 0.59}, 40.0, 0.90, 14.0),
    0.005: (0.33, {"DRNAREA": 0.84, "ELEV": 0.59}, 42.0, 0.89, 13.0),
    0.002: (0.51, {"DRNAREA": 0.85, "ELEV": 0.60}, 46.0, 0.88, 12.0),
}

_WA4_COEFFS: Dict[float, tuple] = {
    0.500: (1.00, {"DRNAREA": 0.70}, 68.0, 0.72, 6.0),
    0.200: (1.32, {"DRNAREA": 0.71}, 64.0, 0.74, 7.0),
    0.100: (1.53, {"DRNAREA": 0.72}, 62.0, 0.75, 7.0),
    0.040: (1.77, {"DRNAREA": 0.72}, 60.0, 0.76, 8.0),
    0.020: (1.92, {"DRNAREA": 0.73}, 59.0, 0.77, 8.0),
    0.010: (2.06, {"DRNAREA": 0.73}, 59.0, 0.77, 8.0),
    0.005: (2.19, {"DRNAREA": 0.74}, 61.0, 0.76, 8.0),
    0.002: (2.37, {"DRNAREA": 0.74}, 64.0, 0.75, 7.0),
}

_REGION_COEFFS = {
    WA_OLYMPIC: _WA1_COEFFS,
    WA_W_CASCADE: _WA2_COEFFS,
    WA_E_CASCADE: _WA3_COEFFS,
    WA_EASTERN: _WA4_COEFFS,
}


def build_washington_table() -> RegressionTable:
    """
    Build a :class:`RegressionTable` for all four Washington hydrologic regions.

    Returns a table with 32 equations (4 regions × 8 standard AEPs).

    Returns
    -------
    RegressionTable

    Notes
    -----
    Coefficient values are representative for illustration.  Use authoritative
    published values from the applicable USGS report for engineering design.
    """
    table = RegressionTable(publication=_PUB)
    for region in WA_REGIONS:
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
