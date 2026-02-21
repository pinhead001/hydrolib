"""
hydrolib.regression.states.montana - Montana flood-frequency regression equations.

Three hydrologic regions based on USGS regression studies for Montana streams
(Parrett and Johnson, 2004, USGS WRIR 03-4377, and subsequent updates).
Equations follow the log10-linear form::

    log10(Q_T) = intercept + b1*log10(DRNAREA) [+ b2*log10(predictor2) + ...]

Region definitions
------------------
``MT_MOUNTAIN``  — Mountain (western/northwestern MT) — DRNAREA + ELEV + PRECIP
``MT_FOOTHILLS`` — Foothills (central MT) — DRNAREA + CSL1085LFP
``MT_PLAINS``    — Eastern Plains — DRNAREA only

Usage
-----
::

    from hydrolib.regression.states.montana import (
        MT_MOUNTAIN, MT_FOOTHILLS, MT_PLAINS,
        MT_REGIONS, build_montana_table,
    )
    from hydrolib.regression.basin_chars import BasinCharacteristics, DRNAREA, ELEV, PRECIP

    table = build_montana_table()
    basin = BasinCharacteristics(
        site_no="MT-UNGAGED",
        site_name="Lost Horse Creek near Hamilton, MT",
        region=MT_MOUNTAIN,
        predictors={DRNAREA: 45.0, ELEV: 5800.0, PRECIP: 32.0},
    )
    q100 = table.estimate(basin, aep=0.01)

Note
----
Coefficient values embedded here are representative for illustration and
testing.  Use authoritative published values for engineering design.
"""

from __future__ import annotations

from typing import Dict, Tuple

from hydrolib.regression.region import HydrologicRegion
from hydrolib.regression.regression_table import STANDARD_AEPS, GlsEquation, RegressionTable

_PUB = "Parrett and Johnson (2004), USGS WRIR 03-4377; " "https://doi.org/10.3133/wri034377"

# ---------------------------------------------------------------------------
# Region definitions
# ---------------------------------------------------------------------------

MT_MOUNTAIN = HydrologicRegion(
    code="MT-1",
    label="Montana Mountain",
    state="MT",
    required_predictors=("DRNAREA", "ELEV", "PRECIP"),
    publication=_PUB,
    notes="ELEV = mean basin elevation (ft); PRECIP = mean annual precipitation (in). "
    "Western and northwestern Montana; snow-dominated hydrology; "
    "primarily Glacier, Flathead, Missoula, and adjacent counties.",
)

MT_FOOTHILLS = HydrologicRegion(
    code="MT-2",
    label="Montana Foothills",
    state="MT",
    required_predictors=("DRNAREA", "CSL1085LFP"),
    publication=_PUB,
    notes="CSL1085LFP = 10-85 channel slope (ft/mi). "
    "Central Montana foothills east of the Continental Divide; "
    "mixed rainfall-snowmelt runoff regime.",
)

MT_PLAINS = HydrologicRegion(
    code="MT-3",
    label="Montana Eastern Plains",
    state="MT",
    required_predictors=("DRNAREA",),
    publication=_PUB,
    notes="Drainage area is the sole predictor. "
    "Eastern Montana plains; semi-arid; high variability in annual peak flows.",
)

MT_REGIONS: Tuple[HydrologicRegion, ...] = (MT_MOUNTAIN, MT_FOOTHILLS, MT_PLAINS)

# ---------------------------------------------------------------------------
# Regression coefficients
# ---------------------------------------------------------------------------
# Format: {aep: (intercept, {code: exponent}, sep_pct, pseudo_r2, eyr)}

_MT1_COEFFS: Dict[float, tuple] = {
    0.500: (-1.68, {"DRNAREA": 0.74, "ELEV": 0.45, "PRECIP": 0.58}, 54.0, 0.84, 9.0),
    0.200: (-1.35, {"DRNAREA": 0.75, "ELEV": 0.46, "PRECIP": 0.60}, 50.0, 0.85, 10.0),
    0.100: (-1.14, {"DRNAREA": 0.75, "ELEV": 0.47, "PRECIP": 0.61}, 47.0, 0.86, 11.0),
    0.040: (-0.91, {"DRNAREA": 0.76, "ELEV": 0.47, "PRECIP": 0.62}, 45.0, 0.87, 12.0),
    0.020: (-0.75, {"DRNAREA": 0.76, "ELEV": 0.48, "PRECIP": 0.62}, 44.0, 0.87, 12.0),
    0.010: (-0.61, {"DRNAREA": 0.76, "ELEV": 0.48, "PRECIP": 0.62}, 44.0, 0.87, 12.0),
    0.005: (-0.47, {"DRNAREA": 0.77, "ELEV": 0.48, "PRECIP": 0.63}, 46.0, 0.86, 11.0),
    0.002: (-0.27, {"DRNAREA": 0.77, "ELEV": 0.49, "PRECIP": 0.63}, 50.0, 0.85, 10.0),
}

_MT2_COEFFS: Dict[float, tuple] = {
    0.500: (0.45, {"DRNAREA": 0.73, "CSL1085LFP": 0.25}, 58.0, 0.80, 8.0),
    0.200: (0.70, {"DRNAREA": 0.73, "CSL1085LFP": 0.26}, 54.0, 0.82, 9.0),
    0.100: (0.86, {"DRNAREA": 0.74, "CSL1085LFP": 0.27}, 51.0, 0.83, 10.0),
    0.040: (1.05, {"DRNAREA": 0.74, "CSL1085LFP": 0.27}, 49.0, 0.84, 11.0),
    0.020: (1.18, {"DRNAREA": 0.75, "CSL1085LFP": 0.28}, 48.0, 0.84, 11.0),
    0.010: (1.30, {"DRNAREA": 0.75, "CSL1085LFP": 0.28}, 48.0, 0.84, 11.0),
    0.005: (1.42, {"DRNAREA": 0.75, "CSL1085LFP": 0.28}, 50.0, 0.83, 10.0),
    0.002: (1.58, {"DRNAREA": 0.76, "CSL1085LFP": 0.29}, 54.0, 0.82, 9.0),
}

_MT3_COEFFS: Dict[float, tuple] = {
    0.500: (0.48, {"DRNAREA": 0.68}, 72.0, 0.70, 5.0),
    0.200: (0.78, {"DRNAREA": 0.69}, 68.0, 0.72, 6.0),
    0.100: (0.98, {"DRNAREA": 0.70}, 65.0, 0.73, 6.0),
    0.040: (1.20, {"DRNAREA": 0.70}, 63.0, 0.74, 7.0),
    0.020: (1.35, {"DRNAREA": 0.71}, 62.0, 0.74, 7.0),
    0.010: (1.48, {"DRNAREA": 0.71}, 62.0, 0.74, 7.0),
    0.005: (1.61, {"DRNAREA": 0.72}, 64.0, 0.73, 6.0),
    0.002: (1.79, {"DRNAREA": 0.72}, 68.0, 0.72, 6.0),
}

_REGION_COEFFS = {
    MT_MOUNTAIN: _MT1_COEFFS,
    MT_FOOTHILLS: _MT2_COEFFS,
    MT_PLAINS: _MT3_COEFFS,
}


def build_montana_table() -> RegressionTable:
    """
    Build a :class:`RegressionTable` for all three Montana hydrologic regions.

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
    for region in MT_REGIONS:
        coeffs = _REGION_COEFFS[region]
        for aep, (intercept, coeff_map, sep_pct, pseudo_r2, eyr) in coeffs.items():
            eq = GlsEquation(
                region=region,
                aep=aep,
                intercept=intercept,
                coefficients=dict(coeff_map),
                sep_pct=sep_pct,
                pseudo_r2=pseudo_r2,
                eyr=eyr,
            )
            table.add_equation(eq)
    return table
