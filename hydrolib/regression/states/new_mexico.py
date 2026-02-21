"""
hydrolib.regression.states.new_mexico - New Mexico flood-frequency regression equations.

Four hydrologic regions based on physiographic and hydrologic similarity in
New Mexico.  Equations follow the log10-linear GLS form::

    log10(Q_T) = intercept + b1*log10(DRNAREA) [+ b2*log10(predictor2) + ...]

Region definitions
------------------
``NM_MOUNTAIN``   — Northern/Northwestern Mountains — DRNAREA + PRECIP
``NM_TRANSITION`` — Transitional Highlands — DRNAREA + CSL1085LFP
``NM_PLAINS``     — Eastern Plains (Great Plains) — DRNAREA only
``NM_DESERT``     — Desert Basin and Range (southwest) — DRNAREA + CSL1085LFP

Usage
-----
::

    from hydrolib.regression.states.new_mexico import (
        NM_MOUNTAIN, NM_TRANSITION, NM_PLAINS, NM_DESERT,
        NM_REGIONS, build_new_mexico_table,
    )
    from hydrolib.regression.basin_chars import BasinCharacteristics, DRNAREA, PRECIP

    table = build_new_mexico_table()
    basin = BasinCharacteristics(
        site_no="08313000",
        site_name="Rio Grande at Otowi Bridge, NM",
        region=NM_MOUNTAIN,
        predictors={DRNAREA: 14357.0, PRECIP: 18.0},
    )
    q100 = table.estimate(basin, aep=0.01)

Note
----
Coefficient values embedded here are representative for illustration and
testing.  Use authoritative published values for engineering design.
Reference: Thomas, B.E., Hjalmarson, H.W., and Waltemeyer, S.D. (1997),
USGS Water-Supply Paper 2386; and subsequent StreamStats updates.
"""

from __future__ import annotations

from typing import Dict, Tuple

from hydrolib.regression.region import HydrologicRegion
from hydrolib.regression.regression_table import GlsEquation, RegressionTable

_PUB = "USGS WSP 2386 / NM StreamStats; https://streamstats.usgs.gov/ss/"

# ---------------------------------------------------------------------------
# Region definitions
# ---------------------------------------------------------------------------

NM_MOUNTAIN = HydrologicRegion(
    code="NM-1",
    label="New Mexico Northern Mountain",
    state="NM",
    required_predictors=("DRNAREA", "PRECIP"),
    publication=_PUB,
    notes="PRECIP = mean annual precipitation (in). Sangre de Cristo, Jemez, "
    "and Sacramento Mountains; snow-dominated to mixed rainfall-snowmelt.",
)

NM_TRANSITION = HydrologicRegion(
    code="NM-2",
    label="New Mexico Transitional Highlands",
    state="NM",
    required_predictors=("DRNAREA", "CSL1085LFP"),
    publication=_PUB,
    notes="CSL1085LFP = 10-85 channel slope (ft/mi). Transition zone between "
    "mountains and plains; mixed runoff regime including convective storms.",
)

NM_PLAINS = HydrologicRegion(
    code="NM-3",
    label="New Mexico Eastern Plains",
    state="NM",
    required_predictors=("DRNAREA",),
    publication=_PUB,
    notes="Drainage area is the sole predictor. Semi-arid Great Plains "
    "in eastern New Mexico; highly variable, convection-dominated runoff.",
)

NM_DESERT = HydrologicRegion(
    code="NM-4",
    label="New Mexico Desert Basin and Range",
    state="NM",
    required_predictors=("DRNAREA", "CSL1085LFP"),
    publication=_PUB,
    notes="CSL1085LFP = 10-85 channel slope (ft/mi). Arid basin-and-range "
    "terrain of southwest NM; sparse runoff with ephemeral streams.",
)

NM_REGIONS: Tuple[HydrologicRegion, ...] = (NM_MOUNTAIN, NM_TRANSITION, NM_PLAINS, NM_DESERT)

# ---------------------------------------------------------------------------
# Regression coefficients (illustrative — use published values for design)
# ---------------------------------------------------------------------------
# Format: {aep: (intercept, {code: exponent}, sep_pct, pseudo_r2, eyr)}

_NM1_COEFFS: Dict[float, tuple] = {
    0.500: (-0.45, {"DRNAREA": 0.72, "PRECIP": 0.65}, 58.0, 0.82, 8.0),
    0.200: (-0.12, {"DRNAREA": 0.73, "PRECIP": 0.68}, 53.0, 0.84, 9.0),
    0.100: (0.10, {"DRNAREA": 0.74, "PRECIP": 0.70}, 50.0, 0.85, 10.0),
    0.040: (0.34, {"DRNAREA": 0.74, "PRECIP": 0.72}, 48.0, 0.86, 11.0),
    0.020: (0.50, {"DRNAREA": 0.75, "PRECIP": 0.73}, 47.0, 0.86, 11.0),
    0.010: (0.64, {"DRNAREA": 0.75, "PRECIP": 0.74}, 47.0, 0.86, 11.0),
    0.005: (0.77, {"DRNAREA": 0.76, "PRECIP": 0.74}, 49.0, 0.85, 10.0),
    0.002: (0.95, {"DRNAREA": 0.76, "PRECIP": 0.75}, 53.0, 0.84, 9.0),
}

_NM2_COEFFS: Dict[float, tuple] = {
    0.500: (0.35, {"DRNAREA": 0.70, "CSL1085LFP": 0.20}, 65.0, 0.77, 6.0),
    0.200: (0.65, {"DRNAREA": 0.70, "CSL1085LFP": 0.21}, 61.0, 0.79, 7.0),
    0.100: (0.85, {"DRNAREA": 0.71, "CSL1085LFP": 0.22}, 58.0, 0.80, 8.0),
    0.040: (1.07, {"DRNAREA": 0.71, "CSL1085LFP": 0.23}, 56.0, 0.81, 8.0),
    0.020: (1.22, {"DRNAREA": 0.72, "CSL1085LFP": 0.23}, 55.0, 0.82, 9.0),
    0.010: (1.35, {"DRNAREA": 0.72, "CSL1085LFP": 0.24}, 55.0, 0.82, 9.0),
    0.005: (1.48, {"DRNAREA": 0.73, "CSL1085LFP": 0.24}, 57.0, 0.81, 8.0),
    0.002: (1.65, {"DRNAREA": 0.73, "CSL1085LFP": 0.25}, 61.0, 0.79, 7.0),
}

_NM3_COEFFS: Dict[float, tuple] = {
    0.500: (0.80, {"DRNAREA": 0.64}, 78.0, 0.65, 4.0),
    0.200: (1.13, {"DRNAREA": 0.65}, 74.0, 0.67, 5.0),
    0.100: (1.36, {"DRNAREA": 0.66}, 71.0, 0.69, 5.0),
    0.040: (1.61, {"DRNAREA": 0.67}, 69.0, 0.70, 6.0),
    0.020: (1.77, {"DRNAREA": 0.67}, 68.0, 0.71, 6.0),
    0.010: (1.92, {"DRNAREA": 0.68}, 68.0, 0.71, 6.0),
    0.005: (2.06, {"DRNAREA": 0.68}, 70.0, 0.70, 6.0),
    0.002: (2.24, {"DRNAREA": 0.69}, 74.0, 0.68, 5.0),
}

_NM4_COEFFS: Dict[float, tuple] = {
    0.500: (0.22, {"DRNAREA": 0.66, "CSL1085LFP": 0.28}, 70.0, 0.72, 5.0),
    0.200: (0.55, {"DRNAREA": 0.67, "CSL1085LFP": 0.29}, 66.0, 0.74, 6.0),
    0.100: (0.76, {"DRNAREA": 0.67, "CSL1085LFP": 0.30}, 63.0, 0.75, 6.0),
    0.040: (1.00, {"DRNAREA": 0.68, "CSL1085LFP": 0.30}, 61.0, 0.76, 7.0),
    0.020: (1.15, {"DRNAREA": 0.68, "CSL1085LFP": 0.31}, 60.0, 0.77, 7.0),
    0.010: (1.29, {"DRNAREA": 0.69, "CSL1085LFP": 0.31}, 60.0, 0.77, 7.0),
    0.005: (1.42, {"DRNAREA": 0.69, "CSL1085LFP": 0.32}, 62.0, 0.76, 7.0),
    0.002: (1.59, {"DRNAREA": 0.70, "CSL1085LFP": 0.32}, 66.0, 0.74, 6.0),
}

_REGION_COEFFS = {
    NM_MOUNTAIN: _NM1_COEFFS,
    NM_TRANSITION: _NM2_COEFFS,
    NM_PLAINS: _NM3_COEFFS,
    NM_DESERT: _NM4_COEFFS,
}


def build_new_mexico_table() -> RegressionTable:
    """
    Build a :class:`RegressionTable` for all four New Mexico hydrologic regions.

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
    for region in NM_REGIONS:
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
