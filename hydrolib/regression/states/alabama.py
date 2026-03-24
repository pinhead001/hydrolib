"""
hydrolib.regression.states.alabama - Alabama flood-frequency regression equations.

Four physiographic hydrologic regions based on USGS SIR 2014-5010.
Equations follow the log10-linear GLS form::

    log10(Q_T) = intercept + b1*log10(DRNAREA) [+ b2*log10(CSL1085LFP)]

Region definitions
------------------
``AL_APPALACHIAN`` — Appalachian Plateaus — DRNAREA + CSL1085LFP
``AL_VALLEY_RIDGE`` — Valley and Ridge — DRNAREA + CSL1085LFP
``AL_PIEDMONT``    — Piedmont — DRNAREA + CSL1085LFP
``AL_COASTAL``     — Coastal Plain — DRNAREA only

Usage
-----
::

    from hydrolib.regression.states.alabama import (
        AL_APPALACHIAN, AL_VALLEY_RIDGE, AL_PIEDMONT, AL_COASTAL,
        AL_REGIONS, build_alabama_table,
    )
    from hydrolib.regression.basin_chars import BasinCharacteristics, DRNAREA, CSL1085LFP

    table = build_alabama_table()
    basin = BasinCharacteristics(
        site_no="02430000",
        site_name="Bear Creek near Hackleburg, AL",
        region=AL_APPALACHIAN,
        predictors={DRNAREA: 232.0, CSL1085LFP: 8.2},
    )
    q100 = table.estimate(basin, aep=0.01)

Note
----
Coefficient values embedded here are representative for illustration and
testing.  Use authoritative published values for engineering design.
Reference: Southard, R.E., and Veilleux, A.G. (2014), USGS SIR 2014-5010.
"""

from __future__ import annotations

from typing import Dict, Tuple

from hydrolib.regression.region import HydrologicRegion
from hydrolib.regression.regression_table import GlsEquation, RegressionTable

_PUB = "Southard and Veilleux (2014), USGS SIR 2014-5010; https://doi.org/10.3133/sir20145010"

# ---------------------------------------------------------------------------
# Region definitions
# ---------------------------------------------------------------------------

AL_APPALACHIAN = HydrologicRegion(
    code="AL-1",
    label="Alabama Appalachian Plateaus",
    state="AL",
    required_predictors=("DRNAREA", "CSL1085LFP"),
    publication=_PUB,
    notes="CSL1085LFP = 10-85 channel slope (ft/mi). "
    "Cumberland Plateau and Warrior coal basin; "
    "sandstone-dominated watersheds with relatively high relief.",
)

AL_VALLEY_RIDGE = HydrologicRegion(
    code="AL-2",
    label="Alabama Valley and Ridge",
    state="AL",
    required_predictors=("DRNAREA", "CSL1085LFP"),
    publication=_PUB,
    notes="CSL1085LFP = 10-85 channel slope (ft/mi). "
    "Valley and Ridge physiographic province in northeast Alabama; "
    "limestone and shale valleys with karst features common.",
)

AL_PIEDMONT = HydrologicRegion(
    code="AL-3",
    label="Alabama Piedmont",
    state="AL",
    required_predictors=("DRNAREA", "CSL1085LFP"),
    publication=_PUB,
    notes="CSL1085LFP = 10-85 channel slope (ft/mi). "
    "Piedmont physiographic province in east-central Alabama; "
    "metamorphic and igneous rocks; moderately steep terrain.",
)

AL_COASTAL = HydrologicRegion(
    code="AL-4",
    label="Alabama Coastal Plain",
    state="AL",
    required_predictors=("DRNAREA",),
    publication=_PUB,
    notes="Drainage area is the sole predictor. "
    "Gulf Coastal Plain of southern Alabama; "
    "low relief, high rainfall, swampy terrain.",
)

AL_REGIONS: Tuple[HydrologicRegion, ...] = (
    AL_APPALACHIAN,
    AL_VALLEY_RIDGE,
    AL_PIEDMONT,
    AL_COASTAL,
)

# ---------------------------------------------------------------------------
# Regression coefficients (illustrative — use published values for design)
# ---------------------------------------------------------------------------
# Format: {aep: (intercept, {code: exponent}, sep_pct, pseudo_r2, eyr)}

_AL1_COEFFS: Dict[float, tuple] = {
    0.500: (1.42, {"DRNAREA": 0.74, "CSL1085LFP": 0.36}, 46.0, 0.89, 13.0),
    0.200: (1.74, {"DRNAREA": 0.74, "CSL1085LFP": 0.37}, 42.0, 0.90, 14.0),
    0.100: (1.95, {"DRNAREA": 0.74, "CSL1085LFP": 0.38}, 40.0, 0.91, 15.0),
    0.040: (2.19, {"DRNAREA": 0.75, "CSL1085LFP": 0.38}, 38.0, 0.91, 16.0),
    0.020: (2.35, {"DRNAREA": 0.75, "CSL1085LFP": 0.39}, 37.0, 0.91, 16.0),
    0.010: (2.49, {"DRNAREA": 0.75, "CSL1085LFP": 0.39}, 37.0, 0.91, 16.0),
    0.005: (2.63, {"DRNAREA": 0.76, "CSL1085LFP": 0.39}, 38.0, 0.90, 15.0),
    0.002: (2.81, {"DRNAREA": 0.76, "CSL1085LFP": 0.40}, 41.0, 0.89, 14.0),
}

_AL2_COEFFS: Dict[float, tuple] = {
    0.500: (1.35, {"DRNAREA": 0.77, "CSL1085LFP": 0.30}, 50.0, 0.87, 11.0),
    0.200: (1.67, {"DRNAREA": 0.77, "CSL1085LFP": 0.31}, 46.0, 0.88, 12.0),
    0.100: (1.88, {"DRNAREA": 0.77, "CSL1085LFP": 0.32}, 44.0, 0.89, 13.0),
    0.040: (2.13, {"DRNAREA": 0.78, "CSL1085LFP": 0.32}, 42.0, 0.89, 13.0),
    0.020: (2.29, {"DRNAREA": 0.78, "CSL1085LFP": 0.33}, 41.0, 0.89, 14.0),
    0.010: (2.43, {"DRNAREA": 0.78, "CSL1085LFP": 0.33}, 41.0, 0.89, 14.0),
    0.005: (2.57, {"DRNAREA": 0.79, "CSL1085LFP": 0.33}, 43.0, 0.88, 13.0),
    0.002: (2.76, {"DRNAREA": 0.79, "CSL1085LFP": 0.34}, 46.0, 0.87, 12.0),
}

_AL3_COEFFS: Dict[float, tuple] = {
    0.500: (1.50, {"DRNAREA": 0.73, "CSL1085LFP": 0.33}, 52.0, 0.85, 10.0),
    0.200: (1.83, {"DRNAREA": 0.73, "CSL1085LFP": 0.34}, 48.0, 0.86, 11.0),
    0.100: (2.04, {"DRNAREA": 0.73, "CSL1085LFP": 0.35}, 46.0, 0.87, 12.0),
    0.040: (2.29, {"DRNAREA": 0.74, "CSL1085LFP": 0.35}, 44.0, 0.87, 12.0),
    0.020: (2.45, {"DRNAREA": 0.74, "CSL1085LFP": 0.36}, 43.0, 0.88, 13.0),
    0.010: (2.59, {"DRNAREA": 0.74, "CSL1085LFP": 0.36}, 43.0, 0.88, 13.0),
    0.005: (2.73, {"DRNAREA": 0.75, "CSL1085LFP": 0.36}, 45.0, 0.87, 12.0),
    0.002: (2.92, {"DRNAREA": 0.75, "CSL1085LFP": 0.37}, 48.0, 0.86, 11.0),
}

_AL4_COEFFS: Dict[float, tuple] = {
    0.500: (1.55, {"DRNAREA": 0.82}, 54.0, 0.84, 10.0),
    0.200: (1.89, {"DRNAREA": 0.82}, 50.0, 0.85, 11.0),
    0.100: (2.11, {"DRNAREA": 0.83}, 47.0, 0.86, 12.0),
    0.040: (2.36, {"DRNAREA": 0.83}, 45.0, 0.87, 12.0),
    0.020: (2.52, {"DRNAREA": 0.83}, 44.0, 0.87, 13.0),
    0.010: (2.67, {"DRNAREA": 0.83}, 44.0, 0.87, 13.0),
    0.005: (2.81, {"DRNAREA": 0.84}, 46.0, 0.86, 12.0),
    0.002: (3.00, {"DRNAREA": 0.84}, 50.0, 0.85, 11.0),
}

_REGION_COEFFS = {
    AL_APPALACHIAN: _AL1_COEFFS,
    AL_VALLEY_RIDGE: _AL2_COEFFS,
    AL_PIEDMONT: _AL3_COEFFS,
    AL_COASTAL: _AL4_COEFFS,
}


def build_alabama_table() -> RegressionTable:
    """
    Build a :class:`RegressionTable` for all four Alabama hydrologic regions.

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
    for region in AL_REGIONS:
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
