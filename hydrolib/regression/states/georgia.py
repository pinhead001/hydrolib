"""
hydrolib.regression.states.georgia - Georgia flood-frequency regression equations.

Four physiographic hydrologic regions based on USGS regression studies for
Georgia streams.  Equations follow the log10-linear form::

    log10(Q_T) = intercept + b1*log10(DRNAREA) [+ b2*log10(predictor2)]

Region definitions
------------------
``GA_BLUE_RIDGE``   — Blue Ridge Mountains — DRNAREA + ELEV
``GA_VALLEY_RIDGE`` — Valley and Ridge — DRNAREA + CSL1085LFP
``GA_PIEDMONT``     — Piedmont — DRNAREA only
``GA_COASTAL``      — Coastal Plain — DRNAREA + PRECIP

Usage
-----
::

    from hydrolib.regression.states.georgia import (
        GA_BLUE_RIDGE, GA_VALLEY_RIDGE, GA_PIEDMONT, GA_COASTAL,
        GA_REGIONS, build_georgia_table,
    )
    from hydrolib.regression.basin_chars import BasinCharacteristics, DRNAREA, ELEV

    table = build_georgia_table()
    basin = BasinCharacteristics(
        site_no="02333500",
        site_name="Chattahoochee River near Helen, GA",
        region=GA_BLUE_RIDGE,
        predictors={DRNAREA: 47.0, ELEV: 2650.0},
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

_PUB = "USGS Georgia StreamStats flood-frequency regression equations"

# ---------------------------------------------------------------------------
# Region definitions
# ---------------------------------------------------------------------------

GA_BLUE_RIDGE = HydrologicRegion(
    code="GA-1",
    label="Georgia Blue Ridge Mountains",
    state="GA",
    required_predictors=("DRNAREA", "ELEV"),
    publication=_PUB,
    notes="ELEV = mean basin elevation (ft). Steep, forested headwater basins "
    "in the Blue Ridge physiographic province of north Georgia.",
)

GA_VALLEY_RIDGE = HydrologicRegion(
    code="GA-2",
    label="Georgia Valley and Ridge",
    state="GA",
    required_predictors=("DRNAREA", "CSL1085LFP"),
    publication=_PUB,
    notes="CSL1085LFP = 10-85 channel slope (ft/mi). "
    "Valley and Ridge physiographic province, northwest Georgia.",
)

GA_PIEDMONT = HydrologicRegion(
    code="GA-3",
    label="Georgia Piedmont",
    state="GA",
    required_predictors=("DRNAREA",),
    publication=_PUB,
    notes="Drainage area is the sole predictor. Piedmont physiographic "
    "province; central Georgia from Atlanta southward.",
)

GA_COASTAL = HydrologicRegion(
    code="GA-4",
    label="Georgia Coastal Plain",
    state="GA",
    required_predictors=("DRNAREA", "PRECIP"),
    publication=_PUB,
    notes="PRECIP = mean annual precipitation (in). "
    "Coastal Plain physiographic province, south Georgia.",
)

GA_REGIONS: Tuple[HydrologicRegion, ...] = (
    GA_BLUE_RIDGE,
    GA_VALLEY_RIDGE,
    GA_PIEDMONT,
    GA_COASTAL,
)

# ---------------------------------------------------------------------------
# Regression coefficients
# ---------------------------------------------------------------------------
# Format: {aep: (intercept, {code: exponent}, sep_pct, pseudo_r2, eyr)}

_GA1_COEFFS: Dict[float, tuple] = {
    0.500: (0.30, {"DRNAREA": 0.79, "ELEV": 0.28}, 48.0, 0.87, 11.0),
    0.200: (0.56, {"DRNAREA": 0.80, "ELEV": 0.29}, 44.0, 0.88, 12.0),
    0.100: (0.73, {"DRNAREA": 0.81, "ELEV": 0.32}, 41.0, 0.89, 13.0),
    0.040: (0.90, {"DRNAREA": 0.81, "ELEV": 0.34}, 39.0, 0.90, 14.0),
    0.020: (1.02, {"DRNAREA": 0.82, "ELEV": 0.35}, 38.0, 0.90, 14.0),
    0.010: (1.14, {"DRNAREA": 0.82, "ELEV": 0.35}, 38.0, 0.90, 14.0),
    0.005: (1.26, {"DRNAREA": 0.83, "ELEV": 0.35}, 40.0, 0.89, 13.0),
    0.002: (1.43, {"DRNAREA": 0.84, "ELEV": 0.36}, 43.0, 0.88, 12.0),
}

_GA2_COEFFS: Dict[float, tuple] = {
    0.500: (1.50, {"DRNAREA": 0.72, "CSL1085LFP": 0.28}, 50.0, 0.86, 10.0),
    0.200: (1.83, {"DRNAREA": 0.72, "CSL1085LFP": 0.29}, 46.0, 0.87, 11.0),
    0.100: (2.03, {"DRNAREA": 0.73, "CSL1085LFP": 0.30}, 43.0, 0.88, 12.0),
    0.040: (2.25, {"DRNAREA": 0.73, "CSL1085LFP": 0.31}, 41.0, 0.89, 13.0),
    0.020: (2.39, {"DRNAREA": 0.74, "CSL1085LFP": 0.32}, 40.0, 0.89, 13.0),
    0.010: (2.52, {"DRNAREA": 0.74, "CSL1085LFP": 0.32}, 40.0, 0.89, 13.0),
    0.005: (2.65, {"DRNAREA": 0.75, "CSL1085LFP": 0.32}, 42.0, 0.88, 12.0),
    0.002: (2.83, {"DRNAREA": 0.75, "CSL1085LFP": 0.33}, 45.0, 0.87, 11.0),
}

_GA3_COEFFS: Dict[float, tuple] = {
    0.500: (1.58, {"DRNAREA": 0.75}, 56.0, 0.82, 9.0),
    0.200: (1.92, {"DRNAREA": 0.76}, 51.0, 0.84, 10.0),
    0.100: (2.13, {"DRNAREA": 0.77}, 48.0, 0.85, 11.0),
    0.040: (2.37, {"DRNAREA": 0.77}, 46.0, 0.86, 12.0),
    0.020: (2.52, {"DRNAREA": 0.77}, 45.0, 0.86, 12.0),
    0.010: (2.65, {"DRNAREA": 0.77}, 44.0, 0.86, 12.0),
    0.005: (2.78, {"DRNAREA": 0.78}, 46.0, 0.85, 11.0),
    0.002: (2.96, {"DRNAREA": 0.78}, 49.0, 0.84, 10.0),
}

_GA4_COEFFS: Dict[float, tuple] = {
    0.500: (0.62, {"DRNAREA": 0.78, "PRECIP": 0.32}, 58.0, 0.80, 8.0),
    0.200: (0.92, {"DRNAREA": 0.78, "PRECIP": 0.35}, 53.0, 0.82, 9.0),
    0.100: (1.10, {"DRNAREA": 0.79, "PRECIP": 0.37}, 50.0, 0.83, 10.0),
    0.040: (1.31, {"DRNAREA": 0.79, "PRECIP": 0.38}, 48.0, 0.84, 11.0),
    0.020: (1.44, {"DRNAREA": 0.79, "PRECIP": 0.39}, 47.0, 0.84, 11.0),
    0.010: (1.56, {"DRNAREA": 0.79, "PRECIP": 0.40}, 47.0, 0.84, 11.0),
    0.005: (1.68, {"DRNAREA": 0.80, "PRECIP": 0.40}, 49.0, 0.83, 10.0),
    0.002: (1.84, {"DRNAREA": 0.80, "PRECIP": 0.41}, 53.0, 0.82, 9.0),
}

_REGION_COEFFS = {
    GA_BLUE_RIDGE: _GA1_COEFFS,
    GA_VALLEY_RIDGE: _GA2_COEFFS,
    GA_PIEDMONT: _GA3_COEFFS,
    GA_COASTAL: _GA4_COEFFS,
}


def build_georgia_table() -> RegressionTable:
    """
    Build a :class:`RegressionTable` for all four Georgia hydrologic regions.

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
    for region in GA_REGIONS:
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
