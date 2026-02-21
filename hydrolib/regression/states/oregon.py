"""
hydrolib.regression.states.oregon - Oregon flood-frequency regression equations.

Four hydrologic regions reflecting Oregon's rainfall-to-snowmelt gradient.
Equations follow the log10-linear GLS form::

    log10(Q_T) = intercept + b1*log10(DRNAREA) [+ b2*log10(predictor2) + ...]

Region definitions
------------------
``OR_COAST``    — Coast Range and Klamath Mountains — DRNAREA + PRECIP
``OR_WILLAMETTE`` — Willamette Valley and Foothills — DRNAREA + CSL1085LFP
``OR_CASCADE``  — Western Cascades and High Cascades — DRNAREA + ELEV + PRECIP
``OR_EASTERN``  — Eastern Oregon (Basin and Range / High Desert) — DRNAREA only

Usage
-----
::

    from hydrolib.regression.states.oregon import (
        OR_COAST, OR_WILLAMETTE, OR_CASCADE, OR_EASTERN,
        OR_REGIONS, build_oregon_table,
    )
    from hydrolib.regression.basin_chars import BasinCharacteristics, DRNAREA, PRECIP

    table = build_oregon_table()
    basin = BasinCharacteristics(
        site_no="14306500",
        site_name="Coquille River near Coquille, OR",
        region=OR_COAST,
        predictors={DRNAREA: 766.0, PRECIP: 65.0},
    )
    q100 = table.estimate(basin, aep=0.01)

Note
----
Coefficient values embedded here are representative for illustration and
testing.  Use authoritative published values for engineering design.
Reference: Risley, J.C., and others (2008), USGS SIR 2008-5126.
"""

from __future__ import annotations

from typing import Dict, Tuple

from hydrolib.regression.region import HydrologicRegion
from hydrolib.regression.regression_table import GlsEquation, RegressionTable

_PUB = "Risley and others (2008), USGS SIR 2008-5126; https://doi.org/10.3133/sir20085126"

# ---------------------------------------------------------------------------
# Region definitions
# ---------------------------------------------------------------------------

OR_COAST = HydrologicRegion(
    code="OR-1",
    label="Oregon Coast Range and Klamath Mountains",
    state="OR",
    required_predictors=("DRNAREA", "PRECIP"),
    publication=_PUB,
    notes="PRECIP = mean annual precipitation (in). Oregon Coast Range and "
    "Klamath Mountains; maritime climate; rainfall-dominated flooding.",
)

OR_WILLAMETTE = HydrologicRegion(
    code="OR-2",
    label="Oregon Willamette Valley and Foothills",
    state="OR",
    required_predictors=("DRNAREA", "CSL1085LFP"),
    publication=_PUB,
    notes="CSL1085LFP = 10-85 channel slope (ft/mi). Willamette Valley floor "
    "and adjacent foothills; mixed rainfall and snowmelt events.",
)

OR_CASCADE = HydrologicRegion(
    code="OR-3",
    label="Oregon Western and High Cascades",
    state="OR",
    required_predictors=("DRNAREA", "ELEV", "PRECIP"),
    publication=_PUB,
    notes="ELEV = mean basin elevation (ft); PRECIP = mean annual precipitation (in). "
    "Cascade Range from Columbia River to California border; "
    "rain-on-snow and snowmelt flood generation.",
)

OR_EASTERN = HydrologicRegion(
    code="OR-4",
    label="Oregon Eastern (Basin and Range / High Desert)",
    state="OR",
    required_predictors=("DRNAREA",),
    publication=_PUB,
    notes="Drainage area is the sole predictor. Great Basin and Columbia "
    "Plateau of eastern Oregon; semi-arid; ephemeral drainages.",
)

OR_REGIONS: Tuple[HydrologicRegion, ...] = (OR_COAST, OR_WILLAMETTE, OR_CASCADE, OR_EASTERN)

# ---------------------------------------------------------------------------
# Regression coefficients (illustrative — use published values for design)
# ---------------------------------------------------------------------------
# Format: {aep: (intercept, {code: exponent}, sep_pct, pseudo_r2, eyr)}

_OR1_COEFFS: Dict[float, tuple] = {
    0.500: (0.28, {"DRNAREA": 0.83, "PRECIP": 0.58}, 43.0, 0.91, 14.0),
    0.200: (0.57, {"DRNAREA": 0.83, "PRECIP": 0.60}, 39.0, 0.92, 16.0),
    0.100: (0.77, {"DRNAREA": 0.84, "PRECIP": 0.61}, 37.0, 0.93, 17.0),
    0.040: (0.99, {"DRNAREA": 0.84, "PRECIP": 0.63}, 35.0, 0.93, 18.0),
    0.020: (1.13, {"DRNAREA": 0.85, "PRECIP": 0.63}, 35.0, 0.93, 18.0),
    0.010: (1.26, {"DRNAREA": 0.85, "PRECIP": 0.64}, 35.0, 0.93, 18.0),
    0.005: (1.39, {"DRNAREA": 0.86, "PRECIP": 0.65}, 37.0, 0.93, 17.0),
    0.002: (1.56, {"DRNAREA": 0.86, "PRECIP": 0.65}, 40.0, 0.92, 16.0),
}

_OR2_COEFFS: Dict[float, tuple] = {
    0.500: (1.42, {"DRNAREA": 0.78, "CSL1085LFP": 0.26}, 50.0, 0.86, 10.0),
    0.200: (1.73, {"DRNAREA": 0.79, "CSL1085LFP": 0.27}, 46.0, 0.88, 12.0),
    0.100: (1.93, {"DRNAREA": 0.79, "CSL1085LFP": 0.28}, 44.0, 0.89, 13.0),
    0.040: (2.15, {"DRNAREA": 0.80, "CSL1085LFP": 0.29}, 42.0, 0.90, 14.0),
    0.020: (2.30, {"DRNAREA": 0.80, "CSL1085LFP": 0.29}, 41.0, 0.90, 14.0),
    0.010: (2.43, {"DRNAREA": 0.81, "CSL1085LFP": 0.30}, 41.0, 0.90, 14.0),
    0.005: (2.56, {"DRNAREA": 0.81, "CSL1085LFP": 0.30}, 43.0, 0.89, 13.0),
    0.002: (2.73, {"DRNAREA": 0.82, "CSL1085LFP": 0.31}, 47.0, 0.88, 12.0),
}

_OR3_COEFFS: Dict[float, tuple] = {
    0.500: (-0.70, {"DRNAREA": 0.80, "ELEV": 0.48, "PRECIP": 0.60}, 48.0, 0.87, 11.0),
    0.200: (-0.39, {"DRNAREA": 0.81, "ELEV": 0.49, "PRECIP": 0.62}, 44.0, 0.89, 12.0),
    0.100: (-0.18, {"DRNAREA": 0.81, "ELEV": 0.50, "PRECIP": 0.63}, 42.0, 0.90, 13.0),
    0.040: (0.06, {"DRNAREA": 0.82, "ELEV": 0.50, "PRECIP": 0.64}, 40.0, 0.91, 14.0),
    0.020: (0.22, {"DRNAREA": 0.82, "ELEV": 0.51, "PRECIP": 0.65}, 40.0, 0.91, 14.0),
    0.010: (0.36, {"DRNAREA": 0.83, "ELEV": 0.51, "PRECIP": 0.65}, 40.0, 0.91, 14.0),
    0.005: (0.50, {"DRNAREA": 0.83, "ELEV": 0.52, "PRECIP": 0.66}, 42.0, 0.90, 13.0),
    0.002: (0.68, {"DRNAREA": 0.84, "ELEV": 0.52, "PRECIP": 0.66}, 46.0, 0.89, 12.0),
}

_OR4_COEFFS: Dict[float, tuple] = {
    0.500: (0.75, {"DRNAREA": 0.68}, 74.0, 0.68, 5.0),
    0.200: (1.08, {"DRNAREA": 0.69}, 70.0, 0.70, 6.0),
    0.100: (1.29, {"DRNAREA": 0.70}, 67.0, 0.72, 6.0),
    0.040: (1.53, {"DRNAREA": 0.70}, 65.0, 0.73, 7.0),
    0.020: (1.69, {"DRNAREA": 0.71}, 64.0, 0.74, 7.0),
    0.010: (1.83, {"DRNAREA": 0.71}, 64.0, 0.74, 7.0),
    0.005: (1.97, {"DRNAREA": 0.72}, 66.0, 0.73, 7.0),
    0.002: (2.15, {"DRNAREA": 0.72}, 70.0, 0.71, 6.0),
}

_REGION_COEFFS = {
    OR_COAST: _OR1_COEFFS,
    OR_WILLAMETTE: _OR2_COEFFS,
    OR_CASCADE: _OR3_COEFFS,
    OR_EASTERN: _OR4_COEFFS,
}


def build_oregon_table() -> RegressionTable:
    """
    Build a :class:`RegressionTable` for all four Oregon hydrologic regions.

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
    for region in OR_REGIONS:
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
