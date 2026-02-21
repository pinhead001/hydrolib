"""
hydrolib.regression.states.colorado - Colorado flood-frequency regression equations.

Four hydrologic regions reflecting Colorado's diverse physiography.
Equations follow the log10-linear GLS form::

    log10(Q_T) = intercept + b1*log10(DRNAREA) [+ b2*log10(predictor2) + ...]

Region definitions
------------------
``CO_FRONT_RANGE``  — Front Range Mountains — DRNAREA + ELEV + PRECIP
``CO_EASTERN``      — Eastern Plains — DRNAREA + CSL1085LFP
``CO_WESTERN``      — Western Slope (Colorado Plateau/Rockies) — DRNAREA + ELEV
``CO_SAN_LUIS``     — San Luis Valley / Upper Rio Grande — DRNAREA only

Usage
-----
::

    from hydrolib.regression.states.colorado import (
        CO_FRONT_RANGE, CO_EASTERN, CO_WESTERN, CO_SAN_LUIS,
        CO_REGIONS, build_colorado_table,
    )
    from hydrolib.regression.basin_chars import BasinCharacteristics, DRNAREA, ELEV, PRECIP

    table = build_colorado_table()
    basin = BasinCharacteristics(
        site_no="09105000",
        site_name="Crystal River above Avalanche Creek, CO",
        region=CO_WESTERN,
        predictors={DRNAREA: 191.0, ELEV: 9200.0},
    )
    q100 = table.estimate(basin, aep=0.01)

Note
----
Coefficient values embedded here are representative for illustration and
testing.  Use authoritative published values for engineering design.
Reference: Vaill, J.E. (2000), USGS WRIR 99-4190; and subsequent updates.
"""

from __future__ import annotations

from typing import Dict, Tuple

from hydrolib.regression.region import HydrologicRegion
from hydrolib.regression.regression_table import GlsEquation, RegressionTable

_PUB = "Vaill (2000), USGS WRIR 99-4190; https://doi.org/10.3133/wri994190"

# ---------------------------------------------------------------------------
# Region definitions
# ---------------------------------------------------------------------------

CO_FRONT_RANGE = HydrologicRegion(
    code="CO-1",
    label="Colorado Front Range Mountains",
    state="CO",
    required_predictors=("DRNAREA", "ELEV", "PRECIP"),
    publication=_PUB,
    notes="ELEV = mean basin elevation (ft); PRECIP = mean annual precipitation (in). "
    "Eastern slope of the Rockies from Wyoming to New Mexico border; "
    "snowmelt-dominated with intense convective summer storms.",
)

CO_EASTERN = HydrologicRegion(
    code="CO-2",
    label="Colorado Eastern Plains",
    state="CO",
    required_predictors=("DRNAREA", "CSL1085LFP"),
    publication=_PUB,
    notes="CSL1085LFP = 10-85 channel slope (ft/mi). Semi-arid High Plains "
    "east of the Front Range; convection-dominated flash flooding.",
)

CO_WESTERN = HydrologicRegion(
    code="CO-3",
    label="Colorado Western Slope",
    state="CO",
    required_predictors=("DRNAREA", "ELEV"),
    publication=_PUB,
    notes="ELEV = mean basin elevation (ft). Colorado Plateau and Western Slope "
    "Rockies draining to the Colorado River; snowmelt-dominated.",
)

CO_SAN_LUIS = HydrologicRegion(
    code="CO-4",
    label="Colorado San Luis Valley and Upper Rio Grande",
    state="CO",
    required_predictors=("DRNAREA",),
    publication=_PUB,
    notes="Drainage area is the sole predictor. Closed-basin and semi-arid "
    "San Luis Valley; upper Rio Grande basin in southern Colorado.",
)

CO_REGIONS: Tuple[HydrologicRegion, ...] = (CO_FRONT_RANGE, CO_EASTERN, CO_WESTERN, CO_SAN_LUIS)

# ---------------------------------------------------------------------------
# Regression coefficients (illustrative — use published values for design)
# ---------------------------------------------------------------------------
# Format: {aep: (intercept, {code: exponent}, sep_pct, pseudo_r2, eyr)}

_CO1_COEFFS: Dict[float, tuple] = {
    0.500: (-2.10, {"DRNAREA": 0.75, "ELEV": 0.52, "PRECIP": 0.70}, 55.0, 0.83, 9.0),
    0.200: (-1.78, {"DRNAREA": 0.76, "ELEV": 0.53, "PRECIP": 0.72}, 51.0, 0.85, 10.0),
    0.100: (-1.56, {"DRNAREA": 0.76, "ELEV": 0.54, "PRECIP": 0.74}, 48.0, 0.86, 11.0),
    0.040: (-1.32, {"DRNAREA": 0.77, "ELEV": 0.55, "PRECIP": 0.75}, 46.0, 0.87, 12.0),
    0.020: (-1.16, {"DRNAREA": 0.77, "ELEV": 0.55, "PRECIP": 0.76}, 45.0, 0.87, 12.0),
    0.010: (-1.01, {"DRNAREA": 0.78, "ELEV": 0.55, "PRECIP": 0.76}, 45.0, 0.87, 12.0),
    0.005: (-0.86, {"DRNAREA": 0.78, "ELEV": 0.56, "PRECIP": 0.77}, 47.0, 0.86, 11.0),
    0.002: (-0.65, {"DRNAREA": 0.79, "ELEV": 0.56, "PRECIP": 0.77}, 51.0, 0.85, 10.0),
}

_CO2_COEFFS: Dict[float, tuple] = {
    0.500: (0.90, {"DRNAREA": 0.65, "CSL1085LFP": 0.22}, 74.0, 0.68, 5.0),
    0.200: (1.22, {"DRNAREA": 0.66, "CSL1085LFP": 0.23}, 70.0, 0.70, 5.0),
    0.100: (1.43, {"DRNAREA": 0.66, "CSL1085LFP": 0.24}, 67.0, 0.72, 6.0),
    0.040: (1.67, {"DRNAREA": 0.67, "CSL1085LFP": 0.24}, 65.0, 0.73, 6.0),
    0.020: (1.83, {"DRNAREA": 0.67, "CSL1085LFP": 0.25}, 64.0, 0.74, 7.0),
    0.010: (1.97, {"DRNAREA": 0.68, "CSL1085LFP": 0.25}, 64.0, 0.74, 7.0),
    0.005: (2.11, {"DRNAREA": 0.68, "CSL1085LFP": 0.26}, 66.0, 0.73, 6.0),
    0.002: (2.30, {"DRNAREA": 0.69, "CSL1085LFP": 0.26}, 70.0, 0.71, 5.0),
}

_CO3_COEFFS: Dict[float, tuple] = {
    0.500: (-1.40, {"DRNAREA": 0.78, "ELEV": 0.58}, 52.0, 0.84, 9.0),
    0.200: (-1.10, {"DRNAREA": 0.78, "ELEV": 0.59}, 48.0, 0.86, 10.0),
    0.100: (-0.89, {"DRNAREA": 0.79, "ELEV": 0.60}, 46.0, 0.87, 11.0),
    0.040: (-0.65, {"DRNAREA": 0.79, "ELEV": 0.61}, 44.0, 0.88, 12.0),
    0.020: (-0.49, {"DRNAREA": 0.80, "ELEV": 0.61}, 43.0, 0.88, 12.0),
    0.010: (-0.35, {"DRNAREA": 0.80, "ELEV": 0.62}, 43.0, 0.88, 12.0),
    0.005: (-0.21, {"DRNAREA": 0.81, "ELEV": 0.62}, 45.0, 0.87, 11.0),
    0.002: (-0.01, {"DRNAREA": 0.81, "ELEV": 0.63}, 49.0, 0.86, 10.0),
}

_CO4_COEFFS: Dict[float, tuple] = {
    0.500: (0.68, {"DRNAREA": 0.66}, 76.0, 0.66, 4.0),
    0.200: (1.00, {"DRNAREA": 0.67}, 72.0, 0.68, 5.0),
    0.100: (1.22, {"DRNAREA": 0.68}, 69.0, 0.70, 5.0),
    0.040: (1.46, {"DRNAREA": 0.68}, 67.0, 0.71, 6.0),
    0.020: (1.62, {"DRNAREA": 0.69}, 66.0, 0.72, 6.0),
    0.010: (1.76, {"DRNAREA": 0.69}, 66.0, 0.72, 6.0),
    0.005: (1.90, {"DRNAREA": 0.70}, 68.0, 0.71, 6.0),
    0.002: (2.09, {"DRNAREA": 0.70}, 72.0, 0.69, 5.0),
}

_REGION_COEFFS = {
    CO_FRONT_RANGE: _CO1_COEFFS,
    CO_EASTERN: _CO2_COEFFS,
    CO_WESTERN: _CO3_COEFFS,
    CO_SAN_LUIS: _CO4_COEFFS,
}


def build_colorado_table() -> RegressionTable:
    """
    Build a :class:`RegressionTable` for all four Colorado hydrologic regions.

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
    for region in CO_REGIONS:
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
