"""
hydrolib.regression.states.california - California flood-frequency regression equations.

Four hydrologic regions reflecting California's extreme climatic diversity.
Equations follow the log10-linear GLS form::

    log10(Q_T) = intercept + b1*log10(DRNAREA) [+ b2*log10(predictor2) + ...]

Region definitions
------------------
``CA_NORTH_COAST``   — North Coast Ranges and Klamath Mountains — DRNAREA + PRECIP
``CA_SIERRA``        — Sierra Nevada — DRNAREA + ELEV + PRECIP
``CA_CENTRAL``       — Central Valley and Foothills — DRNAREA + CSL1085LFP
``CA_SOUTH``         — Southern California (Transverse/Peninsular Ranges) — DRNAREA only

Usage
-----
::

    from hydrolib.regression.states.california import (
        CA_NORTH_COAST, CA_SIERRA, CA_CENTRAL, CA_SOUTH,
        CA_REGIONS, build_california_table,
    )
    from hydrolib.regression.basin_chars import BasinCharacteristics, DRNAREA, ELEV, PRECIP

    table = build_california_table()
    basin = BasinCharacteristics(
        site_no="11522500",
        site_name="Scott River near Fort Jones, CA",
        region=CA_NORTH_COAST,
        predictors={DRNAREA: 633.0, PRECIP: 42.0},
    )
    q100 = table.estimate(basin, aep=0.01)

Note
----
Coefficient values embedded here are representative for illustration and
testing.  Use authoritative published values for engineering design.
Reference: Gotvald, A.J., and others (2012), USGS SIR 2012-5113; and
Waananen, A.O., and Crippen, J.R. (1977), USGS WSP 1962-C.
"""

from __future__ import annotations

from typing import Dict, Tuple

from hydrolib.regression.region import HydrologicRegion
from hydrolib.regression.regression_table import GlsEquation, RegressionTable

_PUB = "Gotvald and others (2012), USGS SIR 2012-5113; https://doi.org/10.3133/sir20125113"

# ---------------------------------------------------------------------------
# Region definitions
# ---------------------------------------------------------------------------

CA_NORTH_COAST = HydrologicRegion(
    code="CA-1",
    label="California North Coast Ranges and Klamath",
    state="CA",
    required_predictors=("DRNAREA", "PRECIP"),
    publication=_PUB,
    notes="PRECIP = mean annual precipitation (in). Northern Coast Ranges "
    "and Klamath/Siskiyou Mountains; winter rainfall-dominated flooding.",
)

CA_SIERRA = HydrologicRegion(
    code="CA-2",
    label="California Sierra Nevada",
    state="CA",
    required_predictors=("DRNAREA", "ELEV", "PRECIP"),
    publication=_PUB,
    notes="ELEV = mean basin elevation (ft); PRECIP = mean annual precipitation (in). "
    "Sierra Nevada from Cascades to Tehachapi; snowmelt and rain-on-snow "
    "events; largest spring runoff in western CONUS.",
)

CA_CENTRAL = HydrologicRegion(
    code="CA-3",
    label="California Central Valley and Foothills",
    state="CA",
    required_predictors=("DRNAREA", "CSL1085LFP"),
    publication=_PUB,
    notes="CSL1085LFP = 10-85 channel slope (ft/mi). Sacramento and San "
    "Joaquin Valley tributaries; mixed winter rainfall and snowmelt.",
)

CA_SOUTH = HydrologicRegion(
    code="CA-4",
    label="California Southern (Transverse and Peninsular Ranges)",
    state="CA",
    required_predictors=("DRNAREA",),
    publication=_PUB,
    notes="Drainage area is the sole predictor. Transverse and Peninsular "
    "Ranges of southern California; Mediterranean climate; highly variable "
    "floods driven by atmospheric rivers; very high SEP%.",
)

CA_REGIONS: Tuple[HydrologicRegion, ...] = (CA_NORTH_COAST, CA_SIERRA, CA_CENTRAL, CA_SOUTH)

# ---------------------------------------------------------------------------
# Regression coefficients (illustrative — use published values for design)
# ---------------------------------------------------------------------------
# Format: {aep: (intercept, {code: exponent}, sep_pct, pseudo_r2, eyr)}

_CA1_COEFFS: Dict[float, tuple] = {
    0.500: (0.18, {"DRNAREA": 0.82, "PRECIP": 0.62}, 45.0, 0.90, 13.0),
    0.200: (0.48, {"DRNAREA": 0.83, "PRECIP": 0.64}, 41.0, 0.91, 15.0),
    0.100: (0.68, {"DRNAREA": 0.83, "PRECIP": 0.65}, 39.0, 0.92, 16.0),
    0.040: (0.91, {"DRNAREA": 0.84, "PRECIP": 0.66}, 37.0, 0.93, 17.0),
    0.020: (1.06, {"DRNAREA": 0.84, "PRECIP": 0.67}, 37.0, 0.93, 17.0),
    0.010: (1.19, {"DRNAREA": 0.85, "PRECIP": 0.67}, 37.0, 0.93, 17.0),
    0.005: (1.32, {"DRNAREA": 0.85, "PRECIP": 0.68}, 39.0, 0.92, 16.0),
    0.002: (1.50, {"DRNAREA": 0.86, "PRECIP": 0.68}, 43.0, 0.91, 15.0),
}

_CA2_COEFFS: Dict[float, tuple] = {
    0.500: (-1.90, {"DRNAREA": 0.78, "ELEV": 0.55, "PRECIP": 0.68}, 50.0, 0.86, 10.0),
    0.200: (-1.58, {"DRNAREA": 0.79, "ELEV": 0.56, "PRECIP": 0.70}, 46.0, 0.88, 11.0),
    0.100: (-1.36, {"DRNAREA": 0.79, "ELEV": 0.57, "PRECIP": 0.72}, 44.0, 0.89, 12.0),
    0.040: (-1.11, {"DRNAREA": 0.80, "ELEV": 0.57, "PRECIP": 0.73}, 42.0, 0.90, 13.0),
    0.020: (-0.95, {"DRNAREA": 0.80, "ELEV": 0.58, "PRECIP": 0.74}, 41.0, 0.90, 14.0),
    0.010: (-0.80, {"DRNAREA": 0.81, "ELEV": 0.58, "PRECIP": 0.75}, 41.0, 0.90, 14.0),
    0.005: (-0.65, {"DRNAREA": 0.81, "ELEV": 0.59, "PRECIP": 0.75}, 43.0, 0.89, 13.0),
    0.002: (-0.44, {"DRNAREA": 0.82, "ELEV": 0.59, "PRECIP": 0.76}, 47.0, 0.88, 12.0),
}

_CA3_COEFFS: Dict[float, tuple] = {
    0.500: (1.30, {"DRNAREA": 0.77, "CSL1085LFP": 0.30}, 52.0, 0.84, 9.0),
    0.200: (1.62, {"DRNAREA": 0.78, "CSL1085LFP": 0.31}, 48.0, 0.86, 10.0),
    0.100: (1.82, {"DRNAREA": 0.78, "CSL1085LFP": 0.32}, 46.0, 0.87, 11.0),
    0.040: (2.05, {"DRNAREA": 0.79, "CSL1085LFP": 0.33}, 44.0, 0.88, 12.0),
    0.020: (2.20, {"DRNAREA": 0.79, "CSL1085LFP": 0.33}, 43.0, 0.89, 12.0),
    0.010: (2.33, {"DRNAREA": 0.80, "CSL1085LFP": 0.34}, 43.0, 0.89, 12.0),
    0.005: (2.46, {"DRNAREA": 0.80, "CSL1085LFP": 0.35}, 45.0, 0.88, 12.0),
    0.002: (2.64, {"DRNAREA": 0.81, "CSL1085LFP": 0.35}, 49.0, 0.87, 11.0),
}

_CA4_COEFFS: Dict[float, tuple] = {
    0.500: (1.65, {"DRNAREA": 0.70}, 82.0, 0.62, 3.0),
    0.200: (2.00, {"DRNAREA": 0.71}, 78.0, 0.64, 4.0),
    0.100: (2.23, {"DRNAREA": 0.72}, 75.0, 0.66, 4.0),
    0.040: (2.49, {"DRNAREA": 0.73}, 73.0, 0.67, 5.0),
    0.020: (2.66, {"DRNAREA": 0.73}, 72.0, 0.68, 5.0),
    0.010: (2.81, {"DRNAREA": 0.74}, 72.0, 0.68, 5.0),
    0.005: (2.96, {"DRNAREA": 0.74}, 75.0, 0.67, 5.0),
    0.002: (3.15, {"DRNAREA": 0.75}, 80.0, 0.65, 4.0),
}

_REGION_COEFFS = {
    CA_NORTH_COAST: _CA1_COEFFS,
    CA_SIERRA: _CA2_COEFFS,
    CA_CENTRAL: _CA3_COEFFS,
    CA_SOUTH: _CA4_COEFFS,
}


def build_california_table() -> RegressionTable:
    """
    Build a :class:`RegressionTable` for all four California hydrologic regions.

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
    for region in CA_REGIONS:
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
