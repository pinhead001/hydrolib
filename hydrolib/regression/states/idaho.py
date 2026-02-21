"""
hydrolib.regression.states.idaho - Idaho flood-frequency regression equations.

Four hydrologic regions reflecting Idaho's diverse hydrology.
Equations follow the log10-linear GLS form::

    log10(Q_T) = intercept + b1*log10(DRNAREA) [+ b2*log10(predictor2) + ...]

Region definitions
------------------
``ID_NORTH``   — Northern Mountains (Panhandle) — DRNAREA + PRECIP
``ID_CENTRAL`` — Central Mountains (Salmon/Clearwater) — DRNAREA + ELEV + PRECIP
``ID_SNAKE``   — Snake River Plain — DRNAREA + CSL1085LFP
``ID_SOUTH``   — Southern Idaho Highlands — DRNAREA only

Usage
-----
::

    from hydrolib.regression.states.idaho import (
        ID_NORTH, ID_CENTRAL, ID_SNAKE, ID_SOUTH,
        ID_REGIONS, build_idaho_table,
    )
    from hydrolib.regression.basin_chars import BasinCharacteristics, DRNAREA, ELEV, PRECIP

    table = build_idaho_table()
    basin = BasinCharacteristics(
        site_no="13206000",
        site_name="Boise River near Twin Springs, ID",
        region=ID_CENTRAL,
        predictors={DRNAREA: 1220.0, ELEV: 6800.0, PRECIP: 38.0},
    )
    q100 = table.estimate(basin, aep=0.01)

Note
----
Coefficient values embedded here are representative for illustration and
testing.  Use authoritative published values for engineering design.
Reference: Hortness, J.E., and Berenbrock, C. (2003), USGS SIR 2003-4278.
"""

from __future__ import annotations

from typing import Dict, Tuple

from hydrolib.regression.region import HydrologicRegion
from hydrolib.regression.regression_table import GlsEquation, RegressionTable

_PUB = "Hortness and Berenbrock (2003), USGS SIR 2003-4278; https://doi.org/10.3133/sir20034278"

# ---------------------------------------------------------------------------
# Region definitions
# ---------------------------------------------------------------------------

ID_NORTH = HydrologicRegion(
    code="ID-1",
    label="Idaho Northern Mountains (Panhandle)",
    state="ID",
    required_predictors=("DRNAREA", "PRECIP"),
    publication=_PUB,
    notes="PRECIP = mean annual precipitation (in). Northern Rockies of the "
    "Panhandle; maritime-influenced; high precipitation; perennial streams.",
)

ID_CENTRAL = HydrologicRegion(
    code="ID-2",
    label="Idaho Central Mountains (Salmon and Clearwater)",
    state="ID",
    required_predictors=("DRNAREA", "ELEV", "PRECIP"),
    publication=_PUB,
    notes="ELEV = mean basin elevation (ft); PRECIP = mean annual precipitation (in). "
    "Salmon River Mountains and Clearwater highlands; "
    "snowmelt-dominated; some of deepest gorges in North America.",
)

ID_SNAKE = HydrologicRegion(
    code="ID-3",
    label="Idaho Snake River Plain",
    state="ID",
    required_predictors=("DRNAREA", "CSL1085LFP"),
    publication=_PUB,
    notes="CSL1085LFP = 10-85 channel slope (ft/mi). Eastern Snake River Plain "
    "and southern basalt plateaus; irrigation-dominated hydrology.",
)

ID_SOUTH = HydrologicRegion(
    code="ID-4",
    label="Idaho Southern Highlands",
    state="ID",
    required_predictors=("DRNAREA",),
    publication=_PUB,
    notes="Drainage area is the sole predictor. Owyhee Plateau, Twin Falls, "
    "and semi-arid southern highlands; sparse, episodic runoff.",
)

ID_REGIONS: Tuple[HydrologicRegion, ...] = (ID_NORTH, ID_CENTRAL, ID_SNAKE, ID_SOUTH)

# ---------------------------------------------------------------------------
# Regression coefficients (illustrative — use published values for design)
# ---------------------------------------------------------------------------
# Format: {aep: (intercept, {code: exponent}, sep_pct, pseudo_r2, eyr)}

_ID1_COEFFS: Dict[float, tuple] = {
    0.500: (-0.20, {"DRNAREA": 0.78, "PRECIP": 0.68}, 48.0, 0.88, 12.0),
    0.200: (0.12, {"DRNAREA": 0.79, "PRECIP": 0.70}, 44.0, 0.89, 13.0),
    0.100: (0.33, {"DRNAREA": 0.79, "PRECIP": 0.71}, 42.0, 0.90, 14.0),
    0.040: (0.56, {"DRNAREA": 0.80, "PRECIP": 0.72}, 40.0, 0.91, 15.0),
    0.020: (0.72, {"DRNAREA": 0.80, "PRECIP": 0.73}, 39.0, 0.91, 15.0),
    0.010: (0.86, {"DRNAREA": 0.81, "PRECIP": 0.73}, 39.0, 0.91, 15.0),
    0.005: (0.99, {"DRNAREA": 0.81, "PRECIP": 0.74}, 41.0, 0.90, 14.0),
    0.002: (1.17, {"DRNAREA": 0.82, "PRECIP": 0.74}, 44.0, 0.89, 13.0),
}

_ID2_COEFFS: Dict[float, tuple] = {
    0.500: (-1.55, {"DRNAREA": 0.77, "ELEV": 0.50, "PRECIP": 0.65}, 50.0, 0.86, 10.0),
    0.200: (-1.23, {"DRNAREA": 0.78, "ELEV": 0.51, "PRECIP": 0.67}, 46.0, 0.87, 11.0),
    0.100: (-1.01, {"DRNAREA": 0.78, "ELEV": 0.52, "PRECIP": 0.68}, 44.0, 0.88, 12.0),
    0.040: (-0.77, {"DRNAREA": 0.79, "ELEV": 0.53, "PRECIP": 0.69}, 42.0, 0.89, 13.0),
    0.020: (-0.61, {"DRNAREA": 0.79, "ELEV": 0.53, "PRECIP": 0.70}, 41.0, 0.89, 13.0),
    0.010: (-0.46, {"DRNAREA": 0.80, "ELEV": 0.54, "PRECIP": 0.70}, 41.0, 0.89, 13.0),
    0.005: (-0.31, {"DRNAREA": 0.80, "ELEV": 0.54, "PRECIP": 0.71}, 43.0, 0.88, 12.0),
    0.002: (-0.11, {"DRNAREA": 0.81, "ELEV": 0.55, "PRECIP": 0.71}, 47.0, 0.87, 11.0),
}

_ID3_COEFFS: Dict[float, tuple] = {
    0.500: (0.55, {"DRNAREA": 0.74, "CSL1085LFP": 0.22}, 60.0, 0.80, 8.0),
    0.200: (0.87, {"DRNAREA": 0.74, "CSL1085LFP": 0.23}, 56.0, 0.82, 9.0),
    0.100: (1.07, {"DRNAREA": 0.75, "CSL1085LFP": 0.24}, 53.0, 0.83, 10.0),
    0.040: (1.30, {"DRNAREA": 0.75, "CSL1085LFP": 0.25}, 51.0, 0.84, 11.0),
    0.020: (1.45, {"DRNAREA": 0.76, "CSL1085LFP": 0.25}, 50.0, 0.85, 11.0),
    0.010: (1.58, {"DRNAREA": 0.76, "CSL1085LFP": 0.26}, 50.0, 0.85, 11.0),
    0.005: (1.71, {"DRNAREA": 0.77, "CSL1085LFP": 0.26}, 52.0, 0.84, 10.0),
    0.002: (1.88, {"DRNAREA": 0.77, "CSL1085LFP": 0.27}, 56.0, 0.82, 9.0),
}

_ID4_COEFFS: Dict[float, tuple] = {
    0.500: (0.62, {"DRNAREA": 0.69}, 72.0, 0.70, 5.0),
    0.200: (0.95, {"DRNAREA": 0.70}, 68.0, 0.72, 6.0),
    0.100: (1.16, {"DRNAREA": 0.71}, 65.0, 0.73, 6.0),
    0.040: (1.41, {"DRNAREA": 0.71}, 63.0, 0.74, 7.0),
    0.020: (1.57, {"DRNAREA": 0.72}, 62.0, 0.75, 7.0),
    0.010: (1.71, {"DRNAREA": 0.72}, 62.0, 0.75, 7.0),
    0.005: (1.85, {"DRNAREA": 0.73}, 64.0, 0.74, 7.0),
    0.002: (2.03, {"DRNAREA": 0.73}, 68.0, 0.72, 6.0),
}

_REGION_COEFFS = {
    ID_NORTH: _ID1_COEFFS,
    ID_CENTRAL: _ID2_COEFFS,
    ID_SNAKE: _ID3_COEFFS,
    ID_SOUTH: _ID4_COEFFS,
}


def build_idaho_table() -> RegressionTable:
    """
    Build a :class:`RegressionTable` for all four Idaho hydrologic regions.

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
    for region in ID_REGIONS:
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
