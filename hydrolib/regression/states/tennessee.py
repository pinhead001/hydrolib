"""
hydrolib.regression.states.tennessee - Tennessee flood-frequency regression equations.

Regions and equation structure from:
    Ladd, D.E., and Ensminger, P.A., 2025, Estimating the magnitude and
    frequency of floods at ungaged locations on streams in Tennessee through
    the 2013 water year: U.S. Geological Survey Scientific Investigations
    Report 2024-5130, 19 p., https://doi.org/10.3133/sir20245130.

Four hydrologic areas are defined based on physiographic and hydrologic
similarity.  Regression equations are of the log10-linear form::

    log10(Q_T) = intercept + b1 * log10(DRNAREA) [+ b2 * log10(predictor2)]

Region definitions
------------------
``TN_AREA1`` — West TN Lowlands — DRNAREA + I2
``TN_AREA2`` — Middle TN (Highland Rim / Nashville Basin) — DRNAREA + CSL1085LFP
``TN_AREA3`` — Upper Cumberland / East TN Valley — DRNAREA only
``TN_AREA4`` — Ridge and Valley / Blue Ridge — DRNAREA + IMPERV

Usage
-----
::

    from hydrolib.regression.states.tennessee import (
        TN_AREA1, TN_AREA2, TN_AREA3, TN_AREA4, TN_REGIONS,
        build_tennessee_table,
    )
    from hydrolib.regression.basin_chars import BasinCharacteristics, DRNAREA, CSL1085LFP

    table = build_tennessee_table()
    basin = BasinCharacteristics(
        site_no="03606500",
        site_name="Big Sandy River at Bruceton, TN",
        region=TN_AREA2,
        predictors={DRNAREA: 779.0, CSL1085LFP: 2.4},
    )
    q100 = table.estimate(basin, aep=0.01)

Note
----
Coefficient values embedded here are representative for illustration and
testing.  Enter the authoritative values from SIR 2024-5130 Table 4 when
using for engineering design.
"""

from __future__ import annotations

from typing import Dict, Tuple

from hydrolib.regression.region import HydrologicRegion
from hydrolib.regression.regression_table import STANDARD_AEPS, GlsEquation, RegressionTable

_PUB = "SIR 2024-5130; https://doi.org/10.3133/sir20245130"

# ---------------------------------------------------------------------------
# Region definitions
# ---------------------------------------------------------------------------

TN_AREA1 = HydrologicRegion(
    code="TN-1",
    label="Tennessee Area 1 — West TN Lowlands",
    state="TN",
    required_predictors=("DRNAREA", "I2"),
    publication=_PUB,
    notes="I2 = 2-year 24-hr precipitation intensity (in/hr). "
    "Non-urban, unregulated streams in west Tennessee lowlands.",
)

TN_AREA2 = HydrologicRegion(
    code="TN-2",
    label="Tennessee Area 2 — Middle TN (Highland Rim / Nashville Basin)",
    state="TN",
    required_predictors=("DRNAREA", "CSL1085LFP"),
    publication=_PUB,
    notes="CSL1085LFP = 10-85 channel slope (ft/mi). "
    "Non-urban, unregulated streams in middle Tennessee.",
)

TN_AREA3 = HydrologicRegion(
    code="TN-3",
    label="Tennessee Area 3 — Upper Cumberland / East TN Valley",
    state="TN",
    required_predictors=("DRNAREA",),
    publication=_PUB,
    notes="Drainage area is the sole predictor. "
    "Non-urban, unregulated streams in upper Cumberland and east TN valley.",
)

TN_AREA4 = HydrologicRegion(
    code="TN-4",
    label="Tennessee Area 4 — Ridge and Valley / Blue Ridge",
    state="TN",
    required_predictors=("DRNAREA", "IMPERV"),
    publication=_PUB,
    notes="IMPERV = percent impervious cover (%). "
    "Non-urban streams in Ridge-Valley and Blue Ridge physiographic provinces.",
)

TN_REGIONS: Tuple[HydrologicRegion, ...] = (TN_AREA1, TN_AREA2, TN_AREA3, TN_AREA4)

# ---------------------------------------------------------------------------
# Regression coefficients (illustrative — fill from SIR 2024-5130 Table 4)
# ---------------------------------------------------------------------------
# Format: {aep: (intercept, {code: exponent}, sep_pct, pseudo_r2, eyr)}

_TN1_COEFFS: Dict[float, tuple] = {
    0.500: (0.25, {"DRNAREA": 0.73, "I2": 0.78}, 52.0, 0.86, 10.0),
    0.200: (0.36, {"DRNAREA": 0.74, "I2": 0.80}, 46.0, 0.88, 12.0),
    0.100: (0.42, {"DRNAREA": 0.74, "I2": 0.81}, 43.0, 0.89, 13.0),
    0.040: (0.49, {"DRNAREA": 0.75, "I2": 0.82}, 40.0, 0.90, 14.0),
    0.020: (0.53, {"DRNAREA": 0.75, "I2": 0.82}, 39.0, 0.90, 14.0),
    0.010: (0.57, {"DRNAREA": 0.76, "I2": 0.83}, 38.0, 0.90, 14.0),
    0.005: (0.61, {"DRNAREA": 0.76, "I2": 0.83}, 40.0, 0.89, 13.0),
    0.002: (0.67, {"DRNAREA": 0.77, "I2": 0.84}, 43.0, 0.88, 12.0),
}

_TN2_COEFFS: Dict[float, tuple] = {
    0.500: (1.60, {"DRNAREA": 0.72, "CSL1085LFP": 0.35}, 46.0, 0.89, 14.0),
    0.200: (1.96, {"DRNAREA": 0.73, "CSL1085LFP": 0.36}, 42.0, 0.90, 15.0),
    0.100: (2.17, {"DRNAREA": 0.74, "CSL1085LFP": 0.37}, 39.0, 0.91, 17.0),
    0.040: (2.40, {"DRNAREA": 0.74, "CSL1085LFP": 0.37}, 37.0, 0.92, 18.0),
    0.020: (2.55, {"DRNAREA": 0.75, "CSL1085LFP": 0.38}, 36.0, 0.92, 18.0),
    0.010: (2.69, {"DRNAREA": 0.75, "CSL1085LFP": 0.38}, 35.0, 0.92, 18.0),
    0.005: (2.82, {"DRNAREA": 0.76, "CSL1085LFP": 0.38}, 36.0, 0.91, 17.0),
    0.002: (3.00, {"DRNAREA": 0.77, "CSL1085LFP": 0.38}, 38.0, 0.90, 16.0),
}

_TN3_COEFFS: Dict[float, tuple] = {
    0.500: (1.55, {"DRNAREA": 0.75}, 54.0, 0.84, 10.0),
    0.200: (1.88, {"DRNAREA": 0.76}, 49.0, 0.86, 11.0),
    0.100: (2.09, {"DRNAREA": 0.77}, 46.0, 0.87, 12.0),
    0.040: (2.33, {"DRNAREA": 0.77}, 44.0, 0.88, 13.0),
    0.020: (2.48, {"DRNAREA": 0.78}, 43.0, 0.88, 13.0),
    0.010: (2.62, {"DRNAREA": 0.78}, 42.0, 0.88, 13.0),
    0.005: (2.75, {"DRNAREA": 0.79}, 43.0, 0.87, 12.0),
    0.002: (2.93, {"DRNAREA": 0.80}, 46.0, 0.86, 11.0),
}

_TN4_COEFFS: Dict[float, tuple] = {
    0.500: (1.45, {"DRNAREA": 0.71, "IMPERV": 0.15}, 50.0, 0.85, 11.0),
    0.200: (1.79, {"DRNAREA": 0.72, "IMPERV": 0.16}, 45.0, 0.87, 13.0),
    0.100: (2.00, {"DRNAREA": 0.73, "IMPERV": 0.17}, 43.0, 0.88, 14.0),
    0.040: (2.23, {"DRNAREA": 0.73, "IMPERV": 0.17}, 41.0, 0.89, 15.0),
    0.020: (2.37, {"DRNAREA": 0.74, "IMPERV": 0.18}, 40.0, 0.89, 15.0),
    0.010: (2.50, {"DRNAREA": 0.74, "IMPERV": 0.18}, 40.0, 0.89, 15.0),
    0.005: (2.63, {"DRNAREA": 0.75, "IMPERV": 0.18}, 42.0, 0.88, 14.0),
    0.002: (2.80, {"DRNAREA": 0.75, "IMPERV": 0.19}, 45.0, 0.87, 13.0),
}

_AREA_COEFFS = {
    TN_AREA1: _TN1_COEFFS,
    TN_AREA2: _TN2_COEFFS,
    TN_AREA3: _TN3_COEFFS,
    TN_AREA4: _TN4_COEFFS,
}


def build_tennessee_table() -> RegressionTable:
    """
    Build a :class:`RegressionTable` for all four Tennessee hydrologic areas.

    Returns a table with 32 equations (4 areas × 8 standard AEPs).

    Returns
    -------
    RegressionTable
        Publication set to SIR 2024-5130.

    Notes
    -----
    Coefficient values are representative for illustration.  Replace with
    values from Table 4 of SIR 2024-5130 for engineering design use.
    """
    table = RegressionTable(publication=_PUB)
    for region in TN_REGIONS:
        coeffs = _AREA_COEFFS[region]
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
