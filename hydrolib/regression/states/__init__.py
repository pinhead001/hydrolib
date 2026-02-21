"""
hydrolib.regression.states - Per-state GLS regression equation modules.

Each sub-module defines the hydrologic regions and coefficient tables for
one state, following a consistent pattern::

    <STATE>_<REGION> : HydrologicRegion   — region definition (hashable)
    <STATE>_REGIONS  : tuple              — all regions for the state
    build_<state>_table() -> RegressionTable  — construct the full table

Supported states
----------------
* ``tennessee`` — SIR 2024-5130 (4 areas, DA + slope/I2/IMPERV)
* ``georgia``   — Blue Ridge, Valley/Ridge, Piedmont, Coastal Plain
* ``montana``   — Mountain (DA + ELEV + PRECIP), Foothills, Plains

Quick start — nationwide merged table
--------------------------------------
::

    from hydrolib.regression.regression_table import RegressionTable
    from hydrolib.regression.states.tennessee import build_tennessee_table
    from hydrolib.regression.states.georgia import build_georgia_table
    from hydrolib.regression.states.montana import build_montana_table

    national = RegressionTable.merge(
        build_tennessee_table(),
        build_georgia_table(),
        build_montana_table(),
    )
    # national.available_states() → ['GA', 'MT', 'TN']
    # national.batch_estimate(site_list, aep=0.01)
"""

from hydrolib.regression.states.georgia import (
    GA_BLUE_RIDGE,
    GA_COASTAL,
    GA_PIEDMONT,
    GA_REGIONS,
    GA_VALLEY_RIDGE,
    build_georgia_table,
)
from hydrolib.regression.states.montana import (
    MT_FOOTHILLS,
    MT_MOUNTAIN,
    MT_PLAINS,
    MT_REGIONS,
    build_montana_table,
)
from hydrolib.regression.states.tennessee import (
    TN_AREA1,
    TN_AREA2,
    TN_AREA3,
    TN_AREA4,
    TN_REGIONS,
    build_tennessee_table,
)

__all__ = [
    # Tennessee
    "TN_AREA1",
    "TN_AREA2",
    "TN_AREA3",
    "TN_AREA4",
    "TN_REGIONS",
    "build_tennessee_table",
    # Georgia
    "GA_BLUE_RIDGE",
    "GA_VALLEY_RIDGE",
    "GA_PIEDMONT",
    "GA_COASTAL",
    "GA_REGIONS",
    "build_georgia_table",
    # Montana
    "MT_MOUNTAIN",
    "MT_FOOTHILLS",
    "MT_PLAINS",
    "MT_REGIONS",
    "build_montana_table",
]
