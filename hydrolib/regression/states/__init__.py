"""
hydrolib.regression.states - Per-state GLS regression equation modules.

Each sub-module defines the hydrologic regions and coefficient tables for
one state, following a consistent pattern::

    <STATE>_<REGION> : HydrologicRegion   — region definition (hashable)
    <STATE>_REGIONS  : tuple              — all regions for the state
    build_<state>_table() -> RegressionTable  — construct the full table

Supported states
----------------
* ``alabama``    — Appalachian Plateaus, Valley/Ridge, Piedmont, Coastal Plain
* ``california`` — North Coast, Sierra Nevada, Central Valley, Southern
* ``colorado``   — Front Range, Eastern Plains, Western Slope, San Luis Valley
* ``georgia``    — Blue Ridge, Valley/Ridge, Piedmont, Coastal Plain
* ``idaho``      — North, Central Mountains, Snake Plain, Southern Highlands
* ``montana``    — Mountain (DA + ELEV + PRECIP), Foothills, Plains
* ``new_mexico`` — Mountain, Transitional, Plains, Desert (4 regions)
* ``oregon``     — Coast, Willamette Valley, Cascades, Eastern
* ``tennessee``  — SIR 2024-5130 (4 areas, DA + slope/I2/IMPERV)
* ``washington`` — Olympic/Coast, W. Cascades, E. Cascades, Eastern
* ``wyoming``    — Mountain, Foothills, Eastern Plains (3 regions)

Quick start — nationwide merged table
--------------------------------------
::

    from hydrolib.regression.regression_table import RegressionTable
    from hydrolib.regression.states import (
        build_alabama_table, build_tennessee_table, build_georgia_table,
        build_montana_table, build_new_mexico_table, build_colorado_table,
        build_wyoming_table, build_idaho_table, build_washington_table,
        build_oregon_table, build_california_table,
    )

    national = RegressionTable.merge(
        build_alabama_table(),
        build_tennessee_table(),
        build_georgia_table(),
        build_montana_table(),
        build_new_mexico_table(),
        build_colorado_table(),
        build_wyoming_table(),
        build_idaho_table(),
        build_washington_table(),
        build_oregon_table(),
        build_california_table(),
    )
    # national.available_states()
    # → ['AL', 'CA', 'CO', 'GA', 'ID', 'MT', 'NM', 'OR', 'TN', 'WA', 'WY']
    # national.batch_estimate(site_list, aep=0.01)
"""

from hydrolib.regression.states.alabama import (
    AL_APPALACHIAN,
    AL_COASTAL,
    AL_PIEDMONT,
    AL_REGIONS,
    AL_VALLEY_RIDGE,
    build_alabama_table,
)
from hydrolib.regression.states.california import (
    CA_CENTRAL,
    CA_NORTH_COAST,
    CA_REGIONS,
    CA_SIERRA,
    CA_SOUTH,
    build_california_table,
)
from hydrolib.regression.states.colorado import (
    CO_EASTERN,
    CO_FRONT_RANGE,
    CO_REGIONS,
    CO_SAN_LUIS,
    CO_WESTERN,
    build_colorado_table,
)
from hydrolib.regression.states.georgia import (
    GA_BLUE_RIDGE,
    GA_COASTAL,
    GA_PIEDMONT,
    GA_REGIONS,
    GA_VALLEY_RIDGE,
    build_georgia_table,
)
from hydrolib.regression.states.idaho import (
    ID_CENTRAL,
    ID_NORTH,
    ID_REGIONS,
    ID_SNAKE,
    ID_SOUTH,
    build_idaho_table,
)
from hydrolib.regression.states.montana import (
    MT_FOOTHILLS,
    MT_MOUNTAIN,
    MT_PLAINS,
    MT_REGIONS,
    build_montana_table,
)
from hydrolib.regression.states.new_mexico import (
    NM_DESERT,
    NM_MOUNTAIN,
    NM_PLAINS,
    NM_REGIONS,
    NM_TRANSITION,
    build_new_mexico_table,
)
from hydrolib.regression.states.oregon import (
    OR_CASCADE,
    OR_COAST,
    OR_EASTERN,
    OR_REGIONS,
    OR_WILLAMETTE,
    build_oregon_table,
)
from hydrolib.regression.states.tennessee import (
    TN_AREA1,
    TN_AREA2,
    TN_AREA3,
    TN_AREA4,
    TN_REGIONS,
    build_tennessee_table,
)
from hydrolib.regression.states.washington import (
    WA_E_CASCADE,
    WA_EASTERN,
    WA_OLYMPIC,
    WA_REGIONS,
    WA_W_CASCADE,
    build_washington_table,
)
from hydrolib.regression.states.wyoming import (
    WY_FOOTHILLS,
    WY_MOUNTAIN,
    WY_PLAINS,
    WY_REGIONS,
    build_wyoming_table,
)

__all__ = [
    # Alabama
    "AL_APPALACHIAN",
    "AL_VALLEY_RIDGE",
    "AL_PIEDMONT",
    "AL_COASTAL",
    "AL_REGIONS",
    "build_alabama_table",
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
    # New Mexico
    "NM_MOUNTAIN",
    "NM_TRANSITION",
    "NM_PLAINS",
    "NM_DESERT",
    "NM_REGIONS",
    "build_new_mexico_table",
    # Colorado
    "CO_FRONT_RANGE",
    "CO_EASTERN",
    "CO_WESTERN",
    "CO_SAN_LUIS",
    "CO_REGIONS",
    "build_colorado_table",
    # Wyoming
    "WY_MOUNTAIN",
    "WY_FOOTHILLS",
    "WY_PLAINS",
    "WY_REGIONS",
    "build_wyoming_table",
    # Idaho
    "ID_NORTH",
    "ID_CENTRAL",
    "ID_SNAKE",
    "ID_SOUTH",
    "ID_REGIONS",
    "build_idaho_table",
    # Washington
    "WA_OLYMPIC",
    "WA_W_CASCADE",
    "WA_E_CASCADE",
    "WA_EASTERN",
    "WA_REGIONS",
    "build_washington_table",
    # Oregon
    "OR_COAST",
    "OR_WILLAMETTE",
    "OR_CASCADE",
    "OR_EASTERN",
    "OR_REGIONS",
    "build_oregon_table",
    # California
    "CA_NORTH_COAST",
    "CA_SIERRA",
    "CA_CENTRAL",
    "CA_SOUTH",
    "CA_REGIONS",
    "build_california_table",
]
