"""
hydrolib.regional_skew - Lookup of USGS state/regional generalized skew studies

USGS regional (generalized) skew is *not* a single nationwide constant looked
up per state. It is produced by separate, periodically-revised state or
multi-state regression studies (frequently Bayesian GLS), each published as
its own USGS Scientific Investigations Report (SIR); some states instead
have a confirmed recommendation to keep using the nationwide Bulletin 17B
map (Interagency Advisory Committee on Water Data, 1982, Plate 1). USGS
maintains the authoritative, continuously-updated index of the current study
for every state at:

    https://www.usgs.gov/streamstats/science/flood-frequency-reports

This module ships a small, explicitly-sourced starter dataset covering a
handful of states, plus the nationwide Bulletin 17B fallback used when a
state has no dedicated entry here. It is **not** a complete substitute for
the USGS index above -- consult the current report for your study area
before using any value from this module in an engineering analysis.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

# Nationwide fallback: Interagency Advisory Committee on Water Data (1982),
# Bulletin 17B, Guidelines for Determining Flood Flow Frequency, Plate 1.
NATIONWIDE_SKEW = -0.302
NATIONWIDE_SKEW_SE = 0.55


@dataclass(frozen=True)
class RegionalSkewEstimate:
    """A generalized (regional) skew estimate for a state or region.

    Attributes
    ----------
    value : float
        Generalized skew coefficient.
    se : float
        Standard error of the generalized skew estimate.
    source : str
        Citation for the study this value comes from.
    source_url : str
        URL of the source report.
    states : List[str]
        USPS state abbreviations this estimate applies to.
    note : Optional[str]
        Additional context, e.g. when a state study recommends falling back
        to the nationwide Bulletin 17B value rather than a state-specific
        constant.
    """

    value: float
    se: float
    source: str
    source_url: str
    states: List[str]
    note: Optional[str] = None

    @property
    def mse(self) -> float:
        """Mean square error, computed as ``se ** 2``."""
        return self.se**2


NATIONWIDE_FALLBACK = RegionalSkewEstimate(
    value=NATIONWIDE_SKEW,
    se=NATIONWIDE_SKEW_SE,
    source=(
        "Interagency Advisory Committee on Water Data, 1982, Bulletin 17B, "
        "Guidelines for Determining Flood Flow Frequency, Plate 1 "
        "(nationwide generalized skew map)"
    ),
    source_url="https://water.usgs.gov/osw/bulletin17b/dl_flow.pdf",
    states=[],
    note="Nationwide default used when no state-specific study is available.",
)

# Starter dataset of state regional-skew studies, verified against published
# USGS reports. Each state not listed here falls back to NATIONWIDE_FALLBACK.
# To add a state, cite the current report from the USGS Flood Frequency
# Reports index (see module docstring) -- do not guess or estimate values.
_STUDIES: List[RegionalSkewEstimate] = [
    RegionalSkewEstimate(
        value=0.44,
        se=0.078**0.5,
        source=(
            "Olson, S.A., 2014, Estimating flood magnitude and frequency at "
            "gaged and ungaged sites on streams in Vermont, based on data "
            "through water year 2011: U.S. Geological Survey Scientific "
            "Investigations Report 2014-5078"
        ),
        source_url="https://pubs.usgs.gov/sir/2014/5078/",
        states=["VT"],
        note=(
            "Constant regional skew developed by generalized least-squares "
            "regression using 145 streamgages. A newer Vermont report "
            "(SIR 2025-5088) may supersede this value -- verify against the "
            "current USGS Flood Frequency Reports index before use."
        ),
    ),
    RegionalSkewEstimate(
        value=0.048,
        se=0.092**0.5,
        source=(
            "Feaster, T.D., and others, 2023, Magnitude and frequency of "
            "floods for rural streams in Georgia, South Carolina, and North "
            "Carolina, 2017: U.S. Geological Survey Scientific "
            "Investigations Report 2023-5006"
        ),
        source_url="https://pubs.usgs.gov/publication/sir20235006",
        states=["GA", "SC", "NC"],
        note="Constant regional skew developed by Bayesian GLS regression.",
    ),
    RegionalSkewEstimate(
        value=NATIONWIDE_SKEW,
        se=NATIONWIDE_SKEW_SE,
        source=(
            "Sando, S.K., 1998, Techniques for estimating peak-flow "
            "magnitude and frequency relations for South Dakota streams: "
            "U.S. Geological Survey Water-Resources Investigations Report "
            "98-4055"
        ),
        source_url="https://pubs.usgs.gov/wri/wri98-4055/",
        states=["SD"],
        note=(
            "This state study tested the Bulletin 17B generalized skew map "
            "against South Dakota gage data and found it adequate; it "
            "recommends using the nationwide value rather than a "
            "state-specific constant."
        ),
    ),
    RegionalSkewEstimate(
        value=NATIONWIDE_SKEW,
        se=0.64,
        source=(
            "Parrett, C., and Johnson, D.R., 2004, Methods for estimating "
            "flood frequency in Montana based on data through water year "
            "1998: U.S. Geological Survey Water-Resources Investigations "
            "Report 03-4308"
        ),
        source_url="https://pubs.usgs.gov/wri/wri03-4308",
        states=["MT"],
        note=(
            "This study found the difference between Montana-specific and "
            "Bulletin 17B nationwide generalized skew to be small and "
            "probably not significant, so it kept the nationwide skew "
            "value but recalibrated its standard error for Montana gages "
            "to 0.64 (versus the published nationwide SE of 0.55) -- that "
            "recalibrated SE is used here. Newer studies (SIR 2018-5046; "
            "SIR 2025-5019, covering MT/ND/SD/WY) may supersede this -- "
            "verify against the current USGS Flood Frequency Reports "
            "index before use."
        ),
    ),
    RegionalSkewEstimate(
        value=NATIONWIDE_SKEW,
        se=NATIONWIDE_SKEW_SE,
        source=(
            "U.S. Geological Survey, 2016, Magnitude, frequency, and "
            "trends of floods at gaged and ungaged sites in Washington, "
            "based on data through water year 2014 (ver. 1.1, October "
            "2016): U.S. Geological Survey Scientific Investigations "
            "Report 2016-5118"
        ),
        source_url="https://pubs.usgs.gov/sir/2016/5118/",
        states=["WA"],
        note=(
            "No dedicated regional-skew study has been completed for "
            "Washington; this report uses the nationwide generalized skew "
            "from the U.S. Water Resources Council (1981) / Bulletin 17B. "
            "Part of the state is also covered by SIR 2020-5073 (Columbia "
            "River Basin regional skew models by flood duration), which "
            "vary spatially and by duration rather than being a single "
            "constant, so are not represented in this module -- verify "
            "against the current USGS Flood Frequency Reports index "
            "before use."
        ),
    ),
    RegionalSkewEstimate(
        value=-0.09,
        se=0.08**0.5,
        source=(
            "Paretti, N.V., and others, 2014, Methods for estimating "
            "magnitude and frequency of floods in Arizona, developed with "
            "unregulated and rural peak-flow data through water year 2010: "
            "U.S. Geological Survey Scientific Investigations Report "
            "2014-5211"
        ),
        source_url="https://pubs.usgs.gov/sir/2014/5211/",
        states=["AZ"],
        note=(
            "Statewide constant regional skew from a Bayesian "
            "generalized least-squares analysis of 448 gages -- no basin "
            "characteristic explained the variation in skew well enough "
            "to justify a spatially-varying model, so a constant was "
            "adopted (MSE 0.08, versus MSE 0.31 for the prior Arizona "
            "regional skew analysis)."
        ),
    ),
    RegionalSkewEstimate(
        value=-0.30,
        se=0.14**0.5,
        source=(
            "U.S. Geological Survey, 2014, Methods for estimating annual "
            "exceedance-probability discharges and largest recorded "
            "floods for unregulated streams in rural Missouri: U.S. "
            "Geological Survey Scientific Investigations Report 2014-5165"
        ),
        source_url="https://pubs.usgs.gov/sir/2014/5165/",
        states=["MO"],
        note=(
            "Statewide constant regional skew from a Bayesian weighted "
            "least-squares/generalized least-squares analysis of 108 "
            "long-term gages (30+ years of record); 35 basin "
            "characteristics were tested and none improved on the "
            "constant model (MSE 0.14)."
        ),
    ),
    RegionalSkewEstimate(
        value=0.029,
        se=0.297,
        source=(
            "U.S. Geological Survey, 1999, Estimating the Magnitude of "
            "Peak Flows for Streams in Maine for Selected Recurrence "
            "Intervals: U.S. Geological Survey Water-Resources "
            "Investigations Report 99-4008"
        ),
        source_url="https://pubs.usgs.gov/publication/wri994008",
        states=["ME"],
        note=(
            "Statewide constant regional skew with standard error of "
            "prediction 0.297, given directly (not as an MSE) in the "
            "source report."
        ),
    ),
    RegionalSkewEstimate(
        value=0.37,
        se=0.14**0.5,
        source=(
            "U.S. Geological Survey, 2017, Methods for estimating "
            "regional coefficient of skewness for unregulated streams in "
            "New England: U.S. Geological Survey Scientific "
            "Investigations Report 2017-5037"
        ),
        source_url="https://pubs.usgs.gov/publication/sir20175037",
        states=["CT", "RI", "MA"],
        note=(
            "Constant regional skew for a Connecticut/Rhode "
            "Island/Massachusetts study area (model error variance 0.13, "
            "average variance of prediction at a new site 0.14, used here "
            "as the MSE). Maine and New Hampshire have their own separate, "
            "older, spatially-varying skew maps and are not part of this "
            "constant."
        ),
    ),
    RegionalSkewEstimate(
        value=0.50,
        se=0.574,
        source=(
            "U.S. Geological Survey, 2025, Estimation of magnitude and "
            "frequency of floods for rural, unregulated streams in and "
            "near Virginia and West Virginia: U.S. Geological Survey "
            "Scientific Investigations Report 2025-5110"
        ),
        source_url="https://pubs.usgs.gov/publication/sir20255110",
        states=["VA", "WV"],
        note=(
            "Constant regional skew for a study area covering parts of "
            "the Mid-Atlantic and South-Atlantic-Gulf regions (average "
            "variance of prediction 0.33, standard error 0.574 given "
            "directly in the source); none of the explanatory variables "
            "tested had enough predictive power to justify a "
            "spatially-varying model. Part of western Maryland falls in "
            "the same study region -- see the Maryland/Delaware entry in "
            "this module's documentation for a separate, more localized "
            "Eastern Coastal Plain estimate. Kentucky and Tennessee, "
            "though grouped with Virginia and West Virginia in some "
            "earlier USGS regional-skew work (Feaster and others, 2023), "
            "are not confirmed to share this specific constant and are "
            "not included here."
        ),
    ),
    RegionalSkewEstimate(
        value=NATIONWIDE_SKEW,
        se=NATIONWIDE_SKEW_SE,
        source=(
            "U.S. Geological Survey, 2015, Regional regression equations "
            "to estimate peak-flow frequency at sites in North Dakota "
            "using data through 2009: U.S. Geological Survey Scientific "
            "Investigations Report 2015-5096"
        ),
        source_url="https://pubs.usgs.gov/publication/sir20155096",
        states=["ND"],
        note=(
            "This state study found the Bulletin 17B nationwide "
            "generalized skew map provides accurate estimates for "
            "natural-flow streams in North Dakota and recommends using "
            "the nationwide value rather than a state-specific constant "
            "(same conclusion as South Dakota's study)."
        ),
    ),
]

# Every other state was investigated but deliberately NOT included above,
# because either no single citable numeric constant could be found, or the
# state's own study explicitly uses a spatially-varying skew (a map, a set
# of regional zones, or a function of elevation/basin characteristics)
# rather than one constant, so it doesn't fit this module's (value, se)
# schema. Adding a single fabricated number for any of them would
# misrepresent the source study, so they fall back to NATIONWIDE_FALLBACK.
#
# A full state-by-state research table -- including states with a
# spatially-varying model, states where a value was found but not a
# matching standard error/MSE (e.g. Arkansas: -0.17; Hawaii: -0.14), and
# states where nothing citable could be confirmed at all -- is maintained
# outside this module (see the project's regional-skew research
# deliverable). Consult the cited reports there, or the live USGS Flood
# Frequency Reports index (https://www.usgs.gov/streamstats/science/
# flood-frequency-reports), directly for a site-specific value.
#
# Representative examples of *why* a state is excluded rather than just
# absent:
#   Idaho, Oregon, Wyoming, California -- each uses several spatially
#       distinct skew zones, or skew as an explicit function of basin
#       elevation, rather than one statewide constant.
#   New Mexico, Texas -- generalized skew is a spatially-varying
#       GAM-based surface (shared multi-state report covering Texas,
#       Oklahoma, and only the eastern part of New Mexico); western New
#       Mexico is not covered by that report at all.
#   New York, Pennsylvania, Maryland, Delaware -- a confirmed constant
#       exists (e.g. 0.32 for eastern New York/Pennsylvania), but only for
#       part of the state, so it is not representative of the whole state.
#   Nevada, Colorado, Utah, New Jersey, Louisiana, Puerto Rico, and others
#       -- a dedicated USGS report exists, but no single statewide
#       constant (and/or its standard error) could be confirmed from
#       available sources in this research pass.

# USPS state/territory abbreviation -> full name, for input normalization
# and for building UI pickers.
STATE_NAMES: Dict[str, str] = {
    "AL": "Alabama",
    "AK": "Alaska",
    "AZ": "Arizona",
    "AR": "Arkansas",
    "CA": "California",
    "CO": "Colorado",
    "CT": "Connecticut",
    "DE": "Delaware",
    "DC": "District of Columbia",
    "FL": "Florida",
    "GA": "Georgia",
    "HI": "Hawaii",
    "ID": "Idaho",
    "IL": "Illinois",
    "IN": "Indiana",
    "IA": "Iowa",
    "KS": "Kansas",
    "KY": "Kentucky",
    "LA": "Louisiana",
    "ME": "Maine",
    "MD": "Maryland",
    "MA": "Massachusetts",
    "MI": "Michigan",
    "MN": "Minnesota",
    "MS": "Mississippi",
    "MO": "Missouri",
    "MT": "Montana",
    "NE": "Nebraska",
    "NV": "Nevada",
    "NH": "New Hampshire",
    "NJ": "New Jersey",
    "NM": "New Mexico",
    "NY": "New York",
    "NC": "North Carolina",
    "ND": "North Dakota",
    "OH": "Ohio",
    "OK": "Oklahoma",
    "OR": "Oregon",
    "PA": "Pennsylvania",
    "PR": "Puerto Rico",
    "RI": "Rhode Island",
    "SC": "South Carolina",
    "SD": "South Dakota",
    "TN": "Tennessee",
    "TX": "Texas",
    "UT": "Utah",
    "VT": "Vermont",
    "VA": "Virginia",
    "WA": "Washington",
    "WV": "West Virginia",
    "WI": "Wisconsin",
    "WY": "Wyoming",
}
_NAME_TO_CODE: Dict[str, str] = {name.upper(): code for code, name in STATE_NAMES.items()}

STATE_SKEW: Dict[str, RegionalSkewEstimate] = {
    state: study for study in _STUDIES for state in study.states
}


def normalize_state(state: str) -> str:
    """Normalize a state name or USPS abbreviation to its two-letter code.

    Parameters
    ----------
    state : str
        Full state/territory name (e.g. ``"Vermont"``) or USPS abbreviation
        (e.g. ``"VT"``), case-insensitive.

    Returns
    -------
    str
        The USPS two-letter code.

    Raises
    ------
    ValueError
        If ``state`` is not a recognized U.S. state, DC, or Puerto Rico.
    """
    key = state.strip().upper()
    if key in STATE_NAMES:
        return key
    if key in _NAME_TO_CODE:
        return _NAME_TO_CODE[key]
    raise ValueError(f"Unrecognized state or territory: {state!r}")


def get_regional_skew(state: str) -> RegionalSkewEstimate:
    """Look up the generalized (regional) skew estimate for a state.

    Returns the state-specific study if one is available in this module's
    starter dataset, otherwise falls back to the nationwide Bulletin 17B
    value. This dataset covers only a handful of states -- check
    ``estimate.states`` and ``estimate.note`` on the result, and consult the
    USGS Flood Frequency Reports index (see module docstring) for the
    current, authoritative study before use in an engineering analysis.

    Parameters
    ----------
    state : str
        Full state/territory name or USPS abbreviation, case-insensitive.

    Returns
    -------
    RegionalSkewEstimate
        The matching state study, or ``NATIONWIDE_FALLBACK`` if none exists.
    """
    code = normalize_state(state)
    estimate = STATE_SKEW.get(code)
    if estimate is None:
        logger.info(
            "No state-specific regional skew study for %s in the starter "
            "dataset; using nationwide Bulletin 17B fallback (%.3f). Check "
            "https://www.usgs.gov/streamstats/science/flood-frequency-reports "
            "for the current study.",
            code,
            NATIONWIDE_SKEW,
        )
        return NATIONWIDE_FALLBACK
    return estimate


def list_available_states() -> List[str]:
    """Return USPS codes for states with a dedicated entry in the starter dataset."""
    return sorted(STATE_SKEW.keys())
