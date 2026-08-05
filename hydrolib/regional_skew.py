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
]

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
