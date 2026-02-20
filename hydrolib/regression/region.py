"""
hydrolib.regression.region - Generic hydrologic region definition.

A :class:`HydrologicRegion` instance describes one region (or area/zone)
for which a separate set of GLS regression equations was developed.  Regions
can be loaded from a CSV, constructed programmatically, or obtained from
publication-specific helpers such as :mod:`hydrolib.regression.sir2024_5130`.

StreamStats predictor codes
---------------------------
The ``required_predictors`` field uses **USGS StreamStats parameter codes**
as canonical keys so that predictor names are consistent across states and
publications.  Common codes:

============  =============================================
Code          Description
============  =============================================
DRNAREA       Drainage area (sq mi)
CSL1085LFP    10-85 channel slope (ft/mi)
I2            2-yr precipitation intensity (in/hr or in)
IMPERV        Impervious cover (%)
ELEV          Mean basin elevation (ft)
FOREST        Forested area (%)
PRECIP        Mean annual precipitation (in)
BFIPCT        Base-flow index (%)
SLOPE         Mean basin slope (%)
============  =============================================

A full list is available at https://streamstats.usgs.gov.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Tuple

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class HydrologicRegion:
    """
    Defines one hydrologic region for which GLS regression equations exist.

    Instances are **hashable** and can be used as dictionary keys.

    Parameters
    ----------
    code : str
        Short machine-readable identifier, e.g. ``"TN_AREA2"`` or
        ``"KY_REGION4"``.  Must be unique within a publication.
    label : str
        Human-readable description, e.g. ``"Tennessee Hydrologic Area 2"``.
    state : str
        Two-letter state abbreviation (or comma-separated list for
        multi-state studies), e.g. ``"TN"`` or ``"TN,AL,GA"``.
    required_predictors : tuple of str
        StreamStats parameter codes that **must** be present in a basin's
        predictor dict before any regression equation can be applied.
        Order does not matter.
    publication : str, optional
        Citation or DOI for the source report, e.g.
        ``"SIR 2024-5130; https://doi.org/10.3133/sir20245130"``.
    notes : str, optional
        Free-text notes (applicability limits, special considerations, etc.).

    Examples
    --------
    >>> tn_area2 = HydrologicRegion(
    ...     code="TN_AREA2",
    ...     label="Tennessee Hydrologic Area 2 — Middle TN",
    ...     state="TN",
    ...     required_predictors=("DRNAREA", "CSL1085LFP"),
    ...     publication="SIR 2024-5130; https://doi.org/10.3133/sir20245130",
    ... )
    >>> ky_reg1 = HydrologicRegion(
    ...     code="KY_REGION1",
    ...     label="Kentucky Flood Region 1",
    ...     state="KY",
    ...     required_predictors=("DRNAREA", "CSL1085LFP"),
    ...     publication="WRIR 03-4180",
    ... )
    """

    code: str
    label: str
    state: str
    required_predictors: Tuple[str, ...] = ()
    publication: str = ""
    notes: str = ""

    # ------------------------------------------------------------------
    # Validation helpers
    # ------------------------------------------------------------------

    def validate_predictors(self, predictors: Dict[str, float]) -> None:
        """
        Raise ``ValueError`` if any required predictor is absent or non-positive.

        Parameters
        ----------
        predictors : dict
            Mapping StreamStats code -> value from a :class:`BasinCharacteristics`.

        Raises
        ------
        ValueError
            If a required predictor is missing or has a non-positive value.
        """
        missing = []
        invalid = []
        for code in self.required_predictors:
            val = predictors.get(code)
            if val is None:
                missing.append(code)
            elif val <= 0:
                invalid.append(f"{code}={val}")
        errors = []
        if missing:
            errors.append(f"missing predictors: {missing}")
        if invalid:
            errors.append(f"non-positive values: {invalid}")
        if errors:
            raise ValueError(f"Region {self.code!r} validation failed — {'; '.join(errors)}")

    def has_predictor(self, code: str) -> bool:
        """Return True if *code* is in ``required_predictors``."""
        return code in self.required_predictors

    # ------------------------------------------------------------------
    # I/O helpers
    # ------------------------------------------------------------------

    @classmethod
    def from_dict(cls, row: Dict[str, str]) -> "HydrologicRegion":
        """
        Construct from a flat dict (e.g., one row of a CSV).

        Expected keys: ``code``, ``label``, ``state``.
        Optional keys: ``required_predictors`` (comma-separated codes),
        ``publication``, ``notes``.
        """
        preds_raw = row.get("required_predictors", "").strip()
        preds: Tuple[str, ...] = (
            tuple(p.strip() for p in preds_raw.split(",") if p.strip()) if preds_raw else ()
        )
        return cls(
            code=row["code"].strip(),
            label=row.get("label", row["code"]).strip(),
            state=row.get("state", "").strip(),
            required_predictors=preds,
            publication=row.get("publication", "").strip(),
            notes=row.get("notes", "").strip(),
        )

    def to_dict(self) -> Dict[str, str]:
        """Serialise to a flat dict (inverse of :meth:`from_dict`)."""
        return {
            "code": self.code,
            "label": self.label,
            "state": self.state,
            "required_predictors": ", ".join(self.required_predictors),
            "publication": self.publication,
            "notes": self.notes,
        }

    def __str__(self) -> str:
        return f"{self.code} ({self.label}, {self.state})"


# ---------------------------------------------------------------------------
# Convenience registry
# ---------------------------------------------------------------------------


class RegionRegistry:
    """
    Simple in-process registry mapping region codes to :class:`HydrologicRegion`.

    Useful when loading equation tables from CSV so that identical region
    objects are shared rather than duplicated.

    Examples
    --------
    >>> registry = RegionRegistry()
    >>> registry.register(HydrologicRegion("TN_AREA2", "TN Area 2", "TN",
    ...                                    required_predictors=("DRNAREA", "CSL1085LFP")))
    >>> region = registry["TN_AREA2"]
    """

    def __init__(self) -> None:
        self._store: Dict[str, HydrologicRegion] = {}

    def register(self, region: HydrologicRegion) -> None:
        """Add *region* to the registry (overwrites if code exists)."""
        self._store[region.code] = region

    def get(self, code: str) -> Optional[HydrologicRegion]:
        """Return region by code, or None."""
        return self._store.get(code)

    def __getitem__(self, code: str) -> HydrologicRegion:
        try:
            return self._store[code]
        except KeyError:
            raise KeyError(f"Region {code!r} not in registry.  Available: {list(self._store)}")

    def __contains__(self, code: str) -> bool:
        return code in self._store

    def __len__(self) -> int:
        return len(self._store)

    def all_regions(self) -> List[HydrologicRegion]:
        return list(self._store.values())

    @classmethod
    def load_from_csv(cls, path: str | Path) -> "RegionRegistry":
        """
        Build a registry from a CSV file.

        CSV must have at minimum a ``code`` column.  Optional columns:
        ``label``, ``state``, ``required_predictors``, ``publication``, ``notes``.

        Parameters
        ----------
        path : str or Path

        Returns
        -------
        RegionRegistry
        """
        path = Path(path)
        registry = cls()
        with path.open(newline="") as fh:
            for row in csv.DictReader(fh):
                region = HydrologicRegion.from_dict(row)
                registry.register(region)
        log.info("Loaded %d regions from %s", len(registry), path)
        return registry
