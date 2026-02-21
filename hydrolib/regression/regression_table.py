"""
hydrolib.regression.regression_table - Generic GLS regression equation table.

:class:`RegressionTable` is a publication-agnostic container for GLS
regression equations of the form::

    log10(Q_T) = intercept + sum_j( coeff_j * log10( predictor_j ) )

where each equation is indexed by a (:class:`HydrologicRegion`, AEP) pair
and the predictor set is arbitrary — enabling any state and any combination
of basin characteristics.

CSV Format
----------
The canonical on-disk format is a "wide" CSV where predictor coefficients
occupy named columns::

    region_code, state, aep, intercept, DRNAREA, CSL1085LFP, ELEV, PRECIP, ...,
    sep_pct, pseudo_r2, eyr

Rules:

* Any column whose name is **not** in the reserved set
  ``{region_code, region_label, state, aep, intercept, sep_pct, pseudo_r2,
  eyr, publication, return_period_yr, notes}`` is interpreted as a predictor
  coefficient column.
* A blank or NaN cell in a predictor column means that predictor is **not
  used** in that equation row.
* A companion ``region_label`` column and ``publication`` column are optional
  but recommended for provenance.

This format is human-readable and mirrors how USGS report tables are typically
typeset.  A single CSV can encode equations for multiple states and regions
simultaneously, enabling a nationwide equation database in one file.

Multi-state example (3 rows from a hypothetical nationwide CSV)::

    region_code,state,aep,intercept,DRNAREA,CSL1085LFP,ELEV,PRECIP,sep_pct
    TN-2,TN,0.01,2.40,0.75,0.38,,,34.0
    GA-1,GA,0.01,1.82,0.80,,0.31,,38.0
    MT-1,MT,0.01,-0.85,0.76,,0.48,0.62,45.0

"""

from __future__ import annotations

import csv
import logging
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

from hydrolib.regression.basin_chars import BasinCharacteristics
from hydrolib.regression.region import HydrologicRegion, RegionRegistry

log = logging.getLogger(__name__)

# Standard AEPs and return periods used by most USGS regional regression studies
STANDARD_AEPS: Tuple[float, ...] = (0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002)
STANDARD_RETURN_PERIODS: Tuple[float, ...] = (2, 5, 10, 25, 50, 100, 200, 500)

# CSV column names reserved for equation metadata (not predictor coefficients)
_RESERVED_COLS: Set[str] = {
    "region_code",
    "region_label",
    "state",
    "aep",
    "return_period_yr",
    "intercept",
    "sep_pct",
    "pseudo_r2",
    "eyr",
    "publication",
    "notes",
}


# ---------------------------------------------------------------------------
# GlsEquation — one equation for one (region, AEP) pair
# ---------------------------------------------------------------------------


@dataclass
class GlsEquation:
    """
    Single GLS regression equation for one hydrologic region and AEP.

    Parameters
    ----------
    region : HydrologicRegion
        Hydrologic region this equation applies to.
    aep : float
        Annual exceedance probability (0 < aep ≤ 1).
    intercept : float
        Constant term (b₀) of the log10-linear regression.
    coefficients : dict
        Mapping StreamStats predictor code → exponent.  For example::

            {"DRNAREA": 0.75, "CSL1085LFP": 0.38}

        An empty dict is valid for intercept-only equations.
    sep_pct : float, optional
        Average standard error of prediction in percent.
    pseudo_r2 : float, optional
        Pseudo-R² (0–1).
    eyr : float, optional
        Equivalent years of record.

    Examples
    --------
    >>> from hydrolib.regression.region import HydrologicRegion
    >>> mt_mtn = HydrologicRegion("MT-1", "Montana Mountain", "MT",
    ...                           required_predictors=("DRNAREA", "ELEV", "PRECIP"))
    >>> eq = GlsEquation(region=mt_mtn, aep=0.01,
    ...                  intercept=-0.85,
    ...                  coefficients={"DRNAREA": 0.76, "ELEV": 0.48, "PRECIP": 0.62},
    ...                  sep_pct=45.0)
    """

    region: HydrologicRegion
    aep: float
    intercept: float
    coefficients: Dict[str, float] = field(default_factory=dict)
    sep_pct: Optional[float] = None
    pseudo_r2: Optional[float] = None
    eyr: Optional[float] = None

    # ------------------------------------------------------------------
    # Derived properties
    # ------------------------------------------------------------------

    @property
    def return_period(self) -> float:
        """Return period in years (1/AEP)."""
        return 1.0 / self.aep

    @property
    def variance_log10(self) -> Optional[float]:
        """
        Variance of the regression estimate in log10(cfs) space.

        Derived from the average standard error of prediction (SEP%)::

            V = [log10(1 + SEP% / 100)]²

        Returns ``None`` if ``sep_pct`` is not set.
        """
        if self.sep_pct is None:
            return None
        return (math.log10(1.0 + self.sep_pct / 100.0)) ** 2

    @property
    def predictor_codes(self) -> List[str]:
        """Sorted list of predictor codes used in this equation."""
        return sorted(self.coefficients)

    # ------------------------------------------------------------------
    # Computation
    # ------------------------------------------------------------------

    def compute(self, basin: BasinCharacteristics) -> float:
        """
        Apply the regression equation to compute peak flow for *basin*.

        Parameters
        ----------
        basin : BasinCharacteristics
            Basin attributes for the site of interest.  The ``predictors``
            dict must contain values for every code in ``coefficients``.

        Returns
        -------
        float
            Estimated peak flow in cfs.

        Raises
        ------
        ValueError
            If the basin's region code does not match this equation's region.
        KeyError
            If a required predictor is missing from *basin*.
        """
        if basin.region.code != self.region.code:
            raise ValueError(
                f"Basin is in region {basin.region.code!r} but equation "
                f"applies to region {self.region.code!r}"
            )
        log_q = self.intercept
        for code, exp in self.coefficients.items():
            val = basin.predictor_value(code)
            if val <= 0:
                raise ValueError(
                    f"Predictor {code!r} must be > 0 for log10 transform; "
                    f"got {val} for site {basin.site_no!r}"
                )
            log_q += exp * math.log10(val)
        return 10**log_q

    # ------------------------------------------------------------------
    # Serialisation
    # ------------------------------------------------------------------

    def to_dict(self) -> Dict[str, object]:
        """Return a flat dict suitable for writing to CSV."""
        d: Dict[str, object] = {
            "region_code": self.region.code,
            "region_label": self.region.label,
            "state": self.region.state,
            "aep": self.aep,
            "return_period_yr": int(round(self.return_period)),
            "intercept": self.intercept,
        }
        d.update(self.coefficients)
        d["sep_pct"] = self.sep_pct if self.sep_pct is not None else ""
        d["pseudo_r2"] = self.pseudo_r2 if self.pseudo_r2 is not None else ""
        d["eyr"] = self.eyr if self.eyr is not None else ""
        return d


# ---------------------------------------------------------------------------
# RegressionTable — the full multi-region, multi-AEP equation set
# ---------------------------------------------------------------------------


class RegressionTable:
    """
    Publication-agnostic container for GLS regression equations.

    Equations are indexed by ``(region_code, aep)`` and can be loaded from
    a Python dict or a wide-format CSV.  A single table can hold equations
    for multiple states, enabling a nationwide database.

    Parameters
    ----------
    publication : str, optional
        Citation or DOI for the source report (for provenance).
    region_registry : RegionRegistry, optional
        Pre-populated registry of :class:`HydrologicRegion` objects.
        When loading from CSV, regions are created automatically if not
        already registered.

    Examples
    --------
    Build from a dict (preferred programmatic API)::

        from hydrolib.regression.states.tennessee import TN_AREA2
        coeff_dict = {
            (TN_AREA2, 0.01): {
                "intercept": 2.40,
                "coefficients": {"DRNAREA": 0.75, "CSL1085LFP": 0.38},
                "sep_pct": 34.0,
            },
        }
        table = RegressionTable.load_from_dict(coeff_dict)
        q100 = table.estimate(basin, aep=0.01)

    Load from CSV::

        table = RegressionTable.load_from_csv("nationwide_equations.csv")

    Merge tables from multiple states::

        from hydrolib.regression.states.tennessee import build_tennessee_table
        from hydrolib.regression.states.georgia import build_georgia_table
        from hydrolib.regression.states.montana import build_montana_table
        national = RegressionTable.merge(
            build_tennessee_table(),
            build_georgia_table(),
            build_montana_table(),
        )
    """

    def __init__(
        self,
        publication: str = "",
        region_registry: Optional[RegionRegistry] = None,
    ) -> None:
        self.publication = publication
        self._registry = region_registry or RegionRegistry()
        # Key: (region_code: str, aep: float) → GlsEquation
        self._eqs: Dict[Tuple[str, float], GlsEquation] = {}

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def load_from_dict(
        cls,
        coeff_dict: Dict[Tuple[HydrologicRegion, float], dict],
        *,
        publication: str = "",
    ) -> "RegressionTable":
        """
        Build from a nested dict of coefficients.

        Parameters
        ----------
        coeff_dict : dict
            Keys are ``(HydrologicRegion, aep)`` tuples.
            Values are dicts with keys:

            * ``"intercept"`` (float) — constant term.
            * ``"coefficients"`` (dict str→float) — predictor exponents.
            * ``"sep_pct"`` (float) — optional standard error of prediction.
            * ``"pseudo_r2"`` (float) — optional.
            * ``"eyr"`` (float) — optional equivalent years of record.

        publication : str, optional

        Returns
        -------
        RegressionTable

        Examples
        --------
        >>> d = {
        ...     (tn_area2, 0.01): {
        ...         "intercept": 2.40,
        ...         "coefficients": {"DRNAREA": 0.75, "CSL1085LFP": 0.38},
        ...         "sep_pct": 34.0,
        ...     },
        ... }
        >>> table = RegressionTable.load_from_dict(d)
        """
        table = cls(publication=publication)
        for (region, aep), kw in coeff_dict.items():
            eq = GlsEquation(
                region=region,
                aep=aep,
                intercept=float(kw["intercept"]),
                coefficients=dict(kw.get("coefficients", {})),
                sep_pct=kw.get("sep_pct"),
                pseudo_r2=kw.get("pseudo_r2"),
                eyr=kw.get("eyr"),
            )
            table._add(eq)
            table._registry.register(region)
        return table

    @classmethod
    def load_from_csv(
        cls,
        path: str | Path,
        *,
        publication: str = "",
        region_registry: Optional[RegionRegistry] = None,
    ) -> "RegressionTable":
        """
        Load the equation table from a wide-format CSV file.

        The CSV must have at minimum the columns ``region_code``, ``aep``,
        and ``intercept``.  Any column whose name is not in the reserved set
        is treated as a predictor coefficient.

        Optional metadata columns: ``region_label``, ``state``,
        ``publication``, ``sep_pct``, ``pseudo_r2``, ``eyr``, ``notes``.

        Parameters
        ----------
        path : str or Path
        publication : str, optional
            If provided, overrides any ``publication`` column in the CSV.
        region_registry : RegionRegistry, optional
            Existing registry; regions are created automatically for
            any ``region_code`` not already registered.

        Returns
        -------
        RegressionTable
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Regression equation CSV not found: {path}")

        registry = region_registry or RegionRegistry()
        table = cls(publication=publication, region_registry=registry)

        with path.open(newline="") as fh:
            reader = csv.DictReader(fh)
            if reader.fieldnames is None:
                raise ValueError(f"CSV file is empty or has no header: {path}")

            predictor_cols = [c for c in reader.fieldnames if c not in _RESERVED_COLS]

            pub_from_csv = ""
            for row in reader:
                region = _region_from_row(row, registry)
                aep = float(row["aep"])
                intercept = float(row["intercept"])

                coefficients: Dict[str, float] = {}
                for col in predictor_cols:
                    val_str = row.get(col, "").strip()
                    if val_str not in ("", "nan", "NaN", "NA", "N/A"):
                        coefficients[col] = float(val_str)

                def _opt_float(key: str) -> Optional[float]:
                    s = row.get(key, "").strip()
                    return float(s) if s not in ("", "nan", "NaN") else None

                eq = GlsEquation(
                    region=region,
                    aep=aep,
                    intercept=intercept,
                    coefficients=coefficients,
                    sep_pct=_opt_float("sep_pct"),
                    pseudo_r2=_opt_float("pseudo_r2"),
                    eyr=_opt_float("eyr"),
                )
                table._add(eq)
                pub_from_csv = row.get("publication", pub_from_csv)

            if not publication and pub_from_csv:
                table.publication = pub_from_csv

        log.info(
            "Loaded %d equations from %s (%d regions)",
            len(table),
            path,
            len(table._registry),
        )
        return table

    @classmethod
    def merge(cls, *tables: "RegressionTable", publication: str = "") -> "RegressionTable":
        """
        Combine multiple :class:`RegressionTable` instances into one.

        Useful for building a nationwide database from per-state tables::

            national = RegressionTable.merge(
                build_tennessee_table(),
                build_georgia_table(),
                build_montana_table(),
            )

        Parameters
        ----------
        *tables : RegressionTable
            Source tables to merge.  Later tables overwrite earlier ones
            for duplicate ``(region_code, aep)`` keys.
        publication : str, optional
            Publication string for the merged table.  If not provided,
            individual publication strings are concatenated.

        Returns
        -------
        RegressionTable
        """
        combined = cls(publication=publication)
        pub_parts: List[str] = []
        for table in tables:
            for eq in table._eqs.values():
                combined._add(eq)
                combined._registry.register(eq.region)
            if table.publication and table.publication not in pub_parts:
                pub_parts.append(table.publication)
        if not publication and pub_parts:
            combined.publication = "; ".join(pub_parts)
        return combined

    # ------------------------------------------------------------------
    # Mutation / internal
    # ------------------------------------------------------------------

    def _add(self, eq: GlsEquation) -> None:
        """Add one equation (overwrites if key already exists)."""
        self._eqs[(eq.region.code, eq.aep)] = eq

    def add_equation(self, eq: GlsEquation) -> None:
        """Add or replace one :class:`GlsEquation`."""
        self._add(eq)
        self._registry.register(eq.region)

    # ------------------------------------------------------------------
    # Access
    # ------------------------------------------------------------------

    def get_equation(self, region: HydrologicRegion, aep: float) -> GlsEquation:
        """
        Retrieve the equation for a given region and AEP.

        Parameters
        ----------
        region : HydrologicRegion
        aep : float

        Returns
        -------
        GlsEquation

        Raises
        ------
        KeyError
        """
        key = (region.code, aep)
        if key not in self._eqs:
            available = sorted(a for r, a in self._eqs if r == region.code)
            raise KeyError(
                f"No equation for region={region.code!r}, AEP={aep}. "
                f"Available AEPs for this region: {available}"
            )
        return self._eqs[key]

    def get_by_region_code(self, region_code: str, aep: float) -> GlsEquation:
        """Retrieve an equation by region code string."""
        key = (region_code, aep)
        if key not in self._eqs:
            raise KeyError(f"No equation for region_code={region_code!r}, AEP={aep}")
        return self._eqs[key]

    def available_aeps(self, region: HydrologicRegion) -> List[float]:
        """Return sorted AEPs available for *region*."""
        return sorted(aep for r, aep in self._eqs if r == region.code)

    def available_regions(self) -> List[HydrologicRegion]:
        """Return all regions that have at least one equation."""
        codes = {r for r, _ in self._eqs}
        regions = []
        for code in sorted(codes):
            r = self._registry.get(code)
            if r:
                regions.append(r)
        return regions

    def available_states(self) -> List[str]:
        """Return sorted list of unique state abbreviations in the table."""
        return sorted({r.state for r in self.available_regions()})

    def regions_for_state(self, state: str) -> List[HydrologicRegion]:
        """Return all regions for a given state abbreviation."""
        return [r for r in self.available_regions() if r.state == state]

    def get_variance(self, region: HydrologicRegion, aep: float) -> Optional[float]:
        """Return the prediction variance (log10 space) for an equation."""
        return self.get_equation(region, aep).variance_log10

    # ------------------------------------------------------------------
    # Computation
    # ------------------------------------------------------------------

    def estimate(self, basin: BasinCharacteristics, aep: float) -> float:
        """
        Estimate peak flow for *basin* at the given *aep*.

        Parameters
        ----------
        basin : BasinCharacteristics
        aep : float

        Returns
        -------
        float
            Peak flow in cfs.
        """
        return self.get_equation(basin.region, aep).compute(basin)

    def estimate_all_aeps(self, basin: BasinCharacteristics) -> Dict[float, float]:
        """
        Estimate peak flows for all AEPs available for the basin's region.

        Parameters
        ----------
        basin : BasinCharacteristics

        Returns
        -------
        dict
            AEP → peak flow (cfs).
        """
        aeps = self.available_aeps(basin.region)
        results: Dict[float, float] = {}
        for aep in aeps:
            try:
                results[aep] = self.estimate(basin, aep)
            except (KeyError, ValueError) as exc:
                log.warning("Skipping AEP=%.4f for %s: %s", aep, basin.site_no, exc)
        return results

    def batch_estimate(
        self,
        basins: Iterable[BasinCharacteristics],
        aep: float,
    ) -> Dict[str, float]:
        """
        Estimate peak flows for multiple basins at one AEP.

        Basins may span multiple regions within this table.  Sites whose
        region or predictor data are insufficient are silently skipped with
        a warning logged.

        Parameters
        ----------
        basins : iterable of BasinCharacteristics
        aep : float

        Returns
        -------
        dict
            ``{site_no: peak_flow_cfs}`` for each successfully estimated site.

        Examples
        --------
        ::

            results = national_table.batch_estimate(site_list, aep=0.01)
            # → {"03606500": 14730.0, "02330450": 8950.0, ...}
        """
        results: Dict[str, float] = {}
        for basin in basins:
            try:
                results[basin.site_no] = self.estimate(basin, aep)
            except (KeyError, ValueError) as exc:
                log.warning("batch_estimate: skipping site %s: %s", basin.site_no, exc)
        return results

    def filter_by_state(self, state: str) -> "RegressionTable":
        """
        Return a new :class:`RegressionTable` containing only equations for *state*.

        Parameters
        ----------
        state : str
            Two-letter state abbreviation (e.g. ``"TN"``, ``"GA"``, ``"MT"``).

        Returns
        -------
        RegressionTable
        """
        filtered = RegressionTable(publication=self.publication)
        for (_, _), eq in self._eqs.items():
            if eq.region.state == state:
                filtered._add(eq)
                filtered._registry.register(eq.region)
        return filtered

    # ------------------------------------------------------------------
    # Output helpers
    # ------------------------------------------------------------------

    def summary_table(self, basin: BasinCharacteristics) -> List[dict]:
        """
        Compute estimates for all AEPs and return as a list of records.

        Parameters
        ----------
        basin : BasinCharacteristics

        Returns
        -------
        list of dict
            Keys: ``aep``, ``return_period``, ``flow_cfs``, ``sep_pct``,
            ``pseudo_r2``, ``eyr``.
        """
        records = []
        for aep, q in self.estimate_all_aeps(basin).items():
            eq = self.get_equation(basin.region, aep)
            records.append(
                {
                    "aep": aep,
                    "return_period": eq.return_period,
                    "flow_cfs": q,
                    "sep_pct": eq.sep_pct,
                    "pseudo_r2": eq.pseudo_r2,
                    "eyr": eq.eyr,
                }
            )
        return sorted(records, key=lambda r: r["aep"], reverse=True)

    def to_csv(self, path: str | Path) -> None:
        """
        Write all equations to a wide-format CSV file.

        The resulting CSV can be reloaded with :meth:`load_from_csv`.

        Parameters
        ----------
        path : str or Path
        """
        path = Path(path)
        if not self._eqs:
            log.warning("RegressionTable is empty; writing empty CSV to %s", path)
            path.write_text("")
            return

        all_pred_codes = sorted({code for eq in self._eqs.values() for code in eq.coefficients})
        fieldnames = (
            [
                "region_code",
                "region_label",
                "state",
                "publication",
                "aep",
                "return_period_yr",
                "intercept",
            ]
            + all_pred_codes
            + ["sep_pct", "pseudo_r2", "eyr"]
        )

        with path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            for (_, _), eq in sorted(self._eqs.items()):
                row = eq.to_dict()
                row["publication"] = self.publication
                for code in all_pred_codes:
                    row.setdefault(code, "")
                writer.writerow(row)

        log.info("Wrote %d equations to %s", len(self._eqs), path)

    # ------------------------------------------------------------------
    # Template helper
    # ------------------------------------------------------------------

    @staticmethod
    def write_template(
        path: str | Path,
        regions: List[HydrologicRegion],
        aeps: Tuple[float, ...] = STANDARD_AEPS,
        extra_predictor_codes: Optional[List[str]] = None,
    ) -> None:
        """
        Write a blank CSV template for users to fill from a published table.

        Parameters
        ----------
        path : str or Path
        regions : list of HydrologicRegion
        aeps : tuple of float
        extra_predictor_codes : list of str, optional
        """
        path = Path(path)
        pred_codes: List[str] = []
        for region in regions:
            for code in region.required_predictors:
                if code not in pred_codes:
                    pred_codes.append(code)
        for code in extra_predictor_codes or []:
            if code not in pred_codes:
                pred_codes.append(code)

        fieldnames = (
            [
                "region_code",
                "region_label",
                "state",
                "publication",
                "aep",
                "return_period_yr",
                "intercept",
            ]
            + pred_codes
            + ["sep_pct", "pseudo_r2", "eyr", "notes"]
        )
        rp_map = dict(zip(aeps, STANDARD_RETURN_PERIODS))

        with path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            for region in regions:
                for aep in aeps:
                    row: Dict[str, object] = {
                        "region_code": region.code,
                        "region_label": region.label,
                        "state": region.state,
                        "publication": region.publication,
                        "aep": aep,
                        "return_period_yr": int(rp_map.get(aep, 1 / aep)),
                        "intercept": "",
                    }
                    for code in pred_codes:
                        row[code] = "" if code in region.required_predictors else "N/A"
                    row.update({"sep_pct": "", "pseudo_r2": "", "eyr": "", "notes": ""})
                    writer.writerow(row)

        log.info("Template with %d rows written to %s", len(regions) * len(aeps), path)

    # ------------------------------------------------------------------
    # Dunder
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return len(self._eqs)

    def __repr__(self) -> str:
        n_regions = len({r for r, _ in self._eqs})
        states = self.available_states()
        return f"RegressionTable(equations={len(self)}, regions={n_regions}, states={states})"


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _region_from_row(row: Dict[str, str], registry: RegionRegistry) -> HydrologicRegion:
    """Obtain or create a HydrologicRegion from a CSV row dict."""
    code = row["region_code"].strip()
    if code in registry:
        return registry[code]

    label = row.get("region_label", code).strip()
    state = row.get("state", "").strip()
    publication = row.get("publication", "").strip()
    notes = row.get("notes", "").strip()
    region = HydrologicRegion(
        code=code,
        label=label,
        state=state,
        # required_predictors not derived from CSV — provide a pre-populated
        # registry if required-predictor enforcement is needed.
        required_predictors=(),
        publication=publication,
        notes=notes,
    )
    registry.register(region)
    return region
