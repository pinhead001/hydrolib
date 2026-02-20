"""
hydrolib.regression.sir2024_5130 - GLS regression equations from SIR 2024-5130.

Reference
---------
Ladd, D.E., and Ensminger, P.A., 2025, Estimating the magnitude and
frequency of floods at ungaged locations on streams in Tennessee through
the 2013 water year: U.S. Geological Survey Scientific Investigations
Report 2024-5130, 19 p., https://doi.org/10.3133/sir20245130.

Table 4 Structure
-----------------
Generalized least-squares regression equations of the form::

    log10(Q_T) = b0 + b1 * log10(DA) [+ b2 * log10(S1085)]
                                       [+ b2 * log10(I2)]
                                       [+ b2 * log10(Imperv)]

where:
    Q_T   = peak flow for annual exceedance probability T (cfs)
    DA    = drainage area (sq mi)
    S1085 = 10-85 channel slope (ft/mi) [Area 2 only]
    I2    = 2-yr recurrence precipitation factor       [Area 1 only]
    Imperv= percent impervious cover (%)               [Area 4 only]

Populating Coefficients
-----------------------
Coefficients are NOT hardcoded because the full Table 4 requires
access to the published PDF/HTML (https://doi.org/10.3133/sir20245130).

Populate them by calling::

    eqs = SIR2024_5130.load_from_dict(coeff_dict)

where ``coeff_dict`` has the structure shown in ``EXAMPLE_COEFF_DICT`` below,
or call::

    eqs = SIR2024_5130.load_from_csv("path/to/table4.csv")

A pre-populated instance is available as ``SIR2024_5130.default()`` once
the user has populated coefficients.

CSV Format
----------
The CSV (or TSV) should have columns::

    hydrologic_area, aep, b0, b1, b2, sep_pct, pseudo_r2, eyr

where ``b2`` is the exponent for the secondary predictor (S1085, I2, or
Imperv) and is ``NaN`` / blank for areas with DA-only equations.
"""

from __future__ import annotations

import csv
import logging
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from hydrolib.regression.basin_chars import BasinCharacteristics, HydrologicArea

log = logging.getLogger(__name__)

# Standard AEPs published in SIR 2024-5130 (Table 4 columns)
STANDARD_AEPS: Tuple[float, ...] = (0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002)

# Corresponding return periods (years)
STANDARD_RETURN_PERIODS: Tuple[float, ...] = (2, 5, 10, 25, 50, 100, 200, 500)


@dataclass
class GlsEquation:
    """
    Single GLS regression equation from SIR 2024-5130 Table 4.

    Parameters
    ----------
    hydrologic_area : HydrologicArea
        Hydrologic area this equation applies to.
    aep : float
        Annual exceedance probability (0 < aep <= 1).
    b0 : float
        Intercept of the log10-linear regression equation.
    b1 : float
        Exponent (coefficient) for log10(DA).
    b2 : float, optional
        Exponent for the secondary predictor:
        * Area 1 → log10(I2)
        * Area 2 → log10(S1085)
        * Area 3 → not used (None)
        * Area 4 → log10(Imperv)
    sep_pct : float, optional
        Average standard error of prediction in percent.
    pseudo_r2 : float, optional
        Pseudo-R² (coefficient of determination, 0-1).
    eyr : float, optional
        Equivalent years of record (measure of equation precision).
    """

    hydrologic_area: HydrologicArea
    aep: float
    b0: float
    b1: float
    b2: Optional[float] = None
    sep_pct: Optional[float] = None
    pseudo_r2: Optional[float] = None
    eyr: Optional[float] = None

    # Variance of prediction in log10 units (derived from sep_pct)
    @property
    def variance_log10(self) -> Optional[float]:
        """
        Variance of the regression estimate in log10(cfs) space.

        Derived from the average standard error of prediction (SEP%):
            V = [log10(1 + SEP%/100)]²

        Returns None if sep_pct is not set.
        """
        if self.sep_pct is None:
            return None
        return (math.log10(1.0 + self.sep_pct / 100.0)) ** 2

    @property
    def return_period(self) -> float:
        """Return period in years (1/AEP)."""
        return 1.0 / self.aep

    def compute(self, basin: BasinCharacteristics) -> float:
        """
        Apply the regression equation to compute peak flow.

        Parameters
        ----------
        basin : BasinCharacteristics
            Basin attributes for the site of interest.

        Returns
        -------
        float
            Estimated peak flow in cfs.

        Raises
        ------
        ValueError
            If a required secondary predictor is missing from *basin*.
        """
        if basin.hydrologic_area != self.hydrologic_area:
            raise ValueError(
                f"Basin is in {basin.hydrologic_area.label} but equation is for "
                f"{self.hydrologic_area.label}"
            )

        log_q = self.b0 + self.b1 * math.log10(basin.drainage_area_sqmi)

        if self.b2 is not None:
            secondary = self._get_secondary(basin)
            log_q += self.b2 * math.log10(secondary)

        return 10**log_q

    def _get_secondary(self, basin: BasinCharacteristics) -> float:
        """Extract the required secondary predictor from basin characteristics."""
        area = self.hydrologic_area
        if area == HydrologicArea.AREA1:
            if basin.climate_factor_2yr is None:
                raise ValueError(
                    "climate_factor_2yr is required for Area 1 equation " f"(site {basin.site_no})"
                )
            return basin.climate_factor_2yr
        elif area == HydrologicArea.AREA2:
            if basin.slope_1085_ftmi is None:
                raise ValueError(
                    "slope_1085_ftmi is required for Area 2 equation " f"(site {basin.site_no})"
                )
            return basin.slope_1085_ftmi
        elif area == HydrologicArea.AREA4:
            if basin.pct_impervious is None:
                raise ValueError(
                    "pct_impervious is required for Area 4 equation " f"(site {basin.site_no})"
                )
            return basin.pct_impervious
        else:
            raise ValueError(f"Unexpected secondary predictor request for {area.label}")


# ---------------------------------------------------------------------------
# SIR2024_5130 – the full equation table
# ---------------------------------------------------------------------------


class SIR2024_5130:
    """
    GLS regression equation table from SIR 2024-5130 (Ladd & Ensminger, 2025).

    Holds one :class:`GlsEquation` per (HydrologicArea, AEP) combination.

    Usage
    -----
    Populate the table from a dict or CSV then call :meth:`estimate`::

        coeff_dict = {
            (HydrologicArea.AREA2, 0.01): dict(b0=2.95, b1=0.77, b2=0.38,
                                                sep_pct=38.2, pseudo_r2=0.93),
            ...
        }
        table = SIR2024_5130.load_from_dict(coeff_dict)
        q100 = table.estimate(basin, aep=0.01)

    Alternatively, load from a CSV file::

        table = SIR2024_5130.load_from_csv("sir2024_5130_table4.csv")
    """

    def __init__(self, equations: Optional[Dict[Tuple[HydrologicArea, float], GlsEquation]] = None):
        """
        Parameters
        ----------
        equations : dict, optional
            Mapping (HydrologicArea, aep) -> GlsEquation.
        """
        self._eqs: Dict[Tuple[HydrologicArea, float], GlsEquation] = equations or {}

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def load_from_dict(
        cls,
        coeff_dict: Dict[Tuple[HydrologicArea, float], dict],
    ) -> "SIR2024_5130":
        """
        Build from a nested dictionary of coefficients.

        Parameters
        ----------
        coeff_dict : dict
            Keys are ``(HydrologicArea, aep)`` tuples.
            Values are dicts with keys: ``b0``, ``b1``, and optionally
            ``b2``, ``sep_pct``, ``pseudo_r2``, ``eyr``.

        Returns
        -------
        SIR2024_5130

        Example
        -------
        >>> coeff_dict = {
        ...     (HydrologicArea.AREA2, 0.01): {
        ...         "b0": 2.95, "b1": 0.77, "b2": 0.38,
        ...         "sep_pct": 38.2, "pseudo_r2": 0.93
        ...     },
        ... }
        >>> table = SIR2024_5130.load_from_dict(coeff_dict)
        """
        eqs: Dict[Tuple[HydrologicArea, float], GlsEquation] = {}
        for (area, aep), kw in coeff_dict.items():
            eqs[(area, aep)] = GlsEquation(
                hydrologic_area=area,
                aep=aep,
                b0=kw["b0"],
                b1=kw["b1"],
                b2=kw.get("b2"),
                sep_pct=kw.get("sep_pct"),
                pseudo_r2=kw.get("pseudo_r2"),
                eyr=kw.get("eyr"),
            )
        return cls(eqs)

    @classmethod
    def load_from_csv(cls, path: str | Path) -> "SIR2024_5130":
        """
        Load the equation table from a CSV file.

        The CSV must have at minimum the columns::

            hydrologic_area, aep, b0, b1

        Optional columns: ``b2``, ``sep_pct``, ``pseudo_r2``, ``eyr``.
        ``hydrologic_area`` should be an integer 1-4 or string "area1" etc.

        Parameters
        ----------
        path : str or Path
            Path to the CSV file.

        Returns
        -------
        SIR2024_5130
        """
        eqs: Dict[Tuple[HydrologicArea, float], GlsEquation] = {}
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Table 4 CSV not found: {path}")

        with path.open(newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                raw_area = row["hydrologic_area"].strip()
                if raw_area.isdigit():
                    area = HydrologicArea.from_int(int(raw_area))
                else:
                    area = HydrologicArea(raw_area.lower())

                aep = float(row["aep"])

                def _float(key: str) -> Optional[float]:
                    val = row.get(key, "").strip()
                    if val in ("", "nan", "NaN", "NA", "N/A"):
                        return None
                    return float(val)

                eqs[(area, aep)] = GlsEquation(
                    hydrologic_area=area,
                    aep=aep,
                    b0=float(row["b0"]),
                    b1=float(row["b1"]),
                    b2=_float("b2"),
                    sep_pct=_float("sep_pct"),
                    pseudo_r2=_float("pseudo_r2"),
                    eyr=_float("eyr"),
                )
        log.info("Loaded %d GLS equations from %s", len(eqs), path)
        return cls(eqs)

    # ------------------------------------------------------------------
    # Access
    # ------------------------------------------------------------------

    def get_equation(self, area: HydrologicArea, aep: float) -> GlsEquation:
        """
        Retrieve the equation for a given area and AEP.

        Parameters
        ----------
        area : HydrologicArea
        aep : float
            Annual exceedance probability (e.g., 0.01 for the 1% / 100-yr).

        Returns
        -------
        GlsEquation

        Raises
        ------
        KeyError
            If the (area, aep) combination is not in the table.
        """
        key = (area, aep)
        if key not in self._eqs:
            available = sorted(f"({a.label}, {p})" for a, p in self._eqs)
            raise KeyError(f"No equation for ({area.label}, AEP={aep}). " f"Available: {available}")
        return self._eqs[key]

    def available_aeps(self, area: HydrologicArea) -> List[float]:
        """Return sorted list of AEPs available for the given area."""
        return sorted(aep for a, aep in self._eqs if a == area)

    def available_areas(self) -> List[HydrologicArea]:
        """Return list of hydrologic areas in the table."""
        return list({a for a, _ in self._eqs})

    # ------------------------------------------------------------------
    # Computation
    # ------------------------------------------------------------------

    def estimate(
        self,
        basin: BasinCharacteristics,
        aep: float,
    ) -> float:
        """
        Estimate peak flow for a basin at a given AEP.

        Parameters
        ----------
        basin : BasinCharacteristics
            Basin characteristics for the ungaged site.
        aep : float
            Annual exceedance probability (e.g., 0.01 for 100-yr).

        Returns
        -------
        float
            Peak flow estimate in cfs.
        """
        eq = self.get_equation(basin.hydrologic_area, aep)
        return eq.compute(basin)

    def estimate_all_aeps(
        self,
        basin: BasinCharacteristics,
    ) -> Dict[float, float]:
        """
        Estimate peak flows for all available AEPs.

        Parameters
        ----------
        basin : BasinCharacteristics

        Returns
        -------
        dict
            Mapping AEP -> estimated peak flow (cfs), sorted by AEP.
        """
        aeps = self.available_aeps(basin.hydrologic_area)
        result: Dict[float, float] = {}
        for aep in aeps:
            try:
                result[aep] = self.estimate(basin, aep)
            except (ValueError, KeyError) as exc:
                log.warning("Skipping AEP=%.4f: %s", aep, exc)
        return result

    def get_variance(self, area: HydrologicArea, aep: float) -> Optional[float]:
        """
        Return the prediction variance (log10 space) for an equation.

        Parameters
        ----------
        area : HydrologicArea
        aep : float

        Returns
        -------
        float or None
            Variance in log10(cfs) units, or None if sep_pct not available.
        """
        return self.get_equation(area, aep).variance_log10

    # ------------------------------------------------------------------
    # Reporting helpers
    # ------------------------------------------------------------------

    def summary_table(
        self,
        basin: BasinCharacteristics,
    ) -> "list[dict]":
        """
        Compute estimates for all AEPs and return as a list of records.

        Parameters
        ----------
        basin : BasinCharacteristics

        Returns
        -------
        list of dict
            Each dict has keys: ``aep``, ``return_period``, ``flow_cfs``,
            ``sep_pct``, ``pseudo_r2``, ``eyr``.
        """
        records = []
        for aep, q in self.estimate_all_aeps(basin).items():
            eq = self.get_equation(basin.hydrologic_area, aep)
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

    def __len__(self) -> int:
        return len(self._eqs)

    def __repr__(self) -> str:
        n = len(self._eqs)
        areas = {a.label for a, _ in self._eqs}
        return f"SIR2024_5130(equations={n}, areas={sorted(areas)})"


# ---------------------------------------------------------------------------
# CSV template helper
# ---------------------------------------------------------------------------


def write_table4_template(path: str | Path) -> None:
    """
    Write a blank CSV template for users to populate from SIR 2024-5130 Table 4.

    Parameters
    ----------
    path : str or Path
        Output path for the template CSV file.

    Notes
    -----
    Open the USGS report PDF at https://doi.org/10.3133/sir20245130
    and fill in the coefficients from Table 4 for each row.
    ``b2`` should be left blank for Area 3 (DA-only equation).
    """
    path = Path(path)
    headers = [
        "hydrologic_area",
        "aep",
        "return_period_yr",
        "b0",
        "b1",
        "b2",
        "sep_pct",
        "pseudo_r2",
        "eyr",
        "notes",
    ]
    rows = []
    for area_int in range(1, 5):
        area_label = f"area{area_int}"
        b2_note = {
            1: "b2=log10(I2) coef",
            2: "b2=log10(S1085) coef",
            3: "leave blank",
            4: "b2=log10(Imperv) coef",
        }[area_int]
        for aep, rp in zip(STANDARD_AEPS, STANDARD_RETURN_PERIODS):
            rows.append(
                {
                    "hydrologic_area": area_label,
                    "aep": aep,
                    "return_period_yr": int(rp),
                    "b0": "",
                    "b1": "",
                    "b2": "",
                    "sep_pct": "",
                    "pseudo_r2": "",
                    "eyr": "",
                    "notes": b2_note,
                }
            )

    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)

    log.info("Table 4 template written to %s", path)
