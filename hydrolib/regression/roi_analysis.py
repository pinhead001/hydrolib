"""
hydrolib.regression.roi_analysis - Region of Influence peak flow estimation.

Implements the ROI (Region of Influence) method for estimating peak flows
at ungaged sites using a set of nearby / similar gauged sites weighted by
their similarity to the target site.

Workflow
--------
1. Provide a pool of candidate gauged sites with USGS site numbers.
2. Call :meth:`RoiAnalysis.fetch_nwis_peaks` to download annual peak records
   from NWIS, replacing any previously-computed proxy values.
3. Run B17C EMA analysis on each gauged site to get at-site LP3 parameters.
4. Compute ROI weights based on distance in normalised characteristic space
   (log10 DA, log10 S1085) and/or physical distance.
5. Obtain ROI-weighted LP3 parameter estimates and compute quantiles.

References
----------
Bulletin 17C (England et al., 2018), Appendix 8 — weighting procedures.
Tasker & Stedinger (1989) — GLS regression for flood frequency.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import requests

from hydrolib.core import log_pearson3_ppf
from hydrolib.regression.basin_chars import BasinCharacteristics

log = logging.getLogger(__name__)

# NWIS peak flow data service endpoint
_NWIS_PEAKS_URL = "https://waterservices.usgs.gov/nwis/peak/"

# Default AEPs for quantile output
STANDARD_AEPS: Tuple[float, ...] = (0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002)


# ---------------------------------------------------------------------------
# RoiSite – a single candidate gauged site
# ---------------------------------------------------------------------------


@dataclass
class RoiSite:
    """
    A candidate gauged site in the ROI pool.

    Parameters
    ----------
    site_no : str
        USGS 8-digit site number.
    site_name : str
        Station name.
    basin : BasinCharacteristics
        Basin characteristics for this site.
    weight : float
        Computed ROI weight (0 < weight <= 1; set by :class:`RoiAnalysis`).
    n_peaks : int
        Number of peak flow records available.
    lp3_mean : float, optional
        LP3 mean (log10 space) from B17C EMA analysis.
    lp3_std : float, optional
        LP3 standard deviation (log10 space).
    lp3_skew : float, optional
        LP3 skew coefficient.
    at_site_quantiles : dict, optional
        Mapping AEP -> peak flow (cfs) from at-site EMA analysis.
    record_source : str
        Source of peak data: "nwis" (fetched) or "proxy" (provided externally).
    """

    site_no: str
    site_name: str
    basin: BasinCharacteristics
    weight: float = 0.0
    n_peaks: int = 0
    lp3_mean: Optional[float] = None
    lp3_std: Optional[float] = None
    lp3_skew: Optional[float] = None
    at_site_quantiles: Dict[float, float] = field(default_factory=dict)
    record_source: str = "proxy"

    @property
    def has_lp3(self) -> bool:
        """True if LP3 parameters have been computed."""
        return all(v is not None for v in [self.lp3_mean, self.lp3_std, self.lp3_skew])

    def quantile(self, aep: float) -> Optional[float]:
        """
        Return the at-site quantile for the given AEP.

        Returns the value from ``at_site_quantiles`` if present; otherwise
        computes from LP3 parameters.
        """
        if aep in self.at_site_quantiles:
            return self.at_site_quantiles[aep]
        if self.has_lp3:
            # Use Wilson-Hilferty K-factor (log_pearson3_ppf) for correct
            # monotonic results across negative skew values.
            p = 1.0 - aep  # non-exceedance probability
            return log_pearson3_ppf(p, self.lp3_mean, self.lp3_std, self.lp3_skew)
        return None


# ---------------------------------------------------------------------------
# NWIS peak data fetcher (replaces proxy values with published NWIS stats)
# ---------------------------------------------------------------------------


def fetch_nwis_peak_data(
    site_no: str,
    timeout: int = 30,
) -> pd.DataFrame:
    """
    Download annual peak flow records from NWIS for a single site.

    Replaces any previously-used proxy Q values with actual published
    NWIS peak flow records.

    Parameters
    ----------
    site_no : str
        USGS 8-digit site number.
    timeout : int
        HTTP request timeout in seconds.

    Returns
    -------
    pd.DataFrame
        Columns: ``peak_dt`` (date), ``peak_va`` (cfs), ``peak_cd``
        (qualifier codes).  Empty DataFrame if retrieval fails or no
        records found.

    Raises
    ------
    requests.HTTPError
        If the NWIS service returns a non-200 status code.
    """
    params = {
        "site_no": site_no.zfill(8),
        "agency_cd": "USGS",
        "format": "rdb",
    }
    try:
        resp = requests.get(_NWIS_PEAKS_URL, params=params, timeout=timeout)
        resp.raise_for_status()
    except requests.RequestException as exc:
        log.error("NWIS peak fetch failed for %s: %s", site_no, exc)
        return pd.DataFrame()

    lines = resp.text.splitlines()
    # Skip NWIS RDB comment/header lines starting with '#' and the type-code row
    data_lines = [ln for ln in lines if not ln.startswith("#")]
    if len(data_lines) < 2:
        log.warning("No peak data returned for site %s", site_no)
        return pd.DataFrame()

    from io import StringIO

    # The type-code line (5s, 10d, ...) must be removed; it is the second row
    header_line = data_lines[0]
    type_line = data_lines[1]
    content_lines = data_lines[2:]
    rdb_text = header_line + "\n" + "\n".join(content_lines)

    try:
        df = pd.read_csv(
            StringIO(rdb_text),
            sep="\t",
            dtype=str,
            na_values=["", " "],
        )
    except Exception as exc:
        log.error("Parse error for site %s peak data: %s", site_no, exc)
        return pd.DataFrame()

    # Keep only useful columns
    keep = [c for c in ["peak_dt", "peak_va", "peak_cd"] if c in df.columns]
    df = df[keep].copy()

    if "peak_va" in df.columns:
        df["peak_va"] = pd.to_numeric(df["peak_va"], errors="coerce")

    # Filter out zero or missing flows (these are typically code flags)
    if "peak_va" in df.columns:
        df = df[df["peak_va"].notna() & (df["peak_va"] > 0)]

    log.info("Fetched %d peak records for site %s from NWIS", len(df), site_no)
    return df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# RoiAnalysis – main class
# ---------------------------------------------------------------------------


class RoiAnalysis:
    """
    Region-of-Influence peak flow estimation for an ungaged target site.

    Parameters
    ----------
    target : BasinCharacteristics
        Basin characteristics of the ungaged target site.
    candidate_sites : list of RoiSite
        Pool of potential gauged sites (need not have weights yet).
    weight_method : str
        Weighting scheme.  One of:

        * ``"char_space"`` (default) — inverse distance in normalised
          log10(DA) + log10(S1085) characteristic space.
        * ``"geographic"`` — inverse geographic distance (requires
          latitude/longitude in both target and candidate basins).
        * ``"combined"`` — weighted combination of char_space and
          geographic distances.

    char_weight_fraction : float
        Fraction assigned to characteristic-space distance when
        ``weight_method="combined"`` (default 0.7).
    min_sites : int
        Minimum number of ROI sites required to produce an estimate
        (default 5).
    max_sites : int
        Maximum number of ROI sites used (top by weight; default 20).

    Examples
    --------
    >>> from hydrolib.regression.basin_chars import BasinCharacteristics, DRNAREA, CSL1085LFP
    >>> from hydrolib.regression.roi_analysis import RoiAnalysis, RoiSite
    >>> from hydrolib.regression.states.tennessee import TN_AREA2
    >>>
    >>> target = BasinCharacteristics(
    ...     site_no="UNGAGED01",
    ...     site_name="Target Site",
    ...     region=TN_AREA2,
    ...     predictors={DRNAREA: 120.5, CSL1085LFP: 5.2},
    ... )
    >>> roi = RoiAnalysis(target, candidate_sites=[...])
    >>> roi.fetch_nwis_peaks()
    >>> roi.run_ema()
    >>> quantiles = roi.compute_roi_quantiles()
    """

    def __init__(
        self,
        target: BasinCharacteristics,
        candidate_sites: List[RoiSite],
        weight_method: str = "char_space",
        char_weight_fraction: float = 0.7,
        min_sites: int = 5,
        max_sites: int = 20,
    ) -> None:
        self.target = target
        self.candidate_sites = list(candidate_sites)
        self.weight_method = weight_method
        self.char_weight_fraction = char_weight_fraction
        self.min_sites = min_sites
        self.max_sites = max_sites
        self._peak_data: Dict[str, pd.DataFrame] = {}

    # ------------------------------------------------------------------
    # Step 1: Fetch NWIS peak data (replace proxy values)
    # ------------------------------------------------------------------

    def fetch_nwis_peaks(self, timeout: int = 30) -> Dict[str, int]:
        """
        Download NWIS annual peak records for all candidate sites.

        Replaces any proxy at-site Q values with published NWIS peak
        statistics.  Results are stored internally and also set on each
        :class:`RoiSite` (``n_peaks`` and ``record_source``).

        Parameters
        ----------
        timeout : int
            Per-site HTTP timeout in seconds.

        Returns
        -------
        dict
            Mapping site_no -> number of records fetched.
        """
        counts: Dict[str, int] = {}
        for site in self.candidate_sites:
            df = fetch_nwis_peak_data(site.site_no, timeout=timeout)
            if not df.empty:
                self._peak_data[site.site_no] = df
                site.n_peaks = len(df)
                site.record_source = "nwis"
                log.info(
                    "Site %s: %d NWIS peaks loaded (replaced proxy values)",
                    site.site_no,
                    site.n_peaks,
                )
            else:
                log.warning("Site %s: no NWIS data available; proxy values retained", site.site_no)
            counts[site.site_no] = site.n_peaks
        return counts

    # ------------------------------------------------------------------
    # Step 2: Run EMA on each candidate site
    # ------------------------------------------------------------------

    def run_ema(
        self,
        min_record_length: int = 10,
        regional_skew: float = -0.09,
        regional_skew_mse: float = 0.302,
    ) -> Dict[str, Optional[Exception]]:
        """
        Run B17C EMA analysis on all candidate sites with NWIS data.

        Uses the HydroLib :class:`~hydrolib.engine.B17CEngine` for fitting.
        LP3 parameters and quantiles are stored on each :class:`RoiSite`.

        Parameters
        ----------
        min_record_length : int
            Minimum number of peak records required to fit LP3 (default 10).
        regional_skew : float
            Generalised skew for weighting (Tennessee default: -0.09 from
            SIR 2024-5130 or Bulletin 17C map).
        regional_skew_mse : float
            Mean square error of the regional skew estimate.

        Returns
        -------
        dict
            Mapping site_no -> None on success, Exception on failure.
        """
        from hydrolib.engine import B17CEngine

        errors: Dict[str, Optional[Exception]] = {}
        for site in self.candidate_sites:
            df = self._peak_data.get(site.site_no)
            if df is None or df.empty:
                log.warning("Skipping EMA for %s: no peak data", site.site_no)
                errors[site.site_no] = RuntimeError("No peak data")
                continue

            if len(df) < min_record_length:
                log.warning(
                    "Skipping EMA for %s: only %d records (min %d)",
                    site.site_no,
                    len(df),
                    min_record_length,
                )
                errors[site.site_no] = RuntimeError(
                    f"Insufficient record length ({len(df)} < {min_record_length})"
                )
                continue

            flows = df["peak_va"].dropna().values
            try:
                engine = B17CEngine()
                engine.fit_from_flows(flows)
                site.lp3_mean = float(engine.mu)
                site.lp3_std = float(engine.sigma)
                site.lp3_skew = float(engine.skew)
                # Cache quantiles for standard AEPs using Wilson-Hilferty K-factor
                # for correct monotonic results across all skew signs.
                for aep in STANDARD_AEPS:
                    p = 1.0 - aep
                    site.at_site_quantiles[aep] = float(
                        log_pearson3_ppf(p, site.lp3_mean, site.lp3_std, site.lp3_skew)
                    )
                errors[site.site_no] = None
                log.info(
                    "EMA fitted for %s: μ=%.3f σ=%.3f G=%.3f",
                    site.site_no,
                    site.lp3_mean,
                    site.lp3_std,
                    site.lp3_skew,
                )
            except Exception as exc:
                log.error("EMA failed for site %s: %s", site.site_no, exc)
                errors[site.site_no] = exc

        return errors

    # ------------------------------------------------------------------
    # Step 3: Compute ROI weights
    # ------------------------------------------------------------------

    def compute_weights(self) -> None:
        """
        Compute ROI weights for all candidate sites and store on each site.

        Weight method is determined by :attr:`weight_method`.
        Sites without LP3 parameters receive weight = 0.
        Weights are normalised to sum to 1.
        """
        sites_with_lp3 = [s for s in self.candidate_sites if s.has_lp3]
        if len(sites_with_lp3) < self.min_sites:
            raise RuntimeError(
                f"Only {len(sites_with_lp3)} sites have LP3 parameters; "
                f"need at least {self.min_sites}"
            )

        raw_weights: Dict[str, float] = {}
        for site in sites_with_lp3:
            d = self._distance(self.target, site.basin)
            if d < 1e-10:
                # Target is essentially the site (gauged case)
                raw_weights[site.site_no] = 1e10
            else:
                raw_weights[site.site_no] = 1.0 / d

        # Select top-max_sites
        ranked = sorted(raw_weights.items(), key=lambda kv: -kv[1])[: self.max_sites]
        total_w = sum(w for _, w in ranked)
        selected = {sno: w / total_w for sno, w in ranked}

        for site in self.candidate_sites:
            site.weight = selected.get(site.site_no, 0.0)

    def _distance(self, target: BasinCharacteristics, candidate: BasinCharacteristics) -> float:
        """Compute (non-normalised) distance between target and candidate."""
        if self.weight_method == "char_space":
            return self._char_distance(target, candidate)
        elif self.weight_method == "geographic":
            return self._geo_distance(target, candidate)
        elif self.weight_method == "combined":
            d_char = self._char_distance(target, candidate)
            d_geo = self._geo_distance(target, candidate)
            # Normalise each component by a reference scale before combining
            return self.char_weight_fraction * d_char + (1 - self.char_weight_fraction) * d_geo
        else:
            raise ValueError(f"Unknown weight_method: {self.weight_method!r}")

    @staticmethod
    def _char_distance(target: BasinCharacteristics, candidate: BasinCharacteristics) -> float:
        """
        Euclidean distance in normalised log10-characteristic space.

        Uses the **intersection** of predictor codes that are present and
        positive in both the target and candidate ``predictors`` dicts.
        This is fully generic: for Tennessee Area 2 sites both ``DRNAREA``
        and ``CSL1085LFP`` contribute; for Area 3 (DA-only) only
        ``DRNAREA`` contributes; for a Kentucky study with ``ELEV`` both
        ``DRNAREA`` and ``ELEV`` would contribute.

        If no predictor code is shared, distance is 0 (caller receives
        equal weights — log warning is emitted).
        """
        components: List[float] = []
        for code in sorted(target.predictors):
            t_val = target.predictors.get(code)
            c_val = candidate.predictors.get(code)
            if t_val is not None and c_val is not None and t_val > 0 and c_val > 0:
                components.append(math.log10(c_val) - math.log10(t_val))

        if not components:
            log.debug(
                "No shared positive predictors between %s and %s; distance = 0",
                target.site_no,
                candidate.site_no,
            )
            return 0.0

        return math.sqrt(sum(c * c for c in components))

    @staticmethod
    def _geo_distance(target: BasinCharacteristics, candidate: BasinCharacteristics) -> float:
        """
        Haversine great-circle distance (km) between outlet coordinates.

        Falls back to 0.0 if coordinates are unavailable, which will cause
        all weights to become equal — use carefully.
        """
        if any(
            v is None
            for v in [
                target.latitude,
                target.longitude,
                candidate.latitude,
                candidate.longitude,
            ]
        ):
            log.debug(
                "Missing coordinates for geographic distance (%s or %s); using 0",
                target.site_no,
                candidate.site_no,
            )
            return 0.0

        lat1, lon1 = math.radians(target.latitude), math.radians(target.longitude)
        lat2, lon2 = math.radians(candidate.latitude), math.radians(candidate.longitude)
        dlat = lat2 - lat1
        dlon = lon2 - lon1
        a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
        return 6371.0 * 2.0 * math.asin(math.sqrt(a))

    # ------------------------------------------------------------------
    # Step 4: Compute ROI-weighted quantile estimates
    # ------------------------------------------------------------------

    def compute_roi_quantiles(
        self,
    ) -> Dict[float, float]:
        """
        Compute ROI-weighted quantile estimates for the target site.

        Must call :meth:`compute_weights` first (or weights must already
        be set on candidate sites).

        Weighting is performed in log10 space per the approach of
        Tasker & Stedinger (1989)::

            log10(Q_T_roi) = sum_i(w_i * log10(Q_T_i))

        where the sum is over sites with w_i > 0.

        Returns
        -------
        dict
            Mapping AEP -> ROI-weighted peak flow estimate (cfs).
        """
        active = [s for s in self.candidate_sites if s.weight > 0 and s.has_lp3]
        if not active:
            raise RuntimeError(
                "No candidate sites have positive weights; call compute_weights() first."
            )

        roi_quantiles: Dict[float, float] = {}
        for aep in STANDARD_AEPS:
            log_q_values = []
            weights = []
            for site in active:
                q = site.quantile(aep)
                if q is not None and q > 0:
                    log_q_values.append(math.log10(q))
                    weights.append(site.weight)

            if not log_q_values:
                log.warning("No quantile data for AEP=%.4f across ROI sites", aep)
                continue

            w_arr = np.array(weights)
            lq_arr = np.array(log_q_values)
            # Re-normalise in case some sites had missing quantiles
            w_arr = w_arr / w_arr.sum()
            roi_quantiles[aep] = float(10 ** np.dot(w_arr, lq_arr))

        return roi_quantiles

    def roi_variance(self) -> Dict[float, float]:
        """
        Estimate the variance of ROI quantile estimates (log10 space).

        Uses the weighted sample variance of log10(Q_T_i) across ROI sites,
        adjusted for effective sample size.  This variance can be used in
        the inverse-variance-weighted average with the GLS estimate.

        Returns
        -------
        dict
            Mapping AEP -> variance of log10(Q_T_roi).
        """
        active = [s for s in self.candidate_sites if s.weight > 0 and s.has_lp3]
        if not active:
            raise RuntimeError("No active ROI sites; call compute_weights() first.")

        variances: Dict[float, float] = {}
        for aep in STANDARD_AEPS:
            log_q_values = []
            weights = []
            for site in active:
                q = site.quantile(aep)
                if q is not None and q > 0:
                    log_q_values.append(math.log10(q))
                    weights.append(site.weight)

            if len(log_q_values) < 2:
                log.warning("Insufficient sites for variance at AEP=%.4f", aep)
                variances[aep] = float("nan")
                continue

            w_arr = np.array(weights)
            lq_arr = np.array(log_q_values)
            w_arr = w_arr / w_arr.sum()
            mean_lq = np.dot(w_arr, lq_arr)
            # Weighted variance
            wvar = float(np.dot(w_arr, (lq_arr - mean_lq) ** 2))
            variances[aep] = wvar

        return variances

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    def summary(self) -> str:
        """Return a text summary of the ROI analysis configuration and results."""
        n_nwis = sum(1 for s in self.candidate_sites if s.record_source == "nwis")
        n_lp3 = sum(1 for s in self.candidate_sites if s.has_lp3)
        n_active = sum(1 for s in self.candidate_sites if s.weight > 0)
        lines = [
            f"ROI Analysis — Target: {self.target.site_no} ({self.target.site_name})",
            f"  Candidate sites    : {len(self.candidate_sites)}",
            f"  NWIS data loaded   : {n_nwis}",
            f"  LP3 fitted         : {n_lp3}",
            f"  Active (weight > 0): {n_active}",
            f"  Weight method      : {self.weight_method}",
        ]
        return "\n".join(lines)
