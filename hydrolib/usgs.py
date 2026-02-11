"""
hydrolib.usgs - USGS data retrieval
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import cached_property
from io import StringIO
from typing import ClassVar, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import requests

from .core import PeakRecord


class USGSgage:
    """Class to handle USGS gage data retrieval and storage."""

    BASE_URL_DAILY: ClassVar[str] = "https://waterservices.usgs.gov/nwis/dv/"
    BASE_URL_PEAKS: ClassVar[str] = "https://nwis.waterdata.usgs.gov/nwis/peak"

    def __init__(self, site_no: str):
        self._site_no = str(site_no).zfill(8)
        self._site_name: Optional[str] = None
        self._drainage_area: Optional[float] = None
        self._daily_data: Optional[pd.DataFrame] = None
        self._peak_data: Optional[pd.DataFrame] = None

    @property
    def site_no(self) -> str:
        return self._site_no

    @property
    def site_name(self) -> Optional[str]:
        return self._site_name

    @site_name.setter
    def site_name(self, value: str):
        self._site_name = value

    @property
    def drainage_area(self) -> Optional[float]:
        return self._drainage_area

    @drainage_area.setter
    def drainage_area(self, value: float):
        self._drainage_area = value

    @property
    def daily_data(self) -> Optional[pd.DataFrame]:
        return self._daily_data

    @daily_data.setter
    def daily_data(self, value: pd.DataFrame):
        self._daily_data = value

    @property
    def peak_data(self) -> Optional[pd.DataFrame]:
        return self._peak_data

    @peak_data.setter
    def peak_data(self, value: pd.DataFrame):
        self._peak_data = value

    @cached_property
    def period_of_record(self) -> Optional[Tuple[int, int]]:
        if self._peak_data is not None:
            return (
                int(self._peak_data["water_year"].min()),
                int(self._peak_data["water_year"].max()),
            )
        return None

    def download_daily_flow(self, start_date: str = None, end_date: str = None) -> pd.DataFrame:
        """Download mean daily streamflow data from USGS."""
        params = {
            "format": "rdb",
            "sites": self._site_no,
            "parameterCd": "00060",
            "statCd": "00003",
        }
        if start_date:
            params["startDT"] = start_date
        if end_date:
            params["endDT"] = end_date

        response = requests.get(self.BASE_URL_DAILY, params=params)
        response.raise_for_status()

        lines = response.text.split("\n")
        data_lines = [l for l in lines if not l.startswith("#") and l.strip()]

        if len(data_lines) < 2:
            raise ValueError(f"No daily data found for site {self._site_no}")

        header_idx = 0
        for i, line in enumerate(data_lines):
            if "datetime" in line.lower():
                header_idx = i
                break

        df = pd.read_csv(StringIO("\n".join(data_lines[header_idx:])), sep="\t", skiprows=[1])

        for line in lines:
            if "#" in line and "TS id" in line:
                name_start = line.find(self._site_no) + len(self._site_no)
                self._site_name = line[name_start:].strip()
                break

        date_col = [c for c in df.columns if "datetime" in c.lower()][0]
        flow_col = [c for c in df.columns if "00060" in c and "cd" not in c.lower()]

        if not flow_col:
            raise ValueError("Flow data column not found")

        flow_col = flow_col[0]

        df["date"] = pd.to_datetime(df[date_col])
        df["flow_cfs"] = pd.to_numeric(df[flow_col], errors="coerce")
        df = df[["date", "flow_cfs"]].dropna()
        df = df.set_index("date")

        self._daily_data = df
        return df

    def download_peak_flow(self) -> pd.DataFrame:
        """Download annual peak streamflow data from USGS."""
        params = {
            "site_no": self._site_no,
            "agency_cd": "USGS",
            "format": "rdb",
        }

        response = requests.get(self.BASE_URL_PEAKS, params=params)
        response.raise_for_status()

        lines = response.text.split("\n")
        data_lines = [l for l in lines if not l.startswith("#") and l.strip()]

        if len(data_lines) < 2:
            raise ValueError(f"No peak flow data found for site {self._site_no}")

        for line in lines:
            if "#" in line:
                if "DRAINAGE AREA" in line.upper():
                    try:
                        parts = line.split(":")[-1].strip()
                        self._drainage_area = float(parts.split()[0])
                    except:
                        pass
                if "STATION NAME" in line.upper():
                    self._site_name = line.split(":")[-1].strip()

        df = pd.read_csv(StringIO("\n".join(data_lines)), sep="\t", skiprows=[1])

        df = df[df["agency_cd"] == "USGS"].copy()
        df["peak_date"] = pd.to_datetime(df["peak_dt"], errors="coerce")
        df["peak_flow_cfs"] = pd.to_numeric(df["peak_va"], errors="coerce")

        df["water_year"] = df["peak_date"].apply(
            lambda x: x.year + 1 if x.month >= 10 else x.year if pd.notna(x) else np.nan
        )

        if "peak_cd" in df.columns:
            df["qualification_code"] = df["peak_cd"].fillna("")
        else:
            df["qualification_code"] = ""

        df = df[["water_year", "peak_date", "peak_flow_cfs", "qualification_code"]].dropna(
            subset=["water_year", "peak_flow_cfs"]
        )
        df["water_year"] = df["water_year"].astype(int)

        self._peak_data = df.reset_index(drop=True)

        if "period_of_record" in self.__dict__:
            del self.__dict__["period_of_record"]

        return self._peak_data

    def __repr__(self) -> str:
        return f"USGSgage(site_no='{self._site_no}', name='{self._site_name}')"

    def get_peak_records(self) -> List[PeakRecord]:
        """
        Get peak flow data as a list of PeakRecord objects.

        Returns
        -------
        list of PeakRecord
            Peak flow records suitable for B17CEngine
        """
        if self._peak_data is None:
            self.download_peak_flow()

        records = []
        for _, row in self._peak_data.iterrows():
            records.append(
                PeakRecord(
                    year=int(row["water_year"]),
                    flow=float(row["peak_flow_cfs"]),
                    source="USGS",
                )
            )
        return records


def fetch_nwis_peaks(site_no: str) -> List[PeakRecord]:
    """
    Fetch peak flow records for a single USGS site.

    Parameters
    ----------
    site_no : str
        USGS site number

    Returns
    -------
    list of PeakRecord
        Peak flow records for the site
    """
    gage = USGSgage(site_no)
    gage.download_peak_flow()
    return gage.get_peak_records()


def fetch_nwis_batch(
    sites: List[str], workers: int = 6
) -> Tuple[Dict[str, List[PeakRecord]], Dict[str, str]]:
    """
    Fetch peak flow records for multiple USGS sites in parallel.

    Parameters
    ----------
    sites : list of str
        USGS site numbers
    workers : int
        Number of parallel workers (default: 6)

    Returns
    -------
    tuple
        (successful_results, errors) where:
        - successful_results: dict mapping site_no to list of PeakRecord
        - errors: dict mapping site_no to error message
    """
    results: Dict[str, List[PeakRecord]] = {}
    errors: Dict[str, str] = {}

    with ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_site = {
            executor.submit(fetch_nwis_peaks, site): site for site in sites
        }

        for future in as_completed(future_to_site):
            site = future_to_site[future]
            try:
                results[site] = future.result()
            except Exception as e:
                errors[site] = str(e)

    return results, errors
