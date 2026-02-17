"""
Subprocess wrapper for PeakfqSA Fortran executable.

Handles file I/O, process execution, and result parsing for running
Bulletin 17C flood frequency analysis via the PeakfqSA binary.

Reference: peakfqr R/fortranWrappers.R for call conventions.
"""

from __future__ import annotations

import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from hydrolib.peakfqsa.config import PeakfqSAConfig, find_peakfqsa
from hydrolib.peakfqsa.io_converters import DataFile, SpecificationFile
from hydrolib.peakfqsa.parsers import PeakfqSAResult

logger = logging.getLogger(__name__)


class PeakfqSAExecutionError(RuntimeError):
    """Raised when PeakfqSA exits with a non-zero return code.

    Parameters
    ----------
    returncode : int
        Process exit code.
    stdout : str
        Standard output from the process.
    stderr : str
        Standard error from the process.
    """

    def __init__(self, returncode: int, stdout: str, stderr: str) -> None:
        msg = (
            f"PeakfqSA exited with code {returncode}.\n" f"stdout:\n{stdout}\n" f"stderr:\n{stderr}"
        )
        super().__init__(msg)
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


class PeakfqSATimeoutError(TimeoutError):
    """Raised when PeakfqSA exceeds the configured timeout."""

    pass


class PeakfqSAParseError(ValueError):
    """Raised when PeakfqSA output cannot be parsed."""

    pass


class PeakfqSAWrapper:
    """Subprocess wrapper for the PeakfqSA Fortran executable.

    Handles file I/O, execution, and result parsing for running
    Bulletin 17C flood frequency analysis.

    Parameters
    ----------
    config : PeakfqSAConfig
        Configuration for executable path, timeout, and temp file handling.
    """

    def __init__(self, config: PeakfqSAConfig) -> None:
        self.config = config
        self._executable: Path | None = config.executable_path

    def run(
        self,
        peaks: dict[int, float],
        historical: dict[int, float] | None = None,
        thresholds: list[dict[str, Any]] | None = None,
        begyear: int = 0,
        endyear: int = 0,
        regional_skew: float = -0.302,
        regional_skew_sd: float = 0.55,
        station_name: str = "",
        **kwargs: Any,
    ) -> PeakfqSAResult:
        """Run PeakfqSA analysis on the given peak flow data.

        Parameters
        ----------
        peaks : dict[int, float]
            Systematic annual peaks as {water_year: discharge_cfs}.
        historical : dict[int, float] or None
            Historical flood peaks as {water_year: discharge_cfs}.
        thresholds : list[dict] or None
            Perception thresholds, each with keys: start, end, lower, upper.
        begyear : int
            Beginning year of analysis period.
        endyear : int
            Ending year of analysis period.
        regional_skew : float
            Regional skew coefficient.
        regional_skew_sd : float
            Standard deviation of regional skew.
        station_name : str
            Station name for labeling output.
        **kwargs
            Additional keyword arguments passed to specification file.

        Returns
        -------
        PeakfqSAResult
            Parsed analysis results.

        Raises
        ------
        PeakfqSANotFoundError
            If PeakfqSA executable is not available.
        PeakfqSAExecutionError
            If PeakfqSA exits with non-zero code.
        PeakfqSATimeoutError
            If PeakfqSA exceeds timeout.
        PeakfqSAParseError
            If output cannot be parsed.
        """
        # TODO: Implement full run workflow
        raise NotImplementedError("PeakfqSAWrapper.run() not yet implemented")

    def _write_input_files(
        self,
        tmp_dir: Path,
        peaks: dict[int, float],
        historical: dict[int, float] | None,
        thresholds: list[dict[str, Any]] | None,
        begyear: int,
        endyear: int,
        regional_skew: float,
        regional_skew_sd: float,
        station_name: str,
        **kwargs: Any,
    ) -> tuple[Path, Path]:
        """Write .psf and .dat input files to the temp directory.

        Returns
        -------
        tuple[Path, Path]
            Paths to the specification file and data file.
        """
        # TODO: Implement input file generation
        raise NotImplementedError

    def _execute(self, spec_file: Path, tmp_dir: Path) -> str:
        """Execute PeakfqSA and return raw output text.

        Returns
        -------
        str
            Raw content of the output file.
        """
        # TODO: Implement subprocess execution
        raise NotImplementedError

    def _parse_output(self, output_file: Path) -> PeakfqSAResult:
        """Parse PeakfqSA .out file into a result object.

        Returns
        -------
        PeakfqSAResult
            Parsed results.
        """
        # TODO: Implement output parsing delegation
        raise NotImplementedError

    def is_available(self) -> bool:
        """Check whether PeakfqSA executable is available.

        Returns
        -------
        bool
            True if the executable can be found and validated.
        """
        try:
            if self._executable is None:
                self._executable = find_peakfqsa()
            return self._executable.is_file()
        except Exception:
            return False
