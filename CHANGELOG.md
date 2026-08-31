# Changelog

All notable changes to HydroLib are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `CONTRIBUTING.md`, covering setup, the `make` targets CI runs, the numerical-change
  rules, and the release procedure
- mypy in CI (`make typecheck`) over `hydrolib/`, with fourteen of twenty modules gated
  and six on a documented ratchet in `pyproject.toml`
- Coverage reporting in CI with a 66% floor (`make coverage`), against 68% measured; the
  HTML report is uploaded as a build artifact
- `.github/workflows/release.yml` — tag-driven, refusing to release unless the tag,
  `pyproject.toml` and `.bumpversion.cfg` agree on the version, and asserting the built
  wheel contains `py.typed` and the packaged benchmark data. PyPI publishing is opt-in
  behind a `PUBLISH_TO_PYPI` repository variable and trusted publishing.
- `.github/dependabot.yml` for pip and GitHub Actions updates
- `.pre-commit-config.yaml` pinned to the same tool versions as the dev extra
- `tests/test_batch.py` — first tests for `hydrolib.batch`, which had none

### Changed
- `mypy` and `types-requests` added to the dev extra, so the checker cannot drift from CI

### Fixed
- `analyze_sites()` returned `{"error": "'dict' object has no attribute 'flow'"}` for
  every site. `fetch_nwis_batch` yields plain dicts, `B17CEngine.fit` reads `record.flow`,
  and `run_multi_site`'s `except Exception` swallowed the `AttributeError`. The records are
  now converted at the boundary.
- `Bulletin17C(historical_peaks=..., perception_thresholds=...)` and
  `USGSgage.download_daily_flow(start_date=..., end_date=...)` declared non-optional types
  with a `None` default, so a type checker rejected the documented calls
- `B17CEngine.frequency_table()` and `batch_summary_table()` annotated a return type of
  `"pd.DataFrame"` without importing `pandas` at module level, so the annotations could not
  be resolved by mypy, an IDE, or `typing.get_type_hints()`

### Removed

---

## [0.2.0] - 2026-08-31

### Added
- Instantaneous (unit-value) flow retrieval from USGS NWIS
- Low-flow frequency analysis module (`hydrolib.lowflow`)
  - Annual n-day low-flow frequency with LP3 or lognormal distribution
  - Climatic/water/calendar year definitions
  - Zero-flow-year handling
  - Analytic and bootstrap confidence intervals
- Flow regime metrics module (`hydrolib.regime`)
  - Richards-Baker flashiness index
  - TQmean metric
  - Baseflow separation (UKIH, Lyne-Hollick, HYSEP variants)
  - Monthly and seasonal flow summaries
- Diel (sub-daily) variation analysis
  - Within-day flow range and coefficient of variation
  - Timezone-correct local-day grouping
- Flow series I/O with Parquet backend (`hydrolib.flowio`)
  - Save/load for daily and instantaneous flow data

### Changed
- Improved EMA algorithm convergence handling for edge cases
- Documentation vignettes reorganized (Low-Flow & Flow Regime guide)

### Fixed
- MGBT outlier detection edge case with small sample sizes
- Flow duration curve calculation precision

---

## [0.1.0] - 2026-01-28

### Added
- **USGS Data Retrieval** — Download mean daily, annual peak, and instantaneous flow from NWIS
- **Bulletin 17C Analysis**
  - Expected Moments Algorithm (EMA) — USGS standard method
  - Method of Moments (MOM) fallback
  - Weighted regional skew (MSE weighting per B17C Appendix 6)
  - Multiple Grubbs-Beck test (MGBT) for low outlier detection
  - 90% confidence intervals (5%/95% limits)
- **Hydrograph Plotting**
  - Daily time series plots
  - Summary hydrographs (day of water year with percentile bands)
  - Flow duration curves
- **Frequency Curve Plotting**
  - Log-probability axis
  - LP3 fitted curve with confidence interval band
  - Multi-skew overlay (station / weighted / regional)
- **Streamlit Web Application**
  - Interactive single/multi-gage analysis
  - Regional skew input controls
  - ZIP export (PNG plots, CSV data, LP3 parameters)
  - Multi-gage comparison tables
- **CLI Tools**
  - `hydrolib validate` — EMA validation against reference fixtures
  - `hydrolib benchmark` — Numerical benchmarking (text/JSON output)
- **Technical Reports** — Automated Markdown report generation
- **Validation Framework** — Parity testing against USGS Fortran reference implementation

### Fixed
- Initial release

---

## Notes on Versioning

- **0.x.x** — Pre-release. API may change without warning.
- **1.0.0** — Stable API. Breaking changes require major version bump.

For upgrade guidance, see the [migration guides](docs/) directory.
