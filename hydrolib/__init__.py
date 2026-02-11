"""
hydrolib - Python library for hydrologic analysis

Includes:
- USGS gage data download (daily and peak flows)
- Summary hydrograph plotting
- Bulletin 17C flood frequency analysis:
  - Method of Moments (MOM)
  - Expected Moments Algorithm (EMA) with full PeakFQ parity
  - B17CEngine for simplified PeakFQ-style analysis
- Multi-site batch processing
- Historical flood information handling
- Technical report generation
"""

from .batch import analyze_sites, batch_summary_table, run_multi_site
from .bulletin17c import (
    Bulletin17C,
    ExpectedMomentsAlgorithm,
    FloodFrequencyAnalysis,
    MethodOfMoments,
)
from .core import (
    AnalysisMethod,
    EMAParameters,
    FlowInterval,
    FrequencyResults,
    PeakRecord,
    SkewMethod,
    compute_ci_lp3,
    grubbs_beck_critical_value,
    kfactor,
    kfactor_array,
    lp3_frequency_factor_peakfq,
    lp3_quantile_peakfq,
)
from .engine import B17CEngine, STANDARD_RETURN_PERIODS
from .hydrograph import Hydrograph
from .plots import (
    apply_b17c_style,
    plot_frequency_curve,
    plot_multi_site_comparison,
    plotting_positions,
)
from .report import HydroReport
from .usgs import USGSgage, fetch_nwis_batch, fetch_nwis_peaks


def analyze_gage(
    site_no: str,
    method: str = "ema",
    regional_skew: float = None,
    regional_skew_mse: float = None,
    historical_peaks: list = None,
    output_dir: str = "./output",
) -> dict:
    """
    Complete flood frequency analysis for a USGS gage.

    Parameters
    ----------
    site_no : str
        USGS site number
    method : str
        'mom' or 'ema' (default: 'ema')
    regional_skew : float, optional
        Regional skew coefficient
    regional_skew_mse : float, optional
        Mean squared error of regional skew
    historical_peaks : list of (year, flow) tuples, optional
        Historical peak observations
    output_dir : str
        Output directory
    """
    import os

    os.makedirs(output_dir, exist_ok=True)

    print(f"Downloading data for USGS {site_no}...")
    gage = USGSgage(site_no)

    try:
        gage.download_daily_flow()
        print(f"  Downloaded {len(gage.daily_data)} days of daily flow data")
    except Exception as e:
        print(f"  Warning: Could not download daily flow data: {e}")

    gage.download_peak_flow()
    print(f"  Downloaded {len(gage.peak_data)} annual peak flow records")
    print(f"  Site name: {gage.site_name}")

    print(f"\nRunning Bulletin 17C analysis (method={method.upper()})...")

    water_years = gage.peak_data["water_year"].values

    analysis = Bulletin17C(
        gage.peak_data["peak_flow_cfs"].values,
        water_years=water_years,
        regional_skew=regional_skew,
        regional_skew_mse=regional_skew_mse,
        historical_peaks=historical_peaks,
    )

    results = analysis.run_analysis(method=method)

    print(f"  Station skew: {results.skew_station:.4f}")
    if results.skew_weighted is not None:
        print(f"  Weighted skew: {results.skew_weighted:.4f}")
    print(f"  Low outlier threshold: {results.low_outlier_threshold:,.0f} cfs")

    if results.method == AnalysisMethod.EMA:
        print(f"  EMA iterations: {results.ema_iterations}")
        print(f"  EMA converged: {results.ema_converged}")

    print("\nGenerating report and figures...")
    report = HydroReport(gage, analysis)
    figures = report.generate_all_figures(output_dir)

    report_path = os.path.join(output_dir, "flood_frequency_report.md")
    report.save_report(report_path)
    print(f"  Report saved to: {report_path}")

    return {
        "gage": gage,
        "analysis": analysis,
        "results": results,
        "figures": figures,
        "report_path": report_path,
    }


__version__ = "0.0.3"
__author__ = "HydroLib"

__all__ = [
    # Core data structures
    "AnalysisMethod",
    "SkewMethod",
    "FlowInterval",
    "EMAParameters",
    "FrequencyResults",
    "PeakRecord",
    # Core functions
    "kfactor",
    "kfactor_array",
    "grubbs_beck_critical_value",
    "lp3_frequency_factor_peakfq",
    "lp3_quantile_peakfq",
    "compute_ci_lp3",
    # USGS data retrieval
    "USGSgage",
    "fetch_nwis_peaks",
    "fetch_nwis_batch",
    # Hydrograph
    "Hydrograph",
    # Bulletin 17C analysis
    "Bulletin17C",
    "MethodOfMoments",
    "ExpectedMomentsAlgorithm",
    "FloodFrequencyAnalysis",
    # B17C Engine (PeakFQ-style)
    "B17CEngine",
    "STANDARD_RETURN_PERIODS",
    # Batch processing
    "run_multi_site",
    "analyze_sites",
    "batch_summary_table",
    # Plotting
    "apply_b17c_style",
    "plot_frequency_curve",
    "plot_multi_site_comparison",
    "plotting_positions",
    # Report
    "HydroReport",
    # Convenience
    "analyze_gage",
]
