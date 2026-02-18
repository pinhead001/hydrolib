"""
Streamlit app for USGS Daily Flow Hydrograph Generator

Run with: streamlit run app/streamlit_app.py
"""

import io
import sys
import zipfile
from datetime import datetime
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

matplotlib.use("Agg")  # Use non-interactive backend

# Add parent directory to path for hydrolib import (needed for Streamlit Cloud)
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.ffa_export import export_comparison_csv, export_ffa_to_zip
from app.ffa_runner import (
    B17C_DEFAULT_SKEW,
    SKEW_OPTIONS,
    build_skew_curves_dict,
    build_station_summary_df,
    compute_skew_tables,
    format_parameters_df,
    format_quantile_df,
    run_ffa,
)

# Import hydrolib
from hydrolib import Hydrograph
from hydrolib.freq_plot import plot_frequency_curve_streamlit, plot_peak_flows_with_thresholds
from hydrolib.usgs import USGSgage

st.set_page_config(
    page_title="USGS Hydrograph-erator",
    page_icon=":ocean:",
    layout="wide",
)

st.title("USGS Hydrograph-erator")
st.markdown(
    "Generate daily flow plots for USGS gages using [hydrolib](https://github.com/pinhead001/hydrolib)"
)

# ---------------------------------------------------------------------------
# Session state — initialize before any sidebar widgets reference it
# ---------------------------------------------------------------------------
if "gage_data" not in st.session_state:
    st.session_state.gage_data = {}
if "gage_info" not in st.session_state:
    st.session_state.gage_info = {}
if "figures" not in st.session_state:
    st.session_state.figures = {}
if "expanded_gage" not in st.session_state:
    st.session_state.expanded_gage = None
if "expanded_plot_idx" not in st.session_state:
    st.session_state.expanded_plot_idx = 0
if "peak_data" not in st.session_state:
    st.session_state.peak_data = {}
if "ffa_results" not in st.session_state:
    st.session_state.ffa_results = {}
if "ffa_results_with_thresholds" not in st.session_state:
    st.session_state.ffa_results_with_thresholds = {}
if "perception_thresholds" not in st.session_state:
    st.session_state.perception_thresholds = {}
if "ffa_year_range" not in st.session_state:
    st.session_state.ffa_year_range = {}  # {site_no: (start_yr, end_yr)}

# ---------------------------------------------------------------------------
# Sidebar — Input Parameters
# ---------------------------------------------------------------------------
st.sidebar.header("Input Parameters")

input_mode = st.sidebar.radio("Input Mode", ["Single Gage", "Multiple Gages"])

if input_mode == "Single Gage":
    gage_input = st.sidebar.text_input("USGS Gage Number", value="09355500")
    gage_list = [gage_input.strip()] if gage_input else []
else:
    gage_input = st.sidebar.text_area(
        "USGS Gage Numbers (one per line or comma-separated)",
        value="09355500\n11066460",
        height=100,
    )
    gage_list = [g.strip() for g in gage_input.replace(",", "\n").split("\n") if g.strip()]

st.sidebar.markdown(f"**{len(gage_list)} gage(s) selected**")

# ---------------------------------------------------------------------------
# Sidebar — Plot Options
# ---------------------------------------------------------------------------
st.sidebar.header("Plot Options")
show_timeseries = st.sidebar.checkbox("Daily Time Series", value=True)
show_summary = st.sidebar.checkbox("Summary Hydrograph", value=True)
show_fdc = st.sidebar.checkbox("Flow Duration Curve", value=True)

# ---------------------------------------------------------------------------
# Sidebar — Plot Date Range (appears above FFA when data is loaded)
# ---------------------------------------------------------------------------
if st.session_state.gage_data:
    st.sidebar.markdown("---")
    st.sidebar.header("Plot Date Range")
    st.sidebar.markdown("*Filter daily flow plots without re-downloading*")

    all_starts = [info["por_start"] for info in st.session_state.gage_info.values()]
    all_ends = [info["por_end"] for info in st.session_state.gage_info.values()]
    overall_start = min(all_starts)
    overall_end = max(all_ends)

    col1, col2 = st.sidebar.columns(2)
    plot_start = col1.date_input(
        "Start Date",
        value=overall_start,
        min_value=overall_start,
        max_value=overall_end,
        key="plot_start",
    )
    plot_end = col2.date_input(
        "End Date",
        value=overall_end,
        min_value=overall_start,
        max_value=overall_end,
        key="plot_end",
    )

    if st.sidebar.button("Update Plots", type="secondary"):
        st.session_state.figures = {}

# ---------------------------------------------------------------------------
# Sidebar — Flood Frequency Analysis
# ---------------------------------------------------------------------------
st.sidebar.markdown("---")
st.sidebar.subheader("Flood Frequency Analysis")
enable_ffa = st.sidebar.checkbox("Enable Flood Frequency Analysis", value=False)
if enable_ffa:
    regional_skew = st.sidebar.number_input(
        "Regional Skew",
        value=-0.302,
        step=0.001,
        format="%.3f",
        help="Generalized skew from Bulletin 17C Appendix or regional study. Nationwide default: -0.302",
    )
    regional_skew_se = st.sidebar.number_input(
        "Regional Skew SE",
        value=0.55,
        step=0.01,
        format="%.2f",
        help="Standard error of the regional skew estimate. Nationwide default: 0.55",
    )
    _auto_skew_source = (
        "B17C 2019 (Nationwide)"
        if abs(regional_skew - B17C_DEFAULT_SKEW) < 5e-4
        else "User-defined"
    )
    map_skew_source = st.sidebar.text_input(
        "Map Skew Source",
        value=_auto_skew_source,
        help=(
            "Source for the regional skew value. Auto-detected as "
            "'B17C 2019 (Nationwide)' when using the default (−0.302). "
            "Edit freely to reflect a state or regional study."
        ),
    )
    show_freq_curve = st.sidebar.checkbox("Frequency Curve", value=True)
    st.sidebar.markdown("**Skew Options**")
    skew_station_on = st.sidebar.checkbox("Station Skew", value=False)
    skew_weighted_on = st.sidebar.checkbox("Weighted Skew", value=True)
    skew_regional_on = st.sidebar.checkbox("Regional Skew", value=False)
    show_threshold_curve = st.sidebar.checkbox(
        "Perception Thresholds",
        value=False,
        help="Overlay the EMA curve re-run with user-defined perception thresholds.",
    )
    st.sidebar.markdown("**PILF (Low Outlier) Override**")
    lo_override_input = st.sidebar.number_input(
        "PILF Threshold Override (cfs)",
        value=0.0,
        min_value=0.0,
        step=100.0,
        format="%.0f",
        help=(
            "Leave at 0 to use the automatic Multiple Grubbs-Beck Test (MGBT) result. "
            "Enter a positive value to override the MGBT and censor all peaks below "
            "this threshold."
        ),
    )
    low_outlier_override = float(lo_override_input) if lo_override_input > 0 else None
else:
    regional_skew = -0.302
    regional_skew_se = 0.55
    map_skew_source = "B17C 2019 (Nationwide)"
    show_freq_curve = False
    skew_station_on = False
    skew_weighted_on = True
    skew_regional_on = False
    show_threshold_curve = False
    low_outlier_override = None

# Download button
download_data = st.sidebar.button("Download Data", type="primary")


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def get_plot_list(site_no):
    """Get list of available plots for a gage (excluding stats dataframe)."""
    if site_no not in st.session_state.figures:
        return []
    return [
        name
        for name in st.session_state.figures[site_no].keys()
        if not name.endswith("_stats") and not name.endswith("_data")
    ]


def get_plot_display_name(plot_key):
    """Convert plot key to display name."""
    names = {
        "daily_timeseries": "Daily Time Series",
        "summary_hydrograph": "Summary Hydrograph",
        "flow_duration_curve": "Flow Duration Curve",
        "frequency_curve": "Frequency Curve",
    }
    return names.get(plot_key, plot_key)


def format_date(dt):
    """Format datetime as M/D/YYYY string."""
    return f"{dt.month}/{dt.day}/{dt.year}"


def generate_plots(site_no, plot_data, gage_info, por_start_str=None, por_end_str=None):
    """Generate plots for a gage using the provided data."""
    site_figs = {}

    if show_timeseries:
        fig1 = Hydrograph.plot_daily_timeseries(
            plot_data,
            site_name=gage_info.get("name"),
            site_no=site_no,
            figsize=(10, 4),
            por_start=por_start_str,
            por_end=por_end_str,
        )
        site_figs["daily_timeseries"] = fig1

    if show_summary:
        fig2 = Hydrograph.plot_summary_hydrograph(
            plot_data,
            site_name=gage_info.get("name"),
            site_no=site_no,
            figsize=(10, 4),
            percentiles=[10, 25, 50, 75, 90],
            por_start=por_start_str,
            por_end=por_end_str,
        )
        site_figs["summary_hydrograph"] = fig2
        site_figs["summary_stats"] = Hydrograph.get_summary_stats(plot_data)

    if show_fdc:
        fig3, stats_df = Hydrograph.plot_flow_duration_curve(
            plot_data,
            site_name=gage_info.get("name"),
            site_no=site_no,
            figsize=(10, 4),
            por_start=por_start_str,
            por_end=por_end_str,
        )
        site_figs["flow_duration_curve"] = fig3
        site_figs["fdc_stats"] = stats_df

    return site_figs


# ---------------------------------------------------------------------------
# Download data when button clicked
# ---------------------------------------------------------------------------
if download_data and gage_list:
    st.session_state.gage_data = {}
    st.session_state.gage_info = {}
    st.session_state.figures = {}
    st.session_state.peak_data = {}
    st.session_state.ffa_results = {}
    st.session_state.ffa_results_with_thresholds = {}
    st.session_state.perception_thresholds = {}
    st.session_state.ffa_year_range = {}

    progress_bar = st.progress(0)
    status_text = st.empty()

    for idx, site_no in enumerate(gage_list):
        status_text.text(f"Fetching site info for {site_no}...")

        try:
            gage = USGSgage(site_no)
            gage.fetch_site_info()

            status_text.text(f"Downloading data for {site_no}...")

            start_dt = gage.daily_por_start if gage.daily_por_start else "1900-01-01"
            end_dt = gage.daily_por_end if gage.daily_por_end else None
            daily_data = gage.download_daily_flow(start_date=start_dt, end_date=end_dt)

            site_name_display = gage.site_name or "Unknown"
            por_start_date = daily_data.index.min()
            por_end_date = daily_data.index.max()

            st.session_state.gage_info[site_no] = {
                "name": site_name_display,
                "drainage_area": gage.drainage_area,
                "por_start": por_start_date,
                "por_end": por_end_date,
                "records": len(daily_data),
                "latitude": gage.latitude,
                "longitude": gage.longitude,
            }

            st.session_state.gage_data[site_no] = daily_data

            if enable_ffa:
                with st.spinner(f"Downloading peak flows for {site_no}..."):
                    try:
                        peak_df = gage.download_peak_flow()
                        st.session_state.peak_data[site_no] = peak_df
                    except Exception as e:
                        st.warning(f"Could not download peak flows for {site_no}: {e}")
                        st.session_state.peak_data[site_no] = None

                if st.session_state.peak_data.get(site_no) is not None:
                    peak_df = st.session_state.peak_data[site_no]
                    with st.spinner(f"Running flood frequency analysis for {site_no}..."):
                        ffa_result = run_ffa(
                            peak_flows=peak_df["peak_flow_cfs"].values,
                            water_years=peak_df["water_year"].values.astype(int),
                            regional_skew=regional_skew,
                            regional_skew_se=regional_skew_se,
                            low_outlier_threshold_override=low_outlier_override,
                        )
                        st.session_state.ffa_results[site_no] = ffa_result
                    if ffa_result.get("error"):
                        st.warning(f"FFA error for {site_no}: {ffa_result['error']}")

        except Exception as e:
            st.error(f"Error processing {site_no}: {str(e)}")

        progress_bar.progress((idx + 1) / len(gage_list))

    status_text.text("Download complete!")
    plt.close("all")
    st.rerun()


# ---------------------------------------------------------------------------
# Main content — generate and display plots for each gage
# ---------------------------------------------------------------------------
if st.session_state.gage_data:
    for site_no in st.session_state.gage_data.keys():
        gage_info = st.session_state.gage_info[site_no]
        full_data = st.session_state.gage_data[site_no]

        # Filter daily data for plot date range
        plot_data = full_data[
            (full_data.index >= pd.Timestamp(plot_start))
            & (full_data.index <= pd.Timestamp(plot_end))
        ]

        if len(plot_data) == 0:
            st.warning(f"No data for {site_no} in selected date range")
            continue

        # ── Gage header ───────────────────────────────────────────────────
        st.subheader(f"USGS {site_no} - {gage_info['name']}")

        por_str = f"{format_date(gage_info['por_start'])} - {format_date(gage_info['por_end'])}"
        plot_str = f"{format_date(plot_data.index.min())} - {format_date(plot_data.index.max())}"
        da_str = f"{gage_info['drainage_area']:,.1f} sq mi" if gage_info["drainage_area"] else "N/A"

        info_cols = st.columns(4)
        info_cols[0].markdown(
            f"**Drainage Area**<br><small>{da_str}</small>", unsafe_allow_html=True
        )
        info_cols[1].markdown(f"**POR**<br><small>{por_str}</small>", unsafe_allow_html=True)
        info_cols[2].markdown(
            f"**Plot Range**<br><small>{plot_str}</small>", unsafe_allow_html=True
        )
        info_cols[3].markdown(
            f"**Records**<br><small>{len(plot_data):,} days</small>",
            unsafe_allow_html=True,
        )

        # ── Daily flow plots ──────────────────────────────────────────────
        if site_no not in st.session_state.figures:
            por_start_str = None
            por_end_str = None
            if (
                plot_data.index.min() > gage_info["por_start"]
                or plot_data.index.max() < gage_info["por_end"]
            ):
                por_start_str = format_date(gage_info["por_start"])
                por_end_str = format_date(gage_info["por_end"])

            st.session_state.figures[site_no] = generate_plots(
                site_no, plot_data, gage_info, por_start_str, por_end_str
            )

        site_figs = st.session_state.figures[site_no]

        plot_keys = [k for k in site_figs.keys() if not k.endswith("_stats") and not k.endswith("_data")]
        if plot_keys:
            cols = st.columns(len(plot_keys))
            for i, plot_key in enumerate(plot_keys):
                with cols[i]:
                    st.pyplot(site_figs[plot_key])

        # ── Peak flow analysis separator ──────────────────────────────────
        _peak_raw = st.session_state.peak_data.get(site_no)
        if enable_ffa and _peak_raw is not None:
            st.divider()
            st.markdown("#### Peak Flow Frequency Analysis")

            # FFA date range controls
            min_wy = int(_peak_raw["water_year"].min())
            max_wy = int(_peak_raw["water_year"].max())
            _stored_wy = st.session_state.ffa_year_range.get(site_no, (min_wy, max_wy))

            col_yr1, col_yr2, col_recalc = st.columns([2, 2, 3])
            ffa_start_yr = col_yr1.number_input(
                "FFA Start Year",
                min_value=min_wy,
                max_value=max_wy,
                value=int(_stored_wy[0]),
                step=1,
                key=f"ffa_start_{site_no}",
            )
            ffa_end_yr = col_yr2.number_input(
                "FFA End Year",
                min_value=min_wy,
                max_value=max_wy,
                value=int(_stored_wy[1]),
                step=1,
                key=f"ffa_end_{site_no}",
            )
            with col_recalc:
                st.write("\u00a0")  # vertical spacer to align button with inputs
                recalc_clicked = st.button(
                    "Recalculate FFA",
                    key=f"recalc_ffa_{site_no}",
                    help="Re-run flood frequency analysis using only peaks within the selected year range.",
                )

            if recalc_clicked:
                _r_mask = (
                    (_peak_raw["water_year"] >= ffa_start_yr)
                    & (_peak_raw["water_year"] <= ffa_end_yr)
                )
                _filtered = _peak_raw[_r_mask]
                if len(_filtered) < 5:
                    st.warning(
                        f"Only {len(_filtered)} peak(s) in {ffa_start_yr}–{ffa_end_yr}. "
                        "Results will be unreliable."
                    )
                else:
                    with st.spinner(
                        f"Re-running FFA for {site_no} ({ffa_start_yr}–{ffa_end_yr})…"
                    ):
                        _recalc = run_ffa(
                            peak_flows=_filtered["peak_flow_cfs"].values,
                            water_years=_filtered["water_year"].values.astype(int),
                            regional_skew=regional_skew,
                            regional_skew_se=regional_skew_se,
                            low_outlier_threshold_override=low_outlier_override,
                        )
                    st.session_state.ffa_results[site_no] = _recalc
                    st.session_state.ffa_year_range[site_no] = (int(ffa_start_yr), int(ffa_end_yr))
                    st.session_state.ffa_results_with_thresholds.pop(site_no, None)
                    if _recalc.get("error"):
                        st.error(f"FFA error: {_recalc['error']}")
                    else:
                        st.rerun()

            # Filtered peak_df used for all FFA computations below
            _wy = st.session_state.ffa_year_range.get(site_no)
            if _wy:
                _ffa_mask = (
                    (_peak_raw["water_year"] >= _wy[0])
                    & (_peak_raw["water_year"] <= _wy[1])
                )
                peak_df_for_ffa = _peak_raw[_ffa_mask]
            else:
                peak_df_for_ffa = _peak_raw
        else:
            peak_df_for_ffa = None

        # ── Peak flow record & perception thresholds expander ─────────────
        if enable_ffa and _peak_raw is not None:
            with st.expander("Peak Flow Record & Perception Thresholds", expanded=False):
                st.markdown(
                    "Define time periods when only peaks within a certain flow range were "
                    "observed and recorded. Enter **Lower (cfs)** = minimum detectable flow "
                    "(use 0 for no lower limit) and **Upper (cfs)** = maximum detectable flow "
                    "(enter **inf** for no upper limit / historical threshold). "
                    "Each period adds censored intervals to EMA — click **Apply** to re-run."
                )

                # Show note when a custom FFA year range is active
                _wy_active = st.session_state.ffa_year_range.get(site_no)
                if _wy_active:
                    st.info(
                        f"FFA is using peaks from **{_wy_active[0]}–{_wy_active[1]}** "
                        f"({len(peak_df_for_ffa)} peaks). "
                        "Hollow bars show peaks excluded from the analysis. "
                        "The threshold analysis will also use the filtered year range."
                    )

                # ── Auto-populate helper ─────────────────────────────────
                def _build_auto_thr_rows(peak_df: pd.DataFrame, ffa_res: dict) -> list[dict]:
                    """Build auto-populate rows: systematic, gaps, and PILFs."""
                    if peak_df.empty:
                        return []
                    rows = []
                    wy = sorted(peak_df["water_year"].astype(int).tolist())
                    min_wy, max_wy = wy[0], wy[-1]
                    flows = peak_df.set_index("water_year")["peak_flow_cfs"]

                    # 1. Systematic record (full range, 0/inf)
                    rows.append({
                        "Start Year": min_wy,
                        "End Year": max_wy,
                        "Lower (cfs)": "0",
                        "Upper (cfs)": "inf",
                        "Note": "Systematic record",
                    })

                    # 2. Gaps in the systematic record
                    wy_set = set(wy)
                    in_gap = False
                    gap_start = None
                    for y in range(min_wy, max_wy + 1):
                        if y not in wy_set:
                            if not in_gap:
                                in_gap = True
                                gap_start = y
                                # Max flow just before gap
                                pre_gap = [f for yr, f in flows.items() if yr < y]
                                pre_max = max(pre_gap) if pre_gap else 0
                        else:
                            if in_gap:
                                rows.append({
                                    "Start Year": gap_start,
                                    "End Year": y - 1,
                                    "Lower (cfs)": "0",
                                    "Upper (cfs)": f"{int(round(pre_max))}",
                                    "Note": "Gap in record",
                                })
                                in_gap = False

                    # 3. PILF years (below MGBT/override threshold)
                    if ffa_res and not ffa_res.get("error") and ffa_res.get("b17c"):
                        r_res = ffa_res["b17c"].results
                        lo_thr = r_res.low_outlier_threshold
                        if lo_thr > 0:
                            thr_str = f"{int(round(lo_thr))}"
                            for yr in peak_df.loc[
                                peak_df["peak_flow_cfs"] < lo_thr, "water_year"
                            ].astype(int):
                                rows.append({
                                    "Start Year": int(yr),
                                    "End Year": int(yr),
                                    "Lower (cfs)": "0",
                                    "Upper (cfs)": thr_str,
                                    "Note": "PILF year",
                                })
                    return rows

                # Auto-populate button
                _cur_ffa = st.session_state.ffa_results.get(site_no)
                auto_populate = st.button(
                    "Auto-populate from record",
                    key=f"auto_thr_{site_no}",
                    help=(
                        "Pre-fill the table with the systematic record, gaps, "
                        "and PILF-identified years. Review and edit before applying."
                    ),
                )
                if auto_populate:
                    st.session_state[f"thr_prepop_{site_no}"] = _build_auto_thr_rows(
                        peak_df_for_ffa, _cur_ffa
                    )
                    st.info(
                        "Table auto-populated. **Please review the values** — especially "
                        "the Upper threshold for gaps — before clicking Apply."
                    )

                # Build initial DataFrame for editor
                _prepop = st.session_state.get(f"thr_prepop_{site_no}", [])
                if _prepop:
                    thr_df_init = pd.DataFrame(_prepop)
                else:
                    thr_df_init = pd.DataFrame(
                        {
                            "Start Year": pd.Series([], dtype="float64"),
                            "End Year": pd.Series([], dtype="float64"),
                            "Lower (cfs)": pd.Series([], dtype="str"),
                            "Upper (cfs)": pd.Series([], dtype="str"),
                            "Note": pd.Series([], dtype="str"),
                        }
                    )

                edited_thr = st.data_editor(
                    thr_df_init,
                    num_rows="dynamic",
                    key=f"thr_editor_{site_no}",
                    column_config={
                        "Start Year": st.column_config.NumberColumn(
                            min_value=1700, max_value=2100, step=1, format="%d",
                        ),
                        "End Year": st.column_config.NumberColumn(
                            min_value=1700, max_value=2100, step=1, format="%d",
                        ),
                        "Lower (cfs)": st.column_config.TextColumn(
                            help="Minimum detectable flow in cfs. Use 0 for no lower limit.",
                        ),
                        "Upper (cfs)": st.column_config.TextColumn(
                            help="Maximum detectable flow in cfs. Enter 'inf' for no upper limit.",
                        ),
                        "Note": st.column_config.TextColumn(
                            help="Optional description of this threshold period.",
                        ),
                    },
                    use_container_width=True,
                )

                # Parse valid threshold rows (supports 'inf' for upper)
                thresholds = []
                for _, row in edited_thr.dropna(subset=["Start Year", "End Year"]).iterrows():
                    try:
                        lo_raw = str(row.get("Lower (cfs)", "0") or "0").strip()
                        hi_raw = str(row.get("Upper (cfs)", "inf") or "inf").strip()
                        lo = 0.0 if lo_raw in ("", "0") else float(lo_raw)
                        hi = float("inf") if hi_raw.lower() in ("inf", "infinity", "") else float(hi_raw)
                        if lo > 0 or not np.isinf(hi):
                            thresholds.append(
                                {
                                    "start_year": int(row["Start Year"]),
                                    "end_year": int(row["End Year"]),
                                    "lower_cfs": lo,
                                    "upper_cfs": hi,
                                    # Legacy key for ffa_runner (uses upper as the threshold)
                                    "threshold_cfs": hi if not np.isinf(hi) else 0.0,
                                }
                            )
                    except (ValueError, TypeError, KeyError):
                        pass

                st.session_state.perception_thresholds[site_no] = thresholds

                # MGBT threshold for bar chart annotation
                _mgbt_thr = None
                if _cur_ffa and not _cur_ffa.get("error") and _cur_ffa.get("b17c"):
                    _lo = _cur_ffa["b17c"].results.low_outlier_threshold
                    if _lo and _lo > 0:
                        _mgbt_thr = _lo

                # Bar chart — full record; hollow bars = years excluded from FFA
                peak_fig = plot_peak_flows_with_thresholds(
                    _peak_raw,
                    site_name=gage_info.get("name", ""),
                    site_no=site_no,
                    thresholds=thresholds or None,
                    ffa_year_range=st.session_state.ffa_year_range.get(site_no),
                    mgbt_threshold=_mgbt_thr,
                )
                st.pyplot(peak_fig)
                plt.close(peak_fig)

                # Apply button + status — FFA re-run uses filtered year range
                col_btn, col_status = st.columns([1, 2])
                with col_btn:
                    apply_clicked = st.button(
                        "Apply Thresholds to Analysis",
                        key=f"apply_thr_{site_no}",
                        disabled=(len(thresholds) == 0),
                        help="Re-run EMA incorporating the threshold periods above.",
                    )
                if apply_clicked:
                    with st.spinner("Re-running EMA with perception thresholds…"):
                        thr_run = run_ffa(
                            peak_flows=peak_df_for_ffa["peak_flow_cfs"].values,
                            water_years=peak_df_for_ffa["water_year"].values.astype(int),
                            regional_skew=regional_skew,
                            regional_skew_se=regional_skew_se,
                            perception_thresholds=thresholds,
                            low_outlier_threshold_override=low_outlier_override,
                        )
                        st.session_state.ffa_results_with_thresholds[site_no] = thr_run
                    st.rerun()
                with col_status:
                    thr_stored = st.session_state.ffa_results_with_thresholds.get(site_no)
                    if thr_stored:
                        if not thr_stored.get("error") and thr_stored.get("b17c"):
                            r_t = thr_stored["b17c"].results
                            st.success(
                                f"Applied — {len(st.session_state.perception_thresholds.get(site_no, []))} "
                                f"period(s), {r_t.n_peaks} total peaks "
                                f"({r_t.n_systematic} sys + {r_t.n_censored} censored)"
                            )
                        else:
                            st.error(f"Analysis error: {thr_stored.get('error')}")

        # ── Frequency curve plot ──────────────────────────────────────────
        if show_freq_curve and site_no in st.session_state.ffa_results:
            _base_result = st.session_state.ffa_results[site_no]
            _thr_result = st.session_state.ffa_results_with_thresholds.get(site_no)
            _thr_valid = bool(
                _thr_result and not _thr_result.get("error") and _thr_result.get("b17c")
            )

            # Use the threshold-adjusted result as the primary curve when available
            ffa_result = _thr_result if _thr_valid else _base_result

            if not ffa_result.get("error") and ffa_result.get("b17c"):
                selected_skew_labels = [
                    lbl
                    for lbl, on in [
                        ("Station Skew", skew_station_on),
                        ("Weighted Skew", skew_weighted_on),
                        ("Regional Skew", skew_regional_on),
                    ]
                    if on
                ]
                skew_curves = build_skew_curves_dict(ffa_result, selected_skew_labels) or None

                # When threshold result is primary, show baseline as comparison overlay
                extra_curves = None
                if show_threshold_curve and _thr_valid:
                    r_base = _base_result["b17c"].results
                    extra_curves = {
                        "Baseline (No Thresholds)": (
                            r_base.mean_log,
                            r_base.std_log,
                            r_base.skew_used,
                            r_base.n_systematic or r_base.n_peaks,
                        )
                    }

                freq_fig = plot_frequency_curve_streamlit(
                    ffa_result["b17c"],
                    site_name=gage_info.get("name", ""),
                    site_no=site_no,
                    skew_curves=skew_curves,
                    extra_curves=extra_curves,
                )
                st.pyplot(freq_fig)

        # ── FFA results expander ──────────────────────────────────────────
        if enable_ffa and site_no in st.session_state.ffa_results:
            ffa_result = st.session_state.ffa_results[site_no]
            if not ffa_result.get("error"):
                selected_skew_labels = [
                    lbl
                    for lbl, on in [
                        ("Station Skew", skew_station_on),
                        ("Weighted Skew", skew_weighted_on),
                        ("Regional Skew", skew_regional_on),
                    ]
                    if on
                ]
                skew_tables = compute_skew_tables(ffa_result, selected_skew_labels)

                with st.expander("Flood Frequency Results", expanded=False):
                    # Station summary table (PeakFQ-style)
                    if peak_df_for_ffa is not None and ffa_result.get("b17c"):
                        st.markdown("**Station Summary**")
                        primary_label = (
                            selected_skew_labels[0] if selected_skew_labels else "Weighted Skew"
                        )
                        summary_df = build_station_summary_df(
                            site_no=site_no,
                            peak_df=peak_df_for_ffa,
                            ffa_result=ffa_result,
                            regional_skew=regional_skew,
                            regional_skew_se=regional_skew_se,
                            primary_skew_label=primary_label,
                            latitude=gage_info.get("latitude"),
                            longitude=gage_info.get("longitude"),
                            map_skew_source=map_skew_source,
                        )
                        st.dataframe(summary_df, use_container_width=True, hide_index=True)

                    # Convergence badge
                    if ffa_result.get("converged"):
                        st.success("EMA Converged")
                    else:
                        st.warning("MOM Fallback (EMA did not converge)")

                    st.markdown("**LP3 Parameters**")
                    st.dataframe(
                        format_parameters_df(ffa_result["parameters"]),
                        use_container_width=True,
                    )

                    if skew_tables:
                        for label, tbl in skew_tables.items():
                            st.markdown(
                                f"**Frequency Table — {label}** (Return intervals 1.5–500 years)"
                            )
                            st.dataframe(
                                format_quantile_df(tbl),
                                use_container_width=True,
                            )
                    else:
                        st.markdown("**Frequency Table** (Return intervals 1.5–500 years)")
                        st.dataframe(
                            format_quantile_df(ffa_result["quantile_df"]),
                            use_container_width=True,
                        )

                    # Threshold-aware frequency tables
                    if show_threshold_curve:
                        thr_res = st.session_state.ffa_results_with_thresholds.get(site_no)
                        if thr_res and not thr_res.get("error") and thr_res.get("b17c"):
                            r_thr = thr_res["b17c"].results
                            st.markdown("---")
                            st.markdown(
                                f"**Perception Threshold Analysis** — record extended to "
                                f"{r_thr.n_peaks} peaks "
                                f"({r_thr.n_systematic} systematic + {r_thr.n_censored} censored)"
                            )
                            thr_skew_tables = compute_skew_tables(thr_res, selected_skew_labels)
                            if thr_skew_tables:
                                for label, tbl in thr_skew_tables.items():
                                    st.markdown(
                                        f"**Threshold Frequency Table — {label}** "
                                        "(Return intervals 1.5–500 years)"
                                    )
                                    st.dataframe(
                                        format_quantile_df(tbl), use_container_width=True
                                    )
                            else:
                                st.markdown(
                                    "**Threshold Frequency Table** (Return intervals 1.5–500 years)"
                                )
                                st.dataframe(
                                    format_quantile_df(thr_res["quantile_df"]),
                                    use_container_width=True,
                                )

        st.divider()

    # ── Multi-gage FFA comparison ─────────────────────────────────────────
    if enable_ffa and len(gage_list) > 1:
        sites_with_ffa = [
            s
            for s in gage_list
            if s in st.session_state.ffa_results
            and not st.session_state.ffa_results[s].get("error")
        ]
        if sites_with_ffa:
            st.subheader("Flood Frequency Comparison")
            rows = []
            for sno in sites_with_ffa:
                r = st.session_state.ffa_results[sno]
                q100 = r["quantile_df"][r["quantile_df"]["Return Interval (yr)"] == 100][
                    "Flow (cfs)"
                ].values
                q100_val = int(q100[0]) if len(q100) > 0 else None
                attrs = st.session_state.gage_info.get(sno, {})
                _wy_cmp = st.session_state.ffa_year_range.get(sno)
                yr_range_str = f"{_wy_cmp[0]}–{_wy_cmp[1]}" if _wy_cmp else "Full record"
                rows.append(
                    {
                        "Site No": sno,
                        "Site Name": attrs.get("name", sno),
                        "FFA Year Range": yr_range_str,
                        "100-yr Flow (cfs)": f"{q100_val:,}" if q100_val else "N/A",
                        "Weighted Skew": f"{r['parameters']['skew_weighted']:.3f}",
                        "Method": r.get("method", "ema").upper(),
                        "Converged": "Yes" if r.get("converged") else "No",
                    }
                )
            st.dataframe(pd.DataFrame(rows), use_container_width=True)

    plt.close("all")


# ---------------------------------------------------------------------------
# Export section
# ---------------------------------------------------------------------------
if st.session_state.gage_data:
    st.sidebar.header("Export")

    selected_gages = st.sidebar.multiselect(
        "Select gages to export",
        options=list(st.session_state.gage_data.keys()),
        default=list(st.session_state.gage_data.keys()),
    )

    if selected_gages:
        zip_buffer = io.BytesIO()

        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
            for site_no in selected_gages:
                gage_info = st.session_state.gage_info[site_no]
                full_data = st.session_state.gage_data[site_no]
                site_figs = st.session_state.figures.get(site_no, {})

                # Daily flow CSV
                csv_buffer = io.StringIO()
                export_df = full_data.copy()
                export_df.index.name = "date"
                export_df.to_csv(csv_buffer)
                zf.writestr(f"{site_no}/daily_flow.csv", csv_buffer.getvalue())

                # Plots as PNG / stats as CSV
                for name, fig in site_figs.items():
                    if name.endswith("_stats"):
                        csv_buffer = io.StringIO()
                        fig.to_csv(csv_buffer, index=False)
                        zf.writestr(f"{site_no}/{name}.csv", csv_buffer.getvalue())
                    else:
                        img_buffer = io.BytesIO()
                        fig.savefig(img_buffer, format="png", dpi=300, bbox_inches="tight")
                        img_buffer.seek(0)
                        zf.writestr(f"{site_no}/{name}.png", img_buffer.read())

                if "summary_stats" not in site_figs and show_summary:
                    summary_stats = Hydrograph.get_summary_stats(full_data)
                    csv_buffer = io.StringIO()
                    summary_stats.to_csv(csv_buffer, index=False)
                    zf.writestr(f"{site_no}/summary_stats.csv", csv_buffer.getvalue())

                # FFA export
                if enable_ffa and site_no in st.session_state.ffa_results:
                    ffa_result = st.session_state.ffa_results[site_no]
                    if not ffa_result.get("error"):
                        selected_skew_labels = [
                            lbl
                            for lbl, on in [
                                ("Station Skew", skew_station_on),
                                ("Weighted Skew", skew_weighted_on),
                                ("Regional Skew", skew_regional_on),
                            ]
                            if on
                        ]
                        skew_curves_export = (
                            build_skew_curves_dict(ffa_result, selected_skew_labels) or None
                        )
                        freq_fig_for_export = None
                        if show_freq_curve and ffa_result.get("b17c"):
                            extra_curves_export = None
                            if show_threshold_curve:
                                thr_exp = st.session_state.ffa_results_with_thresholds.get(site_no)
                                if thr_exp and not thr_exp.get("error") and thr_exp.get("b17c"):
                                    r_texp = thr_exp["b17c"].results
                                    extra_curves_export = {
                                        "With Perception Thresholds": (
                                            r_texp.mean_log,
                                            r_texp.std_log,
                                            r_texp.skew_used,
                                            r_texp.n_systematic or r_texp.n_peaks,
                                        )
                                    }
                            freq_fig_for_export = plot_frequency_curve_streamlit(
                                ffa_result["b17c"],
                                site_name=gage_info.get("name", ""),
                                site_no=site_no,
                                skew_curves=skew_curves_export,
                                extra_curves=extra_curves_export,
                            )
                        export_ffa_to_zip(zf, site_no, ffa_result, freq_fig_for_export)

                        # Threshold-aware export
                        if show_threshold_curve:
                            thr_exp = st.session_state.ffa_results_with_thresholds.get(site_no)
                            if thr_exp and not thr_exp.get("error") and thr_exp.get("b17c"):
                                thr_skew_tables_exp = compute_skew_tables(
                                    thr_exp, selected_skew_labels
                                )
                                for lbl, tbl in thr_skew_tables_exp.items():
                                    lbl_clean = lbl.lower().replace(" ", "_")
                                    csv_buf = io.StringIO()
                                    format_quantile_df(tbl).to_csv(csv_buf, index=False)
                                    zf.writestr(
                                        f"{site_no}/freq_table_w_thresholds_{lbl_clean}.csv",
                                        csv_buf.getvalue(),
                                    )
                                csv_buf = io.StringIO()
                                format_parameters_df(thr_exp["parameters"]).to_csv(
                                    csv_buf, index=False
                                )
                                zf.writestr(
                                    f"{site_no}/lp3_parameters_with_thresholds.csv",
                                    csv_buf.getvalue(),
                                )

            # Multi-gage comparison CSV
            if enable_ffa and len(selected_gages) > 1:
                site_results_for_export = {}
                for sno in selected_gages:
                    if sno in st.session_state.ffa_results:
                        attrs = st.session_state.gage_info.get(sno, {})
                        site_results_for_export[sno] = {
                            "site_name": attrs.get("name", sno),
                            "drainage_area_sqmi": attrs.get("drainage_area", 0),
                            "ffa_results": st.session_state.ffa_results[sno],
                        }
                if site_results_for_export:
                    export_comparison_csv(zf, site_results_for_export)

        zip_buffer.seek(0)

        st.sidebar.download_button(
            label=f"Download {len(selected_gages)} Gage(s) as ZIP",
            data=zip_buffer,
            file_name=f"usgs_data_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zip",
            mime="application/zip",
        )

# Footer
st.sidebar.markdown("---")
st.sidebar.markdown("""
    **hydrolib** v0.1.0
    [GitHub](https://github.com/pinhead001/hydrolib)
    """)
