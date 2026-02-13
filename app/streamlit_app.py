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
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st

matplotlib.use("Agg")  # Use non-interactive backend

# Add parent directory to path for hydrolib import (needed for Streamlit Cloud)
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import hydrolib
from hydrolib import Hydrograph
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

# Sidebar inputs
st.sidebar.header("Input Parameters")

# Gage input - single or multiple
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
    # Parse input
    gage_list = [
        g.strip() for g in gage_input.replace(",", "\n").split("\n") if g.strip()
    ]

st.sidebar.markdown(f"**{len(gage_list)} gage(s) selected**")

# Plot options
st.sidebar.header("Plot Options")
show_timeseries = st.sidebar.checkbox("Daily Time Series", value=True)
show_summary = st.sidebar.checkbox("Summary Hydrograph", value=True)
show_fdc = st.sidebar.checkbox("Flow Duration Curve", value=True)

# Download button
download_data = st.sidebar.button("Download Data", type="primary")

# Initialize session state
if "gage_data" not in st.session_state:
    st.session_state.gage_data = {}  # Stores full daily data for each gage
if "gage_info" not in st.session_state:
    st.session_state.gage_info = {}  # Stores site info for each gage
if "figures" not in st.session_state:
    st.session_state.figures = {}
if "expanded_gage" not in st.session_state:
    st.session_state.expanded_gage = None
if "expanded_plot_idx" not in st.session_state:
    st.session_state.expanded_plot_idx = 0


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
    }
    return names.get(plot_key, plot_key)


def format_date(dt):
    """Format datetime as M/D/YYYY string."""
    return f"{dt.month}/{dt.day}/{dt.year}"


def generate_plots(site_no, plot_data, gage_info, por_start_str=None, por_end_str=None):
    """Generate plots for a gage using the provided data."""
    site_figs = {}

    # Daily time series
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

    # Summary hydrograph
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
        # Store summary stats for export
        site_figs["summary_stats"] = Hydrograph.get_summary_stats(plot_data)

    # Flow duration curve
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


# Download data when button clicked
if download_data and gage_list:
    st.session_state.gage_data = {}
    st.session_state.gage_info = {}
    st.session_state.figures = {}

    progress_bar = st.progress(0)
    status_text = st.empty()

    for idx, site_no in enumerate(gage_list):
        status_text.text(f"Fetching site info for {site_no}...")

        try:
            # Create gage object and fetch site info (including POR)
            gage = USGSgage(site_no)
            gage.fetch_site_info()

            status_text.text(f"Downloading data for {site_no}...")

            # Download entire period of record using POR dates from site info
            # USGS API returns limited data without date range, so use POR dates
            start_dt = gage.daily_por_start if gage.daily_por_start else "1900-01-01"
            end_dt = gage.daily_por_end if gage.daily_por_end else None
            daily_data = gage.download_daily_flow(start_date=start_dt, end_date=end_dt)

            # Store gage info
            site_name_display = gage.site_name or "Unknown"
            por_start_date = daily_data.index.min()
            por_end_date = daily_data.index.max()

            st.session_state.gage_info[site_no] = {
                "name": site_name_display,
                "drainage_area": gage.drainage_area,
                "por_start": por_start_date,
                "por_end": por_end_date,
                "records": len(daily_data),
            }

            # Store full daily data
            st.session_state.gage_data[site_no] = daily_data

        except Exception as e:
            st.error(f"Error processing {site_no}: {str(e)}")

        progress_bar.progress((idx + 1) / len(gage_list))

    status_text.text("Download complete!")
    plt.close("all")
    st.rerun()  # Refresh page to display plots with new data


# Show plot date range controls only if data is loaded
if st.session_state.gage_data:
    st.sidebar.header("Plot Date Range")
    st.sidebar.markdown("*Filter plots without re-downloading*")

    # Get overall date range from all loaded gages
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

    # Update plots button
    if st.sidebar.button("Update Plots", type="secondary"):
        st.session_state.figures = {}

    # Generate/display plots for each gage
    for site_no in st.session_state.gage_data.keys():
        gage_info = st.session_state.gage_info[site_no]
        full_data = st.session_state.gage_data[site_no]

        # Filter data for plot date range
        plot_data = full_data[
            (full_data.index >= pd.Timestamp(plot_start))
            & (full_data.index <= pd.Timestamp(plot_end))
        ]

        if len(plot_data) == 0:
            st.warning(f"No data for {site_no} in selected date range")
            continue

        # Display gage header
        st.subheader(f"USGS {site_no} - {gage_info['name']}")

        # Format POR and plot dates
        por_str = f"{format_date(gage_info['por_start'])} - {format_date(gage_info['por_end'])}"
        plot_str = f"{format_date(plot_data.index.min())} - {format_date(plot_data.index.max())}"
        da_str = f"{gage_info['drainage_area']:,.1f} sq mi" if gage_info['drainage_area'] else "N/A"

        # Display info
        info_cols = st.columns(4)
        info_cols[0].markdown(f"**Drainage Area**<br><small>{da_str}</small>", unsafe_allow_html=True)
        info_cols[1].markdown(f"**POR**<br><small>{por_str}</small>", unsafe_allow_html=True)
        info_cols[2].markdown(f"**Plot Range**<br><small>{plot_str}</small>", unsafe_allow_html=True)
        info_cols[3].markdown(f"**Records**<br><small>{len(plot_data):,} days</small>", unsafe_allow_html=True)

        # Generate plots if not cached or need refresh
        if site_no not in st.session_state.figures:
            # Determine if we need to show POR annotation (when plot range differs from POR)
            por_start_str = None
            por_end_str = None
            if (plot_data.index.min() > gage_info['por_start'] or
                plot_data.index.max() < gage_info['por_end']):
                por_start_str = format_date(gage_info['por_start'])
                por_end_str = format_date(gage_info['por_end'])

            st.session_state.figures[site_no] = generate_plots(
                site_no, plot_data, gage_info, por_start_str, por_end_str
            )

        site_figs = st.session_state.figures[site_no]

        # Determine columns based on selected plots
        plot_keys = [k for k in site_figs.keys() if not k.endswith("_stats")]
        if plot_keys:
            cols = st.columns(len(plot_keys))
            for i, plot_key in enumerate(plot_keys):
                with cols[i]:
                    st.pyplot(site_figs[plot_key])

        st.divider()

    plt.close("all")


# Export section
if st.session_state.gage_data:
    st.sidebar.header("Export")

    # Select which gages to export
    selected_gages = st.sidebar.multiselect(
        "Select gages to export",
        options=list(st.session_state.gage_data.keys()),
        default=list(st.session_state.gage_data.keys()),
    )

    if selected_gages:
        # Create ZIP file in memory
        zip_buffer = io.BytesIO()

        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
            for site_no in selected_gages:
                gage_info = st.session_state.gage_info[site_no]
                full_data = st.session_state.gage_data[site_no]
                site_figs = st.session_state.figures.get(site_no, {})

                # Export daily flow data CSV
                csv_buffer = io.StringIO()
                export_df = full_data.copy()
                export_df.index.name = "date"
                export_df.to_csv(csv_buffer)
                zf.writestr(f"{site_no}/daily_flow.csv", csv_buffer.getvalue())

                # Export plots as PNG
                for name, fig in site_figs.items():
                    if name.endswith("_stats"):
                        # Export stats as CSV
                        csv_buffer = io.StringIO()
                        fig.to_csv(csv_buffer, index=False)
                        zf.writestr(f"{site_no}/{name}.csv", csv_buffer.getvalue())
                    else:
                        # Save PNG
                        img_buffer = io.BytesIO()
                        fig.savefig(
                            img_buffer, format="png", dpi=300, bbox_inches="tight"
                        )
                        img_buffer.seek(0)
                        zf.writestr(f"{site_no}/{name}.png", img_buffer.read())

                # Export summary stats if not already in figures
                if "summary_stats" not in site_figs and show_summary:
                    summary_stats = Hydrograph.get_summary_stats(full_data)
                    csv_buffer = io.StringIO()
                    summary_stats.to_csv(csv_buffer, index=False)
                    zf.writestr(f"{site_no}/summary_stats.csv", csv_buffer.getvalue())

        zip_buffer.seek(0)

        st.sidebar.download_button(
            label=f"Download {len(selected_gages)} Gage(s) as ZIP",
            data=zip_buffer,
            file_name=f"usgs_data_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zip",
            mime="application/zip",
        )

# Footer
st.sidebar.markdown("---")
st.sidebar.markdown(
    """
    **hydrolib** v0.0.3
    [GitHub](https://github.com/pinhead001/hydrolib)
    """
)
