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
import streamlit as st

matplotlib.use("Agg")  # Use non-interactive backend

# Add parent directory to path for hydrolib import (needed for Streamlit Cloud)
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import hydrolib
from hydrolib import Hydrograph
from hydrolib.usgs import USGSgage

st.set_page_config(
    page_title="USGS Hydrograph Generator",
    page_icon=":ocean:",
    layout="wide",
)

# Custom CSS for clickable plots and modal
st.markdown(
    """
    <style>
    .plot-container {
        cursor: pointer;
        transition: transform 0.2s, box-shadow 0.2s;
        border-radius: 8px;
        padding: 5px;
    }
    .plot-container:hover {
        transform: scale(1.02);
        box-shadow: 0 4px 12px rgba(0,0,0,0.15);
    }
    .modal-nav-btn {
        font-size: 24px;
        padding: 10px 20px;
    }
    .plot-indicator {
        text-align: center;
        font-size: 14px;
        color: #666;
        margin-top: 10px;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# JavaScript for arrow key navigation
ARROW_KEY_JS = """
<script>
document.addEventListener('keydown', function(e) {
    if (e.key === 'ArrowLeft') {
        // Find and click the Previous button
        const prevBtn = document.querySelector('[data-testid="baseButton-secondary"]:first-of-type');
        if (prevBtn && prevBtn.innerText.includes('Previous')) {
            prevBtn.click();
        }
    } else if (e.key === 'ArrowRight') {
        // Find and click the Next button
        const nextBtn = document.querySelector('[data-testid="baseButton-secondary"]:last-of-type');
        if (nextBtn && nextBtn.innerText.includes('Next')) {
            nextBtn.click();
        }
    } else if (e.key === 'Escape') {
        // Find and click the Close button
        const closeBtn = document.querySelector('[data-testid="baseButton-primary"]');
        if (closeBtn && closeBtn.innerText.includes('Close')) {
            closeBtn.click();
        }
    }
});
</script>
"""

st.title("USGS Daily Flow Hydrograph Generator")
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

# Date range
st.sidebar.header("Date Range")
col1, col2 = st.sidebar.columns(2)
start_date = col1.date_input("Start Date", value=datetime(1970, 10, 1))
end_date = col2.date_input("End Date", value=datetime(2025, 10, 1))

# Plot options
st.sidebar.header("Plot Options")
show_timeseries = st.sidebar.checkbox("Daily Time Series", value=True)
show_summary = st.sidebar.checkbox("Summary Hydrograph", value=True)
show_fdc = st.sidebar.checkbox("Flow Duration Curve", value=True)

# Run button
run_analysis = st.sidebar.button("Generate Plots", type="primary")

# Initialize session state
if "figures" not in st.session_state:
    st.session_state.figures = {}
if "gage_info" not in st.session_state:
    st.session_state.gage_info = {}
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
        if name != "fdc_stats"
    ]


def get_plot_display_name(plot_key):
    """Convert plot key to display name."""
    names = {
        "daily_timeseries": "Daily Time Series",
        "summary_hydrograph": "Summary Hydrograph",
        "flow_duration_curve": "Flow Duration Curve",
    }
    return names.get(plot_key, plot_key)


def show_expanded_plot():
    """Display the expanded plot modal with navigation."""
    site_no = st.session_state.expanded_gage
    if not site_no or site_no not in st.session_state.figures:
        return

    plot_list = get_plot_list(site_no)
    if not plot_list:
        return

    # Ensure index is valid
    idx = st.session_state.expanded_plot_idx % len(plot_list)
    current_plot = plot_list[idx]

    gage_info = st.session_state.gage_info.get(site_no, {})
    gage_name = gage_info.get("name", "")

    # Modal header
    st.markdown("---")
    st.subheader(f"USGS {site_no} - {gage_name}")
    st.markdown(
        f"**{get_plot_display_name(current_plot)}** ({idx + 1} of {len(plot_list)})"
    )

    # Navigation buttons
    nav_cols = st.columns([1, 3, 1])

    with nav_cols[0]:
        if st.button("< Previous", key="prev_btn", use_container_width=True):
            st.session_state.expanded_plot_idx = (idx - 1) % len(plot_list)
            st.rerun()

    with nav_cols[2]:
        if st.button("Next >", key="next_btn", use_container_width=True):
            st.session_state.expanded_plot_idx = (idx + 1) % len(plot_list)
            st.rerun()

    # Display the expanded plot
    fig = st.session_state.figures[site_no][current_plot]
    st.pyplot(fig, use_container_width=True)

    # Show flow duration stats if viewing FDC
    if current_plot == "flow_duration_curve" and "fdc_stats" in st.session_state.figures[site_no]:
        with st.expander("Flow Duration Statistics"):
            st.dataframe(st.session_state.figures[site_no]["fdc_stats"])

    # Close button
    close_cols = st.columns([2, 1, 2])
    with close_cols[1]:
        if st.button("Close", key="close_btn", type="primary", use_container_width=True):
            st.session_state.expanded_gage = None
            st.session_state.expanded_plot_idx = 0
            st.rerun()

    # Plot indicator dots
    dot_str = " ".join(
        ["●" if i == idx else "○" for i in range(len(plot_list))]
    )
    st.markdown(
        f'<div class="plot-indicator">{dot_str}<br><small>Use arrow keys or buttons to navigate, Esc to close</small></div>',
        unsafe_allow_html=True,
    )

    # Inject arrow key JavaScript
    st.components.v1.html(ARROW_KEY_JS, height=0)


# Show expanded view if a plot is selected
if st.session_state.expanded_gage:
    show_expanded_plot()
else:
    # Normal view - generate and display plots

    if run_analysis and gage_list:
        st.session_state.figures = {}
        st.session_state.gage_info = {}

        progress_bar = st.progress(0)
        status_text = st.empty()

        for idx, site_no in enumerate(gage_list):
            status_text.text(f"Processing {site_no}...")

            try:
                # Download data
                gage = USGSgage(site_no)
                daily_data = gage.download_daily_flow(
                    start_date=start_date.strftime("%Y-%m-%d"),
                    end_date=end_date.strftime("%Y-%m-%d"),
                )

                st.subheader(f"USGS {site_no} - {gage.site_name}")

                # Store gage info
                st.session_state.gage_info[site_no] = {
                    "name": gage.site_name,
                    "drainage_area": gage.drainage_area,
                    "records": len(daily_data),
                }

                info_cols = st.columns(3)
                info_cols[0].metric(
                    "Drainage Area", f"{gage.drainage_area or 'N/A'} sq mi"
                )
                info_cols[1].metric("Records", f"{len(daily_data):,} days")
                info_cols[2].metric(
                    "Period",
                    f"{daily_data.index.min().strftime('%Y-%m-%d')} to {daily_data.index.max().strftime('%Y-%m-%d')}",
                )

                # Determine columns based on selected plots
                plots_to_show = []
                if show_timeseries:
                    plots_to_show.append("timeseries")
                if show_summary:
                    plots_to_show.append("summary")
                if show_fdc:
                    plots_to_show.append("fdc")

                if len(plots_to_show) > 0:
                    cols = st.columns(len(plots_to_show))

                site_figs = {}

                # Daily time series
                if show_timeseries:
                    fig1 = Hydrograph.plot_daily_timeseries(
                        daily_data,
                        site_name=gage.site_name,
                        site_no=gage.site_no,
                        figsize=(10, 4),
                    )
                    col_idx = plots_to_show.index("timeseries")
                    with cols[col_idx]:
                        st.pyplot(fig1)
                        if st.button(
                            "Expand",
                            key=f"expand_ts_{site_no}",
                            use_container_width=True,
                        ):
                            st.session_state.expanded_gage = site_no
                            st.session_state.expanded_plot_idx = 0
                            st.rerun()
                    site_figs["daily_timeseries"] = fig1

                # Summary hydrograph
                if show_summary:
                    fig2 = Hydrograph.plot_summary_hydrograph(
                        daily_data,
                        site_name=gage.site_name,
                        site_no=gage.site_no,
                        figsize=(10, 4),
                        percentiles=[10, 25, 50, 75, 90],
                    )
                    col_idx = plots_to_show.index("summary")
                    with cols[col_idx]:
                        st.pyplot(fig2)
                        if st.button(
                            "Expand",
                            key=f"expand_sh_{site_no}",
                            use_container_width=True,
                        ):
                            st.session_state.expanded_gage = site_no
                            # Find the index of summary_hydrograph in plot list
                            st.session_state.expanded_plot_idx = (
                                1 if show_timeseries else 0
                            )
                            st.rerun()
                    site_figs["summary_hydrograph"] = fig2

                # Flow duration curve
                if show_fdc:
                    fig3, stats_df = Hydrograph.plot_flow_duration_curve(
                        daily_data,
                        site_name=gage.site_name,
                        site_no=gage.site_no,
                        figsize=(10, 4),
                    )
                    col_idx = plots_to_show.index("fdc")
                    with cols[col_idx]:
                        st.pyplot(fig3)
                        if st.button(
                            "Expand",
                            key=f"expand_fdc_{site_no}",
                            use_container_width=True,
                        ):
                            st.session_state.expanded_gage = site_no
                            # Find the index of flow_duration_curve in plot list
                            fdc_idx = 0
                            if show_timeseries:
                                fdc_idx += 1
                            if show_summary:
                                fdc_idx += 1
                            st.session_state.expanded_plot_idx = fdc_idx
                            st.rerun()
                    site_figs["flow_duration_curve"] = fig3
                    site_figs["fdc_stats"] = stats_df

                st.session_state.figures[site_no] = site_figs
                st.divider()

            except Exception as e:
                st.error(f"Error processing {site_no}: {str(e)}")

            progress_bar.progress((idx + 1) / len(gage_list))

        status_text.text("Complete!")
        plt.close("all")

    # Show existing figures if available (after initial run)
    elif st.session_state.figures:
        for site_no, site_figs in st.session_state.figures.items():
            gage_info = st.session_state.gage_info.get(site_no, {})

            st.subheader(f"USGS {site_no} - {gage_info.get('name', '')}")

            info_cols = st.columns(3)
            info_cols[0].metric(
                "Drainage Area", f"{gage_info.get('drainage_area', 'N/A')} sq mi"
            )
            info_cols[1].metric("Records", f"{gage_info.get('records', 0):,} days")

            # Display plots with expand buttons
            plot_keys = [k for k in site_figs.keys() if k != "fdc_stats"]
            if plot_keys:
                cols = st.columns(len(plot_keys))
                for i, plot_key in enumerate(plot_keys):
                    with cols[i]:
                        st.pyplot(site_figs[plot_key])
                        if st.button(
                            "Expand",
                            key=f"expand_{plot_key}_{site_no}_existing",
                            use_container_width=True,
                        ):
                            st.session_state.expanded_gage = site_no
                            st.session_state.expanded_plot_idx = i
                            st.rerun()

            st.divider()

# Download section
if st.session_state.figures and not st.session_state.expanded_gage:
    st.sidebar.header("Download")

    # Select which gages to download
    selected_gages = st.sidebar.multiselect(
        "Select gages to download",
        options=list(st.session_state.figures.keys()),
        default=list(st.session_state.figures.keys()),
    )

    if selected_gages:
        # Create ZIP file in memory
        zip_buffer = io.BytesIO()

        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
            for site_no in selected_gages:
                site_figs = st.session_state.figures[site_no]

                for name, fig in site_figs.items():
                    if name == "fdc_stats":
                        # Save CSV
                        csv_buffer = io.StringIO()
                        fig.to_csv(csv_buffer)
                        zf.writestr(f"{site_no}/{name}.csv", csv_buffer.getvalue())
                    else:
                        # Save PNG
                        img_buffer = io.BytesIO()
                        fig.savefig(
                            img_buffer, format="png", dpi=300, bbox_inches="tight"
                        )
                        img_buffer.seek(0)
                        zf.writestr(f"{site_no}/{name}.png", img_buffer.read())

        zip_buffer.seek(0)

        st.sidebar.download_button(
            label=f"Download {len(selected_gages)} Gage(s) as ZIP",
            data=zip_buffer,
            file_name=f"usgs_plots_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zip",
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
