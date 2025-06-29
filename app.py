import streamlit as st
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point, box, LineString
from skyfield.api import load
import numpy as np
import os
from datetime import date as dt_date
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import shutil
import zipfile
import tempfile
from fpdf import FPDF
import io

# --- Helper functions ---

def download_tle(group, save_folder, max_days=1.0):
    os.makedirs(save_folder, exist_ok=True)
    filename = f"{group}.tle"
    filepath = os.path.join(save_folder, filename)
    url = f'https://celestrak.org/NORAD/elements/gp.php?GROUP={group}&FORMAT=TLE'

    # Simple caching by checking if file exists; no built-in load.exists or days_old in skyfield
    if not os.path.exists(filepath):
        import requests
        r = requests.get(url)
        with open(filepath, "w") as f:
            f.write(r.text)
        source = f"Downloaded fresh TLE from: {url}"
    else:
        source = f"Loaded cached TLE file: {filepath}"
    return filepath, source

def create_aoi(min_lon, min_lat, max_lon, max_lat):
    return gpd.GeoDataFrame({'geometry': [box(min_lon, min_lat, max_lon, max_lat)]}, crs="EPSG:4326")

def generate_times(year, month, day, interval_minutes):
    ts = load.timescale()
    minutes = np.arange(0, 24*60, interval_minutes)
    return ts.utc(year, month, day, 0, minutes)

def find_passing_sats(satellites, times, aoi, swath_width_m):
    passing_sats = []
    colors = plt.cm.get_cmap('tab20', len(satellites))
    buffer_radius = swath_width_m / 2
    plot_data = []

    for i, sat in enumerate(satellites):
        subpoint = sat.at(times).subpoint()
        lats = subpoint.latitude.degrees
        lons = subpoint.longitude.degrees
        points = [Point(lon, lat) for lon, lat in zip(lons, lats)]
        sat_gdf = gpd.GeoDataFrame({'geometry': points}, crs="EPSG:4326")
        sat_gdf_3857 = sat_gdf.to_crs("EPSG:3857")
        swath_buffers = sat_gdf_3857.buffer(buffer_radius)
        swath_gdf = gpd.GeoDataFrame(geometry=swath_buffers, crs="EPSG:3857").to_crs("EPSG:4326")

        intersect = gpd.sjoin(swath_gdf, aoi, how='inner', predicate='intersects')

        if not intersect.empty:
            passing_time = times[intersect.index[0]].utc_iso()
            passing_sats.append((sat.name, passing_time))
            trace_line = LineString([(pt.x, pt.y) for pt in points])
            plot_data.append({
                'swath_gdf': swath_gdf,
                'trace_line': trace_line,
                'color': colors(i),
                'name': sat.name
            })

    return passing_sats, plot_data

def plot_results(aoi, plot_data, swath_width_km):
    fig = plt.figure(figsize=(10, 7))
    ax = plt.axes(projection=ccrs.PlateCarree())

    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, edgecolor='black', alpha=0.3)
    ax.add_feature(cfeature.OCEAN, alpha=0.1)

    aoi.boundary.plot(ax=ax, color='black', linewidth=2, label='AOI')

    for item in plot_data:
        item['swath_gdf'].plot(ax=ax, color='red', alpha=0.4)
        ax.plot(*item['trace_line'].xy, color=item['color'], linewidth=1.5, label=item['name'])

    bounds = aoi.total_bounds
    ax.set_extent([bounds[0] - 2, bounds[2] + 2, bounds[1] - 2, bounds[3] + 2], crs=ccrs.PlateCarree())

    ax.set_title("Satellite Passes", fontsize=14)
    ax.legend(fontsize=8, loc='lower left')
    plt.tight_layout()
    return fig

def create_pdf_report_text_and_image(sat_type, year, month, day, swath_km, tle_source, passing_sats, fig):
    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()

    pdf.set_font("Arial", "B", 16)
    pdf.cell(0, 10, "Satellite Pass Daily Report", ln=True, align='C')
    pdf.ln(10)
    pdf.set_font("Arial", "", 12)
    pdf.multi_cell(0, 10, f"TLE Source: {tle_source}")
    pdf.multi_cell(0, 10, f"Satellite Group: {sat_type}")
    pdf.multi_cell(0, 10, f"Date: {year}-{month:02d}-{day:02d}")
    pdf.multi_cell(0, 10, f"Swath Width: {swath_km} km")
    pdf.ln(5)

    if passing_sats:
        pdf.set_font("Arial", "B", 12)
        pdf.cell(0, 10, f"{len(passing_sats)} satellite(s) passed over the AOI:", ln=True)
        pdf.set_font("Arial", "", 12)
        for name, t in passing_sats:
            pdf.multi_cell(0, 10, f"{name} at {t}")
    else:
        pdf.multi_cell(0, 10, "No satellites passed over the area.")
    pdf.ln(10)

    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmpfile:
        temp_img_path = tmpfile.name
        fig.savefig(temp_img_path, format='PNG', dpi=300)

    page_width = pdf.w - 2*pdf.l_margin
    pdf.image(temp_img_path, x=pdf.l_margin, w=page_width)
    os.remove(temp_img_path)

    pdf_output = pdf.output(dest='S').encode('latin1')
    pdf_bytes = io.BytesIO(pdf_output)
    pdf_bytes.seek(0)
    return pdf_bytes

def show_data_links(sat_name):
    search_url = f"https://www.google.com/search?q={sat_name.replace(' ', '+')}+satellite+data+download"
    st.markdown(f"🔎 [Search for data from {sat_name}]({search_url})")

# --- Streamlit app ---

st.title("Satellite Pass Finder")
st.markdown("This tool finds satellites passing over your AOI and helps you find their data.")

# Initialize session state flags and storage
if "run_analysis" not in st.session_state:
    st.session_state.run_analysis = False
if "passing_sats" not in st.session_state:
    st.session_state.passing_sats = []
if "plot_data" not in st.session_state:
    st.session_state.plot_data = []
if "fig" not in st.session_state:
    st.session_state.fig = None
if "tle_source" not in st.session_state:
    st.session_state.tle_source = ""

# AOI inputs
st.header("1. Define Search Area")
min_lon = st.number_input("Min Longitude", value=127.5)
min_lat = st.number_input("Min Latitude", value=25.5)
max_lon = st.number_input("Max Longitude", value=129.0)
max_lat = st.number_input("Max Latitude", value=27.0)

aoi = create_aoi(min_lon, min_lat, max_lon, max_lat)

# Satellite group
sat_type = st.selectbox("Satellite Group", options=["active", "resource", "science"])

# Date and parameters
date = st.date_input("Date", value=dt_date.today())
swath_km = st.slider("Swath Width (km)", min_value=10, max_value=100, value=30)
interval = st.slider("Time Interval (minutes)", min_value=1, max_value=60, value=10)

# Button to run analysis
if st.button("Run Analysis"):
    st.session_state.run_analysis = True

# Run analysis when triggered
if st.session_state.run_analysis:
    tle_folder = "./.cache_tle"
    tle_path, tle_source = download_tle(sat_type, tle_folder)

    satellites = load.tle_file(tle_path)
    times = generate_times(date.year, date.month, date.day, interval)
    swath_width_m = swath_km * 1000

    passing_sats, plot_data = find_passing_sats(satellites, times, aoi, swath_width_m)
    fig = plot_results(aoi, plot_data, swath_km)

    st.session_state.passing_sats = passing_sats
    st.session_state.plot_data = plot_data
    st.session_state.fig = fig
    st.session_state.tle_source = tle_source

    st.subheader("Results")
    st.write(tle_source)
    if passing_sats:
        st.success(f"{len(passing_sats)} satellite(s) passed over the AOI.")
        for name, t in passing_sats:
            st.write(f"🛰️ {name} at {t}")
            show_data_links(name)
    else:
        st.warning("No satellites passed over the AOI.")
    st.pyplot(fig)

    # Generate and provide PDF report download
    pdf_buffer = create_pdf_report_text_and_image(sat_type, date.year, date.month, date.day, swath_km, tle_source, passing_sats, fig)
    st.download_button(
        label="📄 Download Daily Report (.pdf)",
        data=pdf_buffer,
        file_name="satellite_daily_report.pdf",
        mime="application/pdf"
    )

    st.session_state.run_analysis = False

# Show last results if available and not currently running analysis
elif st.session_state.fig:
    st.subheader("Last Analysis Results")
    st.write(st.session_state.tle_source)
    if st.session_state.passing_sats:
        st.success(f"{len(st.session_state.passing_sats)} satellite(s) passed over the AOI.")
        for name, t in st.session_state.passing_sats:
            st.write(f"🛰️ {name} at {t}")
            show_data_links(name)
    else:
        st.warning("No satellites passed over the AOI.")
    st.pyplot(st.session_state.fig)

st.markdown(
    """
    ---
    📢 **Disclaimer**

    This application retrieves satellite orbital data exclusively from [CelesTrak](https://celestrak.org/), 
    a public source of satellite TLE (Two-Line Element) data. The accuracy of pass predictions depends on 
    the quality and update frequency of CelesTrak's datasets. Only satellites listed in the selected 
    CelesTrak group will be included in the analysis.

    This tool does **not** query all satellites in orbit — only those published and maintained by CelesTrak.
    """
)
