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
import pandas as pd
import zipfile
import tempfile
from fpdf import FPDF
from PIL import Image
import io

if "run_analysis_trigger" in st.session_state and st.session_state["run_analysis_trigger"]:
    st.session_state["run_analysis_trigger"] = False

# --- Helper Functions ---

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

def download_tle(group, save_folder, max_days=1.0):
    os.makedirs(save_folder, exist_ok=True)
    filename = f"{group}.tle"
    filepath = os.path.join(save_folder, filename)
    base_url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP={group}&FORMAT=TLE'
    url = base_url.format(group=group)

    if not load.exists(filepath) or load.days_old(filepath) >= max_days:
        load.download(url, filename=filepath)
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

def preview_aoi_map(aoi):
    fig, ax = plt.subplots(figsize=(6, 6))
    aoi.boundary.plot(ax=ax, edgecolor='blue', linewidth=2)
    ax.set_title("AOI Preview")
    bounds = aoi.total_bounds
    ax.set_xlim(bounds[0] - 1, bounds[2] + 1)
    ax.set_ylim(bounds[1] - 1, bounds[3] + 1)
    plt.tight_layout()
    return fig

# --- Streamlit app starts here ---

st.title("Satellite Pass Finder")
st.markdown("This tool finds satellites that pass over your AOI and tells you how it got the data.")

# Tabs for layout
tabs = st.tabs(["Search & Parameters", "Results", "Favorites"])

# --- Favorites / Presets management ---
if "favorites" not in st.session_state:
    st.session_state["favorites"] = {}

def save_favorite(name, data):
    st.session_state["favorites"][name] = data

def delete_favorite(name):
    if name in st.session_state["favorites"]:
        del st.session_state["favorites"][name]

with tabs[2]:
    st.header("Manage Favorites / Presets")

    if st.session_state["favorites"]:
        fav_to_load = st.selectbox("Load Favorite Preset", options=list(st.session_state["favorites"].keys()))
        if fav_to_load:
            if st.button("Load Selected Favorite"):
                fav = st.session_state["favorites"][fav_to_load]
                st.session_state["loaded_fav"] = fav
                st.experimental_rerun()

        st.subheader("Saved Favorites")
        for fav_name in list(st.session_state["favorites"].keys()):
            col1, col2 = st.columns([4, 1])
            with col1:
                st.write(fav_name)
            with col2:
                if st.button(f"Delete {fav_name}"):
                    delete_favorite(fav_name)
                    st.experimental_rerun()
    else:
        st.info("No favorites saved yet.")

    st.subheader("Save Current Settings as Favorite")
    fav_name_input = st.text_input("Favorite Name")

    if st.button("Save Favorite"):
        if not fav_name_input:
            st.warning("Please enter a favorite name.")
        else:
            # Save current inputs from session_state or default values
            data = {}
            if "loaded_fav" in st.session_state:
                del st.session_state["loaded_fav"]  # Clear loaded favorite so inputs are not overridden

            # Save AOI bounds & satellite group and parameters
            # If user uploaded a shapefile, store its bounds, else the manual input values
            if "uploaded_shp_bounds" in st.session_state:
                data["aoi_bounds"] = st.session_state["uploaded_shp_bounds"]
            else:
                data["aoi_bounds"] = {
                    "min_lon": st.session_state.get("min_lon", 127.5),
                    "min_lat": st.session_state.get("min_lat", 25.5),
                    "max_lon": st.session_state.get("max_lon", 129.0),
                    "max_lat": st.session_state.get("max_lat", 27.0),
                }
            data["sat_type"] = st.session_state.get("sat_type", "Earth Observation")
            data["date"] = st.session_state.get("date", dt_date.today())
            data["swath_km"] = st.session_state.get("swath_km", 30)
            data["interval"] = st.session_state.get("interval", 10)

            save_favorite(fav_name_input, data)
            st.success(f"Favorite '{fav_name_input}' saved!")

# --- Main Search & Parameters Tab ---
with tabs[0]:
    st.header("1. Define Search Area")
    uploaded_shp = st.file_uploader("Upload AOI Shapefile (.zip)", type=["zip"])

    if uploaded_shp is not None:
        with tempfile.TemporaryDirectory() as tmpdir:
            zip_path = os.path.join(tmpdir, "uploaded.zip")
            with open(zip_path, "wb") as f:
                f.write(uploaded_shp.read())
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(tmpdir)
            for file in os.listdir(tmpdir):
                if file.endswith(".shp"):
                    shp_path = os.path.join(tmpdir, file)
                    aoi = gpd.read_file(shp_path).to_crs("EPSG:4326")
                    bounds = aoi.total_bounds
                    min_lon, min_lat, max_lon, max_lat = bounds[0], bounds[1], bounds[2], bounds[3]
                    st.session_state["uploaded_shp_bounds"] = {
                        "min_lon": min_lon, "min_lat": min_lat, "max_lon": max_lon, "max_lat": max_lat
                    }
                    break
    else:
        # If a favorite is loaded, override input values with it
        if "loaded_fav" in st.session_state:
            fav = st.session_state["loaded_fav"]
            min_lon = fav["aoi_bounds"]["min_lon"]
            min_lat = fav["aoi_bounds"]["min_lat"]
            max_lon = fav["aoi_bounds"]["max_lon"]
            max_lat = fav["aoi_bounds"]["max_lat"]
            sat_type = fav["sat_type"]
            date = fav["date"]
            swath_km = fav["swath_km"]
            interval = fav["interval"]

            # Set inputs with loaded favorite values, store in session state
            st.session_state["min_lon"] = min_lon
            st.session_state["min_lat"] = min_lat
            st.session_state["max_lon"] = max_lon
            st.session_state["max_lat"] = max_lat
            st.session_state["sat_type"] = sat_type
            st.session_state["date"] = date
            st.session_state["swath_km"] = swath_km
            st.session_state["interval"] = interval
            del st.session_state["loaded_fav"]  # Clear so it doesn't reload repeatedly

        min_lon = st.number_input("Min Longitude", value=st.session_state.get("min_lon", 127.5))
        min_lat = st.number_input("Min Latitude", value=st.session_state.get("min_lat", 25.5))
        max_lon = st.number_input("Max Longitude", value=st.session_state.get("max_lon", 129.0))
        max_lat = st.number_input("Max Latitude", value=st.session_state.get("max_lat", 27.0))

        aoi = create_aoi(min_lon, min_lat, max_lon, max_lat)

    with st.expander("🗺️ Preview AOI Bounding Box"):
        aoi_fig = preview_aoi_map(aoi)
        st.pyplot(aoi_fig)

    st.header("2. Select Satellite Group and Parameters")

    group_options = {
        "Active": "active",
        "Earth Observation": "resource",
        "Scientific": "science",
        "CubeSats": "cubesat",
        "Weather": "weather",
        "GOES": "goes",
        "NOAA": "noaa",
        "Planet": "planet",
        "Last 30 Days": "last-30-days"
    }
    group_descriptions = {
        "Earth Observation": "Satellites used for imaging, environmental monitoring, and Earth resource mapping.",
        "Weather": "Satellites that provide meteorological data and atmospheric monitoring.",
        "CubeSats": "Miniaturized satellites often used for scientific and academic purposes.",
        "Scientific": "idk science stuff probably.",
        "NOAA": "National Oceanic and Atmospheric Administration satellites, mainly used for weather and ocean monitoring.",
        "GOES": "Geostationary Operational Environmental Satellites for continuous weather observation over the Americas.",
        "Planet": "planet.com",
        "Active": "All currently operational satellites tracked by CelesTrak."
    }

    sat_type = st.selectbox("Satellite Group", options=list(group_options.keys()), index=list(group_options.keys()).index(st.session_state.get("sat_type", "Earth Observation")))
    st.caption(f"📘 **Description:** {group_descriptions.get(sat_type, 'No description available.')}")

    date = st.date_input("Date", value=st.session_state.get("date", dt_date.today()))
    swath_km = st.slider("Swath Width (km)", min_value=10, max_value=100, value=st.session_state.get("swath_km", 30))
    interval = st.slider("Time Interval (minutes)", min_value=1, max_value=60, value=st.session_state.get("interval", 10))

    # Store inputs in session_state for favorites saving
    st.session_state["min_lon"] = min_lon
    st.session_state["min_lat"] = min_lat
    st.session_state["max_lon"] = max_lon
    st.session_state["max_lat"] = max_lat
    st.session_state["sat_type"] = sat_type
    st.session_state["date"] = date
    st.session_state["swath_km"] = swath_km
    st.session_state["interval"] = interval

    if st.button("Run Analysis"):
        with st.spinner("Running analysis..."):
            year, month, day = date.year, date.month, date.day
            tle_folder = "./.cache_tle"
            tle_path, tle_source = download_tle(group_options[sat_type], tle_folder)

            satellites = load.tle_file(tle_path)
            times = generate_times(year, month, day, interval)
            swath_width_m = swath_km * 1000

            passing_sats, plot_data = find_passing_sats(satellites, times, aoi, swath_width_m)
            fig = plot_results(aoi, plot_data, swath_km)

            # Save last run summary in session state for Results tab
            st.session_state["last_run"] = {
                "count": len(passing_sats),
                "sats": passing_sats,
                "tle_source": tle_source,
                "fig": fig,
                "passing_sats": passing_sats,
                "plot_data": plot_data,
                "params": {
                    "sat_type": sat_type,
                    "year": year,
                    "month": month,
                    "day": day,
                    "swath_km": swath_km,
                },
                "aoi": aoi
            }
            st.session_state["run_analysis_trigger"] = True
            st.experimental_rerun()

# --- Results Tab ---
with tabs[1]:
    st.header("3. Results")
    if "last_run" in st.session_state and st.session_state["last_run"]:
        last = st.session_state["last_run"]
        st.write(last["tle_source"])
        if last["count"] > 0:
            st.success(f"{last['count']} satellite(s) passed over the AOI.")
            for name, t in last["sats"]:
                st.write(f"🛰️ {name} at {t}")
                show_data_links(name)
        else:
            st.warning("No satellites passed over the AOI.")

        st.pyplot(last["fig"])

        # PDF Download
        pdf_buffer = create_pdf_report_text_and_image(
            last["params"]["sat_type"],
            last["params"]["year"],
            last["params"]["month"],
            last["params"]["day"],
            last["params"]["swath_km"],
            last["tle_source"],
            last["passing_sats"],
            last["fig"]
        )
        st.download_button(
            label="📄 Download Daily Report (.pdf)",
            data=pdf_buffer,
            file_name="satellite_daily_report.pdf",
            mime="application/pdf"
        )

        # Shapefile Export
        lines = []
        for item in last["plot_data"]:
            lines.append({'geometry': item['trace_line'], 'satellite': item['name']})
        tracks_gdf = gpd.GeoDataFrame(lines, crs="EPSG:4326")

        aoi_boundary = gpd.GeoDataFrame({'geometry': last["aoi"].geometry.boundary}, crs="EPSG:4326")

        shp_export_folder = "./temp_shp_export"
        if os.path.exists(shp_export_folder):
            shutil.rmtree(shp_export_folder)
        os.makedirs(shp_export_folder)

        tracks_path = os.path.join(shp_export_folder, "satellite_ground_tracks.shp")
        aoi_path = os.path.join(shp_export_folder, "aoi_boundary.shp")
        tracks_gdf.to_file(tracks_path)
        aoi_boundary.to_file(aoi_path)

        zip_path = os.path.join(shp_export_folder, "satellite_passes_and_aoi.zip")
        with zipfile.ZipFile(zip_path, 'w') as zipf:
            for base_path in [tracks_path, aoi_path]:
                for ext in [".shp", ".shx", ".dbf", ".prj", ".cpg"]:
                    file = base_path.replace(".shp", ext)
                    if os.path.exists(file):
                        zipf.write(file, arcname=os.path.basename(file))

        with open(zip_path, "rb") as f:
            st.download_button(
                label="📥 Download Satellite Passes & AOI Shapefile (.zip)",
                data=f,
                file_name="satellite_passes_and_aoi.zip",
                mime="application/zip"
            )

    else:
        st.info("No analysis run yet. Please go to the 'Search & Parameters' tab and run analysis.")

# --- Disclaimer ---
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
