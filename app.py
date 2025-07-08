import streamlit as st
import geopandas as gpd
from shapely.geometry import Point, box, LineString
from skyfield.api import load, EarthSatellite
import numpy as np
import os
from datetime import datetime, timedelta
import folium
from streamlit_folium import st_folium
import pandas as pd
import tempfile
import zipfile
import shutil

# --- Caching TLE download with expiration ---
@st.cache_data(ttl=24*3600)  # cache for 24 hours
def download_tle_cached(groups):
    base_url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP={group}&FORMAT=TLE'
    tle_data = {}
    for group in groups:
        url = base_url.format(group=group)
        try:
            response = load.download(url)
            tle_data[group] = response.decode('utf-8')
        except Exception as e:
            st.error(f"Error downloading TLE for group {group}: {e}")
            tle_data[group] = None
    return tle_data

# --- Parse TLE data string into Skyfield satellite objects ---
def parse_tle_string(tle_string):
    sats = []
    lines = tle_string.strip().splitlines()
    for i in range(0, len(lines), 3):
        if i+2 >= len(lines):
            break
        name = lines[i].strip()
        line1 = lines[i+1].strip()
        line2 = lines[i+2].strip()
        try:
            sat = EarthSatellite(line1, line2, name)
            sats.append(sat)
        except Exception as e:
            st.warning(f"Skipping satellite due to TLE parse error: {name} ({e})")
    return sats

# --- Create AOI GeoDataFrame from bounding box or polygon ---
def create_aoi_from_bounds(min_lon, min_lat, max_lon, max_lat):
    return gpd.GeoDataFrame({'geometry': [box(min_lon, min_lat, max_lon, max_lat)]}, crs="EPSG:4326")

def create_aoi_from_file(uploaded_file):
    try:
        gdf = gpd.read_file(uploaded_file)
        gdf = gdf.to_crs("EPSG:4326")
        return gdf
    except Exception as e:
        st.error(f"Error reading uploaded file: {e}")
        return None

# --- Generate times with custom start/end and interval ---
def generate_times_custom(start_dt, end_dt, interval_minutes):
    ts = load.timescale()
    total_minutes = int((end_dt - start_dt).total_seconds() / 60)
    minutes = np.arange(0, total_minutes + 1, interval_minutes)
    times = ts.utc(start_dt.year, start_dt.month, start_dt.day, start_dt.hour, start_dt.minute + minutes)
    return times

# --- Find satellites passing over AOI ---
def find_passing_sats(satellites, times, aoi, swath_width_m):
    passing_sats = []
    plot_lines = []
    buffer_radius = swath_width_m / 2

    for sat in satellites:
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
            first_pass_time = times[intersect.index[0]].utc_iso()
            passing_sats.append((sat.name, first_pass_time))
            trace_line = LineString([(pt.x, pt.y) for pt in points])
            plot_lines.append({'satellite': sat.name, 'geometry': trace_line})
    return passing_sats, gpd.GeoDataFrame(plot_lines, crs="EPSG:4326")

# --- Create Folium map showing AOI and satellite passes ---
def create_folium_map(aoi_gdf, passes_gdf):
    m = folium.Map(location=[aoi_gdf.geometry.centroid.y.values[0], aoi_gdf.geometry.centroid.x.values[0]], zoom_start=6)
    folium.GeoJson(aoi_gdf.geometry, name="AOI", style_function=lambda x: {'color':'black','weight':3,'fill':False}).add_to(m)
    for _, row in passes_gdf.iterrows():
        folium.GeoJson(row.geometry, name=row.satellite, style_function=lambda x: {'color':'red','weight':2}).add_to(m)
    folium.LayerControl().add_to(m)
    return m

# --- Export passes to CSV and GeoJSON ---
def export_passes(passes_gdf):
    csv_data = passes_gdf.to_csv(index=False)
    geojson_data = passes_gdf.to_json()
    return csv_data, geojson_data

# --- App UI ---
st.title("Satellite Pass Finder Plus")

st.header("1. Upload or Define AOI")
upload_option = st.radio("AOI input method", ["Bounding Box", "Upload Shapefile/GeoJSON"])

if upload_option == "Bounding Box":
    col1, col2 = st.columns(2)
    with col1:
        min_lon = st.number_input("Min Longitude", value=127.5)
        min_lat = st.number_input("Min Latitude", value=25.5)
    with col2:
        max_lon = st.number_input("Max Longitude", value=129.0)
        max_lat = st.number_input("Max Latitude", value=27.0)
    aoi = create_aoi_from_bounds(min_lon, min_lat, max_lon, max_lat)
else:
    uploaded_file = st.file_uploader("Upload AOI (shapefile zip or GeoJSON)", type=["zip", "geojson", "json"])
    if uploaded_file:
        if uploaded_file.name.endswith(".zip"):
            with tempfile.TemporaryDirectory() as tmpdir:
                zf = zipfile.ZipFile(uploaded_file)
                zf.extractall(tmpdir)
                files = [os.path.join(tmpdir, f) for f in os.listdir(tmpdir) if f.endswith(".shp")]
                if files:
                    aoi = gpd.read_file(files[0]).to_crs("EPSG:4326")
                else:
                    st.error("No shapefile found inside the ZIP.")
                    aoi = None
        else:
            aoi = create_aoi_from_file(uploaded_file)
    else:
        aoi = None

if aoi is not None:
    st.subheader("AOI Preview Map")
    fol_map = create_folium_map(aoi, gpd.GeoDataFrame(columns=['geometry'], crs="EPSG:4326"))
    st_data = st_folium(fol_map, width=700)

st.header("2. Select Satellite Groups")
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
selected_groups = st.multiselect("Satellite Groups (choose one or more)", options=list(group_options.keys()), default=["Earth Observation"])

if not selected_groups:
    st.warning("Please select at least one satellite group.")

st.header("3. Set Analysis Parameters")
date_start = st.date_input("Start Date", value=datetime.utcnow().date())
time_start = st.time_input("Start Time", value=datetime.utcnow().time())
date_end = st.date_input("End Date", value=date_start)
time_end = st.time_input("End Time", value=(datetime.utcnow() + timedelta(hours=1)).time())
swath_km = st.slider("Swath Width (km)", 10, 100, 30)
interval_min = st.slider("Time Interval (minutes)", 1, 60, 10)

run_analysis = st.button("Run Analysis")

if run_analysis:
    if aoi is None:
        st.error("AOI is not defined. Please upload or set a bounding box.")
    elif not selected_groups:
        st.error("Select at least one satellite group.")
    else:
        start_dt = datetime.combine(date_start, time_start)
        end_dt = datetime.combine(date_end, time_end)
        if end_dt <= start_dt:
            st.error("End datetime must be after start datetime.")
        else:
            with st.spinner("Downloading TLE data..."):
                tle_data_dict = download_tle_cached([group_options[g] for g in selected_groups if g in group_options])
            satellites = []
            for tle_string in tle_data_dict.values():
                if tle_string:
                    satellites.extend(parse_tle_string(tle_string))
            times = generate_times_custom(start_dt, end_dt, interval_min)
            swath_m = swath_km * 1000
            with st.spinner("Analyzing satellite passes..."):
                passing_sats, pass_lines_gdf = find_passing_sats(satellites, times, aoi, swath_m)

            st.subheader("Pass Results")
            if passing_sats:
                st.success(f"{len(passing_sats)} satellite(s) passed over the AOI in the selected time window.")
                for name, t in passing_sats:
                    st.write(f"🛰️ {name} at {t}")
            else:
                st.warning("No satellite passes found for your AOI and time window.")

            # Show passes on map
            fol_map = create_folium_map(aoi, pass_lines_gdf)
            st_folium(fol_map, width=700)

            # Export options
            csv_data, geojson_data = export_passes(pass_lines_gdf)
            st.download_button("Download Passes CSV", csv_data, "sat_passes.csv", "text/csv")
            st.download_button("Download Passes GeoJSON", geojson_data, "sat_passes.geojson", "application/geo+json")

