# Benchmark_29JUN
# This version represents a turning point in the project.
# Stable feature set with multiple AOI support, TLE parsing, pass detection with swath width, and Streamlit UI foundations.

import streamlit as st  # Build interactive UI
import geopandas as gpd  # Handle geospatial data
import matplotlib.pyplot as plt  # Plot maps
from shapely.geometry import Point, box, LineString  # Geometry tools
from skyfield.api import load  # Satellite orbital data
import numpy as np  # Numerical tools
import os  # File handling
from datetime import date as dt_date  # Handle dates

# --- Downloader with caching ---
def download_tle(group, save_folder, max_days=1.0):
    """Download or load cached TLE data for the specified group."""
    os.makedirs(save_folder, exist_ok=True)  # Ensure folder exists
    filename = f"{group}.tle"
    filepath = os.path.join(save_folder, filename)
    base_url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP={group}&FORMAT=TLE'
    url = base_url.format(group=group)

    if not load.exists(filepath) or load.days_old(filepath) >= max_days:
        load.download(url, filename=filepath)  # Fetch new data
        source = f"Downloaded fresh TLE from: {url}"
    else:
        source = f"Loaded cached TLE file: {filepath}"
    return filepath, source

# --- Create AOI from bounding box ---
def create_aoi(min_lon, min_lat, max_lon, max_lat):
    """Create a GeoDataFrame representing the AOI bounding box."""
    return gpd.GeoDataFrame({'geometry': [box(min_lon, min_lat, max_lon, max_lat)]}, crs="EPSG:4326")

# --- Generate time intervals ---
def generate_times(year, month, day, interval_minutes):
    """Generate a Skyfield time array for the given day at specified intervals."""
    ts = load.timescale()
    minutes = np.arange(0, 24*60, interval_minutes)
    return ts.utc(year, month, day, 0, minutes)

# --- Analyze satellite passes ---
def find_passing_sats(satellites, times, aoi, swath_width_m):
    """Find satellites passing over the AOI with given swath width."""
    passing_sats = []
    colors = plt.cm.get_cmap('tab20', len(satellites))  # Assign colors for plotting
    buffer_radius = swath_width_m / 2  # Buffer half swath width for intersection
    plot_data = []

    for i, sat in enumerate(satellites):
        subpoint = sat.at(times).subpoint()  # Compute subpoints for satellite over time
        lats = subpoint.latitude.degrees
        lons = subpoint.longitude.degrees
        points = [Point(lon, lat) for lon, lat in zip(lons, lats)]  # Convert to Point geometries

        sat_gdf = gpd.GeoDataFrame({'geometry': points}, crs="EPSG:4326")  # Create GeoDataFrame
        sat_gdf_3857 = sat_gdf.to_crs("EPSG:3857")  # Project to metric CRS for buffering
        swath_buffers = sat_gdf_3857.buffer(buffer_radius)  # Buffer points by swath radius
        swath_gdf = gpd.GeoDataFrame(geometry=swath_buffers, crs="EPSG:3857").to_crs("EPSG:4326")  # Back to lat/lon

        intersect = gpd.sjoin(swath_gdf, aoi, how='inner', predicate='intersects')  # Check intersection with AOI

        if not intersect.empty:
            passing_time = times[intersect.index[0]].utc_iso()  # Get first intersection time
            passing_sats.append((sat.name, passing_time))
            trace_line = LineString([(pt.x, pt.y) for pt in points])  # Create ground track line
            plot_data.append({
                'swath_gdf': swath_gdf,
                'trace_line': trace_line,
                'color': colors(i),
                'name': sat.name
            })

    return passing_sats, plot_data

# --- Plot map of satellite passes ---
def plot_results(aoi, plot_data, swath_width_km):
    """Plot AOI boundary and satellite ground tracks with swath coverage."""
    fig, ax = plt.subplots(figsize=(10, 7))
    aoi.boundary.plot(ax=ax, color='black', linewidth=2, label='AOI')

    for item in plot_data:
        item['swath_gdf'].plot(ax=ax, color='red', alpha=0.4)  # Plot swath buffers
        ax.plot(*item['trace_line'].xy, color=item['color'], linewidth=1.5, label=item['name'])  # Plot ground track

    bounds = aoi.total_bounds
    ax.set_xlim(bounds[0] - 2, bounds[2] + 2)
    ax.set_ylim(bounds[1] - 2, bounds[3] + 2)
    ax.set_title("Satellite Passes", fontsize=14)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(True)
    ax.legend(fontsize=8, loc='lower left', frameon=True)
    ax.text(bounds[0], bounds[3] + 1, f"Swath width: {swath_width_km} km", fontsize=12, fontweight='bold', color='navy')
    plt.tight_layout()
    return fig

# --- Streamlit UI setup ---
st.title("Satellite Pass Finder")
st.markdown("This tool finds satellites that pass over your AOI and tells you how it got the data.")

st.header("1. Define Search Area")
col1, col2 = st.columns(2)
with col1:
    min_lon = st.number_input("Min Longitude", value=123.0)
    min_lat = st.number_input("Min Latitude", value=20.0)
with col2:
    max_lon = st.number_input("Max Longitude", value=125.0)
    max_lat = st.number_input("Max Latitude", value=22.0)

aoi = create_aoi(min_lon, min_lat, max_lon, max_lat)

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
    "GPS": "Navigation satellites in the Global Positioning System constellation.",
    "Planet": "planet.com",
    "Iridium": "Communications satellites providing global voice and data coverage.",
    "Geodetic": "Satellites used for measuring Earth's shape, gravity, and geophysical phenomena.",
    "Last 30 Days": "Satellites Launched in the last 30 Days.",
    "Active": "All currently operational satellites tracked by CelesTrak."
}

sat_type = st.selectbox("Satellite Group", options=list(group_options.keys()))
st.caption(f"📘 **Description:** {group_descriptions.get(sat_type, 'No description available.')}")

tle_group = group_options[sat_type]

date = st.date_input("Date", value=dt_date.today())
swath_km = st.slider("Swath Width (km)", min_value=10, max_value=100, value=30)
interval = st.slider("Time Interval (minutes)", min_value=1, max_value=60, value=10)

if st.button("Run Analysis"):
    year, month, day = date.year, date.month, date.day
    tle_folder = "./.cache_tle"
    tle_path, tle_source = download_tle(tle_group, tle_folder)

    satellites = load.tle_file(tle_path)
    times = generate_times(year, month, day, interval)
    swath_width_m = swath_km * 1000

    passing_sats, plot_data = find_passing_sats(satellites, times, aoi, swath_width_m)
    fig = plot_results(aoi, plot_data, swath_km)

    st.subheader("3. Results")
    st.write(tle_source)
    if passing_sats:
        st.success(f"{len(passing_sats)} satellite(s) passed over the AOI.")
        for name, t in passing_sats:
            st.write(f"🛰️ {name} at {t}")
    else:
        st.warning("No satellites passed over the AOI.")

    st.pyplot(fig)

    # === File outputs ===
    output_text = "satellite_passes.txt"
    with open(output_text, "w") as f:
        f.write(f"TLE Source: {tle_source}\n\n")
        f.write(f"{sat_type} satellites over AOI on {year}-{month:02d}-{day:02d} ({swath_km} km swath):\n\n")
        if passing_sats:
            for name, t in passing_sats:
                f.write(f"{name} passes at {t}\n")
        else:
            f.write("No satellites passed over the area.\n")

    output_img = "satellite_passes_map.png"
    fig.savefig(output_img, dpi=300)

    with open(output_text, "r") as f:
        st.download_button("📄 Download Pass List", data=f, file_name=output_text, mime="text/plain")

    with open(output_img, "rb") as f:
        st.download_button("🖼️ Download Map Image", data=f, file_name=output_img, mime="image/png")

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
