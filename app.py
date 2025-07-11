import streamlit as st
import geopandas as gpd
from shapely.geometry import Point, box, LineString
from skyfield.api import load
import numpy as np
import os
from datetime import date as dt_date
import folium
from streamlit_folium import st_folium
import tempfile
import shutil
import zipfile
from fpdf import FPDF
import io

# --- Helper functions ---

def create_aoi(min_lon, min_lat, max_lon, max_lat):
    return gpd.GeoDataFrame({'geometry': [box(min_lon, min_lat, max_lon, max_lat)]}, crs="EPSG:4326")

def generate_times(year, month, day, interval_minutes):
    ts = load.timescale()
    minutes = np.arange(0, 24*60, interval_minutes)
    return ts.utc(year, month, day, 0, minutes)

def download_tle(group, save_folder, max_days=1.0):
    os.makedirs(save_folder, exist_ok=True)
    filename = f"{group}.tle"
    filepath = os.path.join(save_folder, filename)
    url = f'https://celestrak.org/NORAD/elements/gp.php?GROUP={group}&FORMAT=TLE'

    # Simplified check - download always for demo, add caching logic if desired
    try:
        import urllib.request
        urllib.request.urlretrieve(url, filepath)
        source = f"Downloaded fresh TLE from: {url}"
    except Exception as e:
        source = f"Failed to download TLE, fallback to existing if any: {filepath}"
    return filepath, source

def find_passing_sats(satellites, times, aoi, swath_width_m):
    passing_sats = []
    plot_data = []
    buffer_radius = swath_width_m / 2

    for sat in satellites:
        subpoint = sat.at(times).subpoint()
        lats = subpoint.latitude.degrees
        lons = subpoint.longitude.degrees
        alt_km = subpoint.elevation.km  # altitude for tooltips
        speed_km_s = sat.at(times).velocity.km_per_s  # velocity vector magnitude approx
        
        points = [Point(lon, lat) for lon, lat in zip(lons, lats)]
        sat_gdf = gpd.GeoDataFrame({'geometry': points}, crs="EPSG:4326")
        sat_gdf_3857 = sat_gdf.to_crs("EPSG:3857")
        swath_buffers = sat_gdf_3857.buffer(buffer_radius)
        swath_gdf = gpd.GeoDataFrame(geometry=swath_buffers, crs="EPSG:3857").to_crs("EPSG:4326")

        intersect = gpd.sjoin(swath_gdf, aoi, how='inner', predicate='intersects')

        if not intersect.empty:
            idx = intersect.index[0]
            passing_time = times[idx].utc_iso()
            passing_sats.append((sat.name, passing_time, alt_km[idx], speed_km_s[idx]))
            trace_line = LineString([(pt.x, pt.y) for pt in points])
            plot_data.append({
                'swath_gdf': swath_gdf,
                'trace_line': trace_line,
                'name': sat.name,
                'altitude_km': alt_km,
                'speed_km_s': speed_km_s
            })
    return passing_sats, plot_data

def create_folium_map(aoi, plot_data, show_swath):
    bounds = aoi.total_bounds  # [minx, miny, maxx, maxy]
    center_lat = (bounds[1] + bounds[3]) / 2
    center_lon = (bounds[0] + bounds[2]) / 2

    fmap = folium.Map(location=[center_lat, center_lon], zoom_start=7)

    # Add AOI boundary
    folium.GeoJson(aoi.geometry, name="AOI", style_function=lambda x: {'color':'blue', 'weight':3, 'fill': False}).add_to(fmap)

    for item in plot_data:
        # Add ground track line with tooltip (satellite name)
        coords = [(y, x) for x, y in item['trace_line'].coords]  # folium expects lat, lon
        folium.PolyLine(
            coords,
            color="red",
            weight=2,
            opacity=0.7,
            tooltip=f"Satellite: {item['name']}"
        ).add_to(fmap)

        # Add swath buffer polygons if toggle is on
        if show_swath:
            for geom in item['swath_gdf'].geometry:
                folium.GeoJson(
                    geom,
                    style_function=lambda feat: {'fillColor': 'orange', 'color': 'orange', 'weight':1, 'fillOpacity':0.2},
                    tooltip=f"Swath area of {item['name']}"
                ).add_to(fmap)
    return fmap

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
        for name, t, alt, spd in passing_sats:
            pdf.multi_cell(0, 10, f"{name} at {t} (Alt: {alt:.1f} km, Speed: {spd:.2f} km/s)")
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
st.markdown("This tool finds satellites passing over your AOI and displays an interactive map with details.")

# Session state for run flag and map toggle
if "run_analysis" not in st.session_state:
    st.session_state.run_analysis = False
if "show_swath" not in st.session_state:
    st.session_state.show_swath = True

tabs = st.tabs(["Main", "About"])

with tabs[0]:
    st.header("1. Define AOI")
    col1, col2 = st.columns(2)
    with col1:
        min_lon = st.number_input("Min Longitude", value=127.5)
        min_lat = st.number_input("Min Latitude", value=25.5)
    with col2:
        max_lon = st.number_input("Max Longitude", value=129.0)
        max_lat = st.number_input("Max Latitude", value=27.0)

    aoi = create_aoi(min_lon, min_lat, max_lon, max_lat)
    st.markdown("#### AOI Preview")
    st.write(aoi)  # Simple display; could add a static plot if desired

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

    sat_type = st.selectbox("Satellite Group", options=list(group_options.keys()))
    st.caption(f"📘 **Description:** {group_descriptions.get(sat_type, 'No description available.')}")

    tle_group = group_options[sat_type]

    date = st.date_input("Date", value=dt_date.today())
    swath_km = st.slider("Swath Width (km)", min_value=10, max_value=100, value=30)
    interval = st.slider("Time Interval (minutes)", min_value=1, max_value=60, value=10)

    # Toggle to show/hide swath buffers on map
    st.session_state.show_swath = st.checkbox("Show Satellite Swath Buffers", value=st.session_state.show_swath)

    if st.button("Run Analysis"):
        st.session_state.run_analysis = True

    if st.session_state.run_analysis:
        progress_bar = st.progress(0)
        status_text = st.empty()

        year, month, day = date.year, date.month, date.day
        tle_folder = "./.cache_tle"

        status_text.text("Downloading TLE data...")
        tle_path, tle_source = download_tle(tle_group, tle_folder)
        progress_bar.progress(20)

        status_text.text("Loading satellite data...")
        satellites = load.tle_file(tle_path)
        progress_bar.progress(40)

        status_text.text("Generating time intervals...")
        times = generate_times(year, month, day, interval)
        progress_bar.progress(60)

        status_text.text("Finding satellites passing over AOI...")
        swath_width_m = swath_km * 1000
        passing_sats, plot_data = find_passing_sats(satellites, times, aoi, swath_width_m)
        progress_bar.progress(80)

        status_text.text("Creating interactive map...")
        folium_map = create_folium_map(aoi, plot_data, st.session_state.show_swath)
        progress_bar.progress(100)

        status_text.empty()
        progress_bar.empty()

        st.subheader("3. Results")
        st.write(tle_source)
        if passing_sats:
            st.success(f"{len(passing_sats)} satellite(s) passed over the AOI.")
            for name, t, alt, spd in passing_sats:
                st.write(f"🛰️ {name} at {t} — Altitude: {alt:.1f} km, Speed: {spd:.2f} km/s")
                show_data_links(name)
        else:
            st.warning("No satellites passed over the AOI.")

        # Display Folium map with interactive passes
        st_folium(folium_map, width=700, height=500)

        # PDF and Shapefile downloads
        with st.spinner("Preparing PDF and Shapefile downloads..."):
            import matplotlib.pyplot as plt
            # Create a simple matplotlib plot for PDF (fallback)
            fig, ax = plt.subplots(figsize=(6,4))
            ax.set_title("Satellite Passes")
            ax.set_xlim(min_lon, max_lon)
            ax.set_ylim(min_lat, max_lat)
            for item in plot_data:
                xs, ys = zip(*item['trace_line'].coords)
                ax.plot(xs, ys, label=item['name'])
            ax.legend(fontsize=8)
            pdf_buffer = create_pdf_report_text_and_image(tle_group, year, month, day, swath_km, tle_source, passing_sats, fig)
            plt.close(fig)

            st.download_button(
                label="📄 Download Daily Report (.pdf)",
                data=pdf_buffer,
                file_name="satellite_daily_report.pdf",
                mime="application/pdf"
            )

            # Shapefile export
            lines = [{'geometry': item['trace_line'], 'satellite': item['name']} for item in plot_data]
            tracks_gdf = gpd.GeoDataFrame(lines, crs="EPSG:4326")
            aoi_boundary = gpd.GeoDataFrame({'geometry': aoi.geometry.boundary}, crs="EPSG:4326")

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

with tabs[1]:
    st.header("About This App")
    st.markdown(
        """
        ### APIs and Python Packages Used

        **1. CelesTrak (Satellite TLE Data API)**  
        - Website: [https://celestrak.org/](https://celestrak.org/)  
        - Provides publicly accessible Two-Line Element (TLE) data for tracking satellites. The app downloads TLE files from CelesTrak groups to compute satellite orbits and positions.
        
        **2. Skyfield**  
        - Website: [https://rhodesmill.org/skyfield/](https://rhodesmill.org/skyfield/)  
        - A Python library for astronomy and satellite position calculations. Skyfield loads TLE data and calculates precise satellite locations over time using orbital mechanics.
        
        **3. GeoPandas**  
        - Website: [https://geopandas.org/](https://geopandas.org/)  
        - Extends pandas to support geographic data. Used to handle geospatial data, create and manipulate geometries such as AOI bounding boxes, satellite ground tracks, and spatial queries.
        
        **4. Shapely**  
        - Website: [https://shapely.readthedocs.io/](https://shapely.readthedocs.io/)  
        - A Python package for manipulation and analysis of planar geometric objects. Used here to create points, lines, and buffers for satellite ground tracks and AOI intersection tests.
        
        **5. Folium & streamlit-folium**  
        - Websites: [https://python-visualization.github.io/folium/](https://python-visualization.github.io/folium/), [https://github.com/randyzwitch/streamlit-folium](https://github.com/randyzwitch/streamlit-folium)  
        - Folium creates interactive Leaflet.js maps. streamlit-folium embeds these maps in Streamlit apps.
        
        **6. FPDF (Python FPDF)**  
        - Website: [https://pyfpdf.github.io/fpdf2/](https://pyfpdf.github.io/fpdf2/)  
        - Generates PDF documents including satellite reports with images.
        
        **...and others used for plotting, file handling, UI, etc.**

        ### How the Analysis Works
        - This tool focuses on Low Earth Orbit (LEO) satellites that pass over your area throughout the day. It computes their ground tracks and swath coverage based on TLE orbital data.
        - The interactive map lets you explore satellite passes and their coverage areas with tooltips showing altitude and speed.

        ### Contact
        Created by: Steven Littel  
        Email: scl323@nau.edu
        """
    )

st.markdown(
    """
    ---
    📢 **Disclaimer**

    The accuracy of pass predictions depends on the quality and update frequency of CelesTrak's datasets. Only satellites listed in the selected CelesTrak group will be included in the analysis.

    This tool does **not** query all satellites in orbit — only those published and maintained by CelesTrak.
    """
)
