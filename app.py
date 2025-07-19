import os
import shutil
import zipfile
import tempfile
import io
import streamlit as st
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from fpdf import FPDF
from PIL import Image
from shapely.geometry import Point, box, LineString
from skyfield.api import load
from datetime import date as dt_date, timedelta

# Preview AOI map
def preview_aoi_map(aoi):
    fig = plt.figure(figsize=(6,6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, edgecolor='black', alpha=0.3)
    ax.add_feature(cfeature.OCEAN, alpha=0.1)
    aoi.boundary.plot(ax=ax, color='blue', linewidth=2, label='AOI')
    bounds = aoi.total_bounds
    ax.set_extent([bounds[0]-2, bounds[2]+2, bounds[1]-2, bounds[3]+2], crs=ccrs.PlateCarree())
    ax.set_title("AOI Preview", fontsize=12)
    plt.tight_layout()
    return fig

# Create PDF report with text and image
def create_pdf_report_text_and_image(sat_type, date_label, swath_km, tle_source, passing_sats, fig):
    pdf = FPDF()
    pdf.set_auto_page_break(True, 15)
    pdf.add_page()
    pdf.set_fill_color(30,144,255)
    pdf.rect(0,0,pdf.w,20,'F')
    pdf.set_text_color(255,255,255)
    pdf.set_font("Arial", "B", 18)
    pdf.cell(0,15,"Satellite Pass Daily Report", ln=True, align='C')
    pdf.ln(10)
    pdf.set_text_color(0,0,0)
    pdf.set_font("Arial", "I", 11)
    pdf.multi_cell(0,8,f"TLE Source: {tle_source}")
    pdf.multi_cell(0,8,f"Satellite Group: {sat_type}")
    pdf.multi_cell(0,8,f"Date(s): {date_label}")
    pdf.multi_cell(0,8,f"Swath Width: {swath_km} km")
    pdf.ln(5)
    pdf.set_draw_color(200,200,200)
    pdf.line(pdf.l_margin, pdf.get_y(), pdf.w - pdf.r_margin, pdf.get_y())
    pdf.ln(5)
    pdf.set_font("Arial", "B", 14)
    pdf.cell(0,10,"Satellites Passing Over the AOI:", ln=True)
    pdf.set_font("Arial", "", 12)
    if passing_sats:
        for name, t, d in passing_sats:
            pdf.cell(5)
            pdf.cell(0,8,f"- {name} on {d} at {t}", ln=True)
    else:
        pdf.multi_cell(0,10,"No satellites passed over the area.")
    pdf.ln(10)
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmpfile:
        temp_img_path = tmpfile.name
        fig.savefig(temp_img_path, format='PNG', dpi=300)
    page_width = pdf.w - 2*pdf.l_margin
    pdf.image(temp_img_path, x=pdf.l_margin, w=page_width)
    os.remove(temp_img_path)
    pdf.alias_nb_pages()
    pdf.set_y(-15)
    pdf.set_font("Arial", "I", 8)
    pdf.set_text_color(128,128,128)
    pdf.cell(0,10,f"Page {pdf.page_no()}/{{nb}}", align='C')
    pdf_output = pdf.output(dest='S').encode('latin1')
    pdf_bytes = io.BytesIO(pdf_output)
    pdf_bytes.seek(0)
    return pdf_bytes

# Show search link for satellite data
def show_data_links(sat_name):
    search_url = f"https://www.google.com/search?q={sat_name.replace(' ', '+')}+satellite+data+download"
    st.markdown(f"🔎 [Search for data from {sat_name}]({search_url})")

# Download or load cached TLE data
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

# Create AOI bounding box GeoDataFrame
def create_aoi(min_lon, min_lat, max_lon, max_lat):
    return gpd.GeoDataFrame({'geometry':[box(min_lon, min_lat, max_lon, max_lat)]}, crs="EPSG:4326")

# Generate time intervals for Skyfield
def generate_times(year, month, day, interval_minutes):
    ts = load.timescale()
    minutes = np.arange(0, 24*60, interval_minutes)
    return ts.utc(year, month, day, 0, minutes)

# Find satellites passing over AOI
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
            plot_data.append({'swath_gdf': swath_gdf, 'trace_line': trace_line, 'color': colors(i), 'name': sat.name})
    return passing_sats, plot_data

# Plot satellite passes and AOI
def plot_results(aoi, plot_data, swath_width_km):
    fig = plt.figure(figsize=(10,7))
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
    ax.set_extent([bounds[0]-2, bounds[2]+2, bounds[1]-2, bounds[3]+2], crs=ccrs.PlateCarree())
    ax.set_title("Satellite Passes", fontsize=14)
    ax.legend(fontsize=8, loc='lower left')
    plt.tight_layout()
    return fig

st.title("Satellite Pass Finder")
st.markdown("This tool finds satellites that pass over your AOI and tells you how it got the data.")
st.markdown("## Select a Tab Below to Begin")

tabs = st.tabs(["Main", "About"])

with tabs[0]:
    st.header("1. Define AOI")

    use_advanced = st.checkbox("Advanced: Upload AOI Shapefiles (for GIS users)")

    if use_advanced:
        uploaded_files = st.file_uploader("Upload zipped shapefiles (.zip)", type="zip", accept_multiple_files=True)
        aoi_list = []

        if uploaded_files:
            for uploaded_zip in uploaded_files:
                with tempfile.TemporaryDirectory() as tmpdir:
                    zip_path = os.path.join(tmpdir, "aoi.zip")
                    with open(zip_path, "wb") as f:
                        f.write(uploaded_zip.getvalue())
                    with zipfile.ZipFile(zip_path, "r") as zip_ref:
                        zip_ref.extractall(tmpdir)
                    shp_files = [f for f in os.listdir(tmpdir) if f.endswith(".shp")]
                    for shp in shp_files:
                        try:
                            aoi_gdf = gpd.read_file(os.path.join(tmpdir, shp)).to_crs("EPSG:4326")
                            aoi_list.append(aoi_gdf)
                        except:
                            pass
            if aoi_list:
                aoi = gpd.GeoDataFrame(pd.concat(aoi_list, ignore_index=True), crs="EPSG:4326")
                st.success(f"Loaded {len(aoi_list)} AOI shapefile(s).")
                st.pyplot(preview_aoi_map(aoi))
            else:
                st.warning("Please upload at least one valid shapefile.")
                st.stop()
        else:
            st.info("Upload zipped shapefiles to proceed.")
            st.stop()

    else:
        use_multi_box = st.checkbox("Define multiple bounding boxes")

        if use_multi_box:
            if "boxes" not in st.session_state:
                st.session_state.boxes = []

            with st.form("add_box_form"):
                col1, col2 = st.columns(2)
                with col1:
                    min_lon = st.number_input("Min Longitude", key="min_lon_new")
                    min_lat = st.number_input("Min Latitude", key="min_lat_new")
                with col2:
                    max_lon = st.number_input("Max Longitude", key="max_lon_new")
                    max_lat = st.number_input("Max Latitude", key="max_lat_new")
                submitted = st.form_submit_button("Add Bounding Box")
                if submitted:
                    st.session_state.boxes.append((min_lon, min_lat, max_lon, max_lat))

            if st.session_state.boxes:
                st.markdown("### Current Boxes:")
                aoi_list = []
                for i, (min_lon, min_lat, max_lon, max_lat) in enumerate(st.session_state.boxes):
                    st.markdown(f"Box {i+1}: ({min_lon}, {min_lat}) to ({max_lon}, {max_lat})")
                    aoi_list.append(create_aoi(min_lon, min_lat, max_lon, max_lat))
                aoi = gpd.GeoDataFrame(pd.concat(aoi_list, ignore_index=True), crs="EPSG:4326")
                st.pyplot(preview_aoi_map(aoi))
            else:
                st.info("Add at least one box to continue.")
                st.stop()

        else:
            col1, col2 = st.columns(2)
            with col1:
                min_lon = st.number_input("Min Longitude", value=127.5)
                min_lat = st.number_input("Min Latitude", value=25.5)
            with col2:
                max_lon = st.number_input("Max Longitude", value=129.0)
                max_lat = st.number_input("Max Latitude", value=27.0)
            aoi = create_aoi(min_lon, min_lat, max_lon, max_lat)
            st.pyplot(preview_aoi_map(aoi))

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
        "Scientific": "Science satellites.",
        "NOAA": "NOAA weather and ocean monitoring satellites.",
        "GOES": "Geostationary weather satellites over the Americas.",
        "Planet": "planet.com satellites.",
        "Last 30 Days": "Satellites launched in the last 30 days.",
        "Active": "All currently operational satellites tracked by CelesTrak."
    }

    sat_type = st.selectbox("Satellite Group", options=list(group_options.keys()))
    st.caption(f"Description: {group_descriptions.get(sat_type, 'No description available.')}")
    tle_group = group_options[sat_type]

    from datetime import timedelta, date as dt_date

    date_range = st.date_input(
        "Select 1–3 Days",
        value=[dt_date.today()],
        min_value=dt_date.today() - timedelta(days=1),
        max_value=dt_date.today() + timedelta(days=3)
    )
    if isinstance(date_range, dt_date):
        date_range = [date_range]
    if len(date_range) > 3:
        st.warning("Limit range to 3 days or fewer.")
        st.stop()

    st.caption(f"TLE predictions use data from: {dt_date.today()}.")
    swath_km = st.slider("Swath Width (km)", 10, 100, 30)
    interval = st.slider("Time Interval (minutes)", 1, 60, 10)

    if st.button("Run Analysis"):
        progress_bar = st.progress(0)
        status_text = st.empty()

        tle_folder = "./.cache_tle"
        status_text.text("Downloading TLE data...")
        tle_path, tle_source = download_tle(tle_group, tle_folder)
        progress_bar.progress(20)

        status_text.text("Loading satellite data...")
        satellites = load.tle_file(tle_path)
        progress_bar.progress(40)

        swath_width_m = swath_km * 1000
        all_passes = []
        all_plot_data = []

        for d in date_range:
            status_text.text(f"Analyzing passes for {d.isoformat()}...")
            times = generate_times(d.year, d.month, d.day, interval)
            day_passes, day_plot_data = find_passing_sats(satellites, times, aoi, swath_width_m)
            for name, t in day_passes:
                all_passes.append((name, t, d.isoformat()))
            all_plot_data.extend(day_plot_data)

        progress_bar.progress(80)
        status_text.text("Plotting results...")
        fig = plot_results(aoi, all_plot_data, swath_km)
        progress_bar.progress(100)
        status_text.empty()
        progress_bar.empty()

        st.subheader("3. Results")
        st.write(tle_source)

        if all_passes:
            st.success(f"{len(all_passes)} satellite pass(es) found.")
            for name, t, d in all_passes:
                st.write(f"{name} on {d} at {t}")
                show_data_links(name)
        else:
            st.warning("No satellites passed over the AOI.")

        st.pyplot(fig)

        with st.spinner("Preparing downloads..."):
            date_label = f"{date_range[0]} to {date_range[-1]}" if len(date_range) > 1 else date_range[0].isoformat()
            pdf_buffer = create_pdf_report_text_and_image(tle_group, date_label, swath_km, tle_source, all_passes, fig)

            st.download_button(
                "Download Daily Report (.pdf)",
                data=pdf_buffer,
                file_name="satellite_daily_report.pdf",
                mime="application/pdf"
            )

            lines = [{'geometry': item['trace_line'], 'satellite': item['name']} for item in all_plot_data]
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
                    "Download Satellite Passes & AOI Shapefile (.zip)",
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
        Provides publicly accessible Two-Line Element (TLE) data for tracking satellites.

        **2. Skyfield**  
        A Python library for astronomy and satellite position calculations.

        **3. GeoPandas**  
        For handling geospatial data, creating and manipulating geometries.

        **4. Shapely**  
        For planar geometric object manipulation.

        **5. Cartopy**  
        For creating geographic maps.

        **6. Matplotlib**  
        For plotting geographic data.

        **7. Streamlit**  
        For building the web app UI.

        **8. FPDF**  
        For generating PDF reports.

        **9. Pillow**  
        For image handling.

        **10. Pandas**  
        For data management.

        **11. Zipfile**  
        For compressing shapefile components.

        **12. tempfile**  
        For temporary file management.

        ### How the Analysis Works
        Focuses on Low Earth Orbit satellites whose passes over an area can be predicted using TLE data. GEO satellites remain fixed relative to Earth's surface and are excluded.

        ### Contact
        Created by: Steven Littel  
        Email: scl323@nau.edu
        """
    )
st.markdown(
    """
    ---
    📢 **Disclaimer**

    Pass prediction accuracy depends on CelesTrak data quality and update frequency. Only satellites listed in the selected CelesTrak group are analyzed.

    This tool does **not** query all satellites in orbit.
    """
)
