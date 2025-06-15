import streamlit as st
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point, box, LineString
from skyfield.api import load
import numpy as np

st.title("Satellite Pass Finder")
st.markdown("This tool finds Earth Observation satellites that pass over your area of interest.")

# === User Inputs ===
st.header("1. Define Search Area")

col1, col2 = st.columns(2)
with col1:
    min_lon = st.number_input("Min Longitude", value=123.0)
    min_lat = st.number_input("Min Latitude", value=20.0)
with col2:
    max_lon = st.number_input("Max Longitude", value=125.0)
    max_lat = st.number_input("Max Latitude", value=22.0)

aoi = gpd.GeoDataFrame({'geometry': [box(min_lon, min_lat, max_lon, max_lat)]}, crs="EPSG:4326")

st.header("2. Select Satellite Type and Parameters")

sat_type = st.selectbox("Choose Satellite Type", options=[
    "Earth Observation", "Scientific"
])

sat_group_urls = {
    "Earth Observation": 'https://celestrak.org/NORAD/elements/gp.php?GROUP=resource&FORMAT=tle',
    "Scientific": 'https://celestrak.org/NORAD/elements/gp.php?GROUP=scientific&FORMAT=tle'
}

tle_url = sat_group_urls[sat_type]

date = st.date_input("Select Date", value=None)
swath_km = st.slider("Swath Width (km)", min_value=10, max_value=100, value=30)
interval = st.slider("Time Interval (minutes)", min_value=1, max_value=60, value=10)

if st.button("Run Analysis"):
    if not date:
        st.warning("Please select a date.")
    else:
        # === Load Satellite Data ===
        year, month, day = date.year, date.month, date.day
        tle_url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP=resource&FORMAT=tle'
        satellites = load.tle_file(tle_url)
        ts = load.timescale()

        minutes = np.arange(0, 24 * 60, interval)
        times = ts.utc(year, month, day, 0, minutes)

        # === Plot Setup ===
        fig, ax = plt.subplots(figsize=(10, 7))
        aoi.boundary.plot(ax=ax, color='black', linewidth=2, label='AOI Bounding Box')

        passing_sats = []
        colors = plt.cm.get_cmap('tab20', len(satellites))
        swath_width_m = swath_km * 1000
        buffer_radius = swath_width_m / 2

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
                swath_gdf.plot(ax=ax, color='red', alpha=0.4)
                trace_line = LineString([(pt.x, pt.y) for pt in points])
                ax.plot(*trace_line.xy, color=colors(i), linewidth=1.5, label=sat.name)

        bounds = aoi.total_bounds
        ax.set_xlim(bounds[0] - 2, bounds[2] + 2)
        ax.set_ylim(bounds[1] - 2, bounds[3] + 2)
        ax.set_title("Satellite Passes", fontsize=14)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.grid(True)
        ax.legend(fontsize=8, loc='lower left', frameon=True)

        # === Show Results ===
        st.subheader("Satellites That Passed Over AOI")
        if passing_sats:
            for name, t in passing_sats:
                st.write(f"🛰️ {name} at {t}")
        else:
            st.write("❌ No satellites passed over the area.")

        st.pyplot(fig)

        output_file = "satellite_passes.txt"

        with open(output_file, "w") as f:
            f.write(f"Earth Observation Satellites that passed over the search area on {year}-{month:02d}-{day:02d} (with {swath_width_m/1000:.0f} km swath):\n\n")
            if passing_sats:
                for name, pass_time in passing_sats:
                    f.write(f"{name} passes over AOI at {pass_time}\n")
            else:
                f.write("No satellites passed over the over the search area on {year}-{month:02d}-{day:02d} (with {swath_width_m/1000:.0f} km swath.\n")

        print(f"\nResults saved to {output_file}")

        image_file = "satellite_passes_map.png"
        plt.savefig(image_file, dpi=300)
        print(f"Map saved as {image_file}")

        # --- Download the text file ---
        with open("satellite_passes.txt", "r") as f:
            st.download_button(
                label="📄 Download Pass List",
                data=f,
                file_name="satellite_passes.txt",
                mime="text/plain"
            )

        # --- Download the image file ---
        with open("satellite_passes_map.png", "rb") as f:
            st.download_button(
                label="🖼️ Download Map Image",
                data=f,
                file_name="satellite_passes_map.png",
                mime="image/png"
            )