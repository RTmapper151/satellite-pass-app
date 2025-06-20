# tle_analysis.py
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point, box, LineString
from skyfield.api import load
import numpy as np
import os

def analyze_tle_passes(year, month, day, tle_folder, tle_group, aoi_bbox, swath_width_m=50000):
    """
    Analyze satellite passes over an AOI bounding box on a given date.
    
    Parameters:
    - year, month, day: int date for analysis
    - tle_folder: folder path where tle files are saved
    - tle_group: string group name, e.g. 'resource', 'weather'
    - aoi_bbox: tuple (minx, miny, maxx, maxy) in EPSG:4326
    - swath_width_m: integer swath width in meters (default 50 km)
    
    Returns:
    - List of tuples (satellite_name, passing_time_iso)
    - Displays a matplotlib plot of satellite tracks and swaths
    """
    
    # Create AOI GeoDataFrame
    aoi = gpd.GeoDataFrame({'geometry': [box(*aoi_bbox)]}, crs="EPSG:4326")

    # Load satellites from TLE file
    tle_path = os.path.join(tle_folder, f"{tle_group}.tle")
    satellites = load.tle_file(tle_path)
    ts = load.timescale()

    # Generate time range every 10 minutes for the day
    minutes = np.arange(0, 24 * 60, 10)
    times = ts.utc(year, month, day, 0, minutes)

    # Prepare plot
    fig, ax = plt.subplots(figsize=(12, 8))
    aoi.boundary.plot(ax=ax, color='black', linewidth=2, label="AOI")

    passing_sats = []
    colors = plt.cm.get_cmap('tab20', len(satellites))
    buffer_radius = swath_width_m / 2

    for i, sat in enumerate(satellites):
        subpoint = sat.at(times).subpoint()
        lats = subpoint.latitude.degrees
        lons = subpoint.longitude.degrees

        points = [Point(lon, lat) for lon, lat in zip(lons, lats)]
        sat_gdf = gpd.GeoDataFrame({'geometry': points}, crs="EPSG:4326")

        # Buffer the satellite points to create swath polygons
        sat_gdf_3857 = sat_gdf.to_crs("EPSG:3857")
        swath_buffers = sat_gdf_3857.buffer(buffer_radius)

        swath_gdf = gpd.GeoDataFrame(geometry=swath_buffers, crs="EPSG:3857").to_crs("EPSG:4326")

        # Spatial join/intersection with AOI
        intersect = gpd.sjoin(swath_gdf, aoi, how='inner', predicate='intersects')

        if not intersect.empty:
            passing_time = times[intersect.index[0]].utc_iso()
            passing_sats.append((sat.name, passing_time))

            # Plot swath (semi-transparent red)
            swath_gdf.plot(ax=ax, color='red', alpha=0.4)

            # Plot ground track line
            trace_line = LineString([(pt.x, pt.y) for pt in points])
            ax.plot(*trace_line.xy, color=colors(i), linewidth=1.5, label=sat.name)

    # Finalize plot limits and labels
    bounds = aoi.total_bounds
    ax.set_xlim(bounds[0] - 2, bounds[2] + 2)
    ax.set_ylim(bounds[1] - 2, bounds[3] + 2)

    ax.set_title(f"Satellites Passing AOI on {year}-{month:02d}-{day:02d}", fontsize=14)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(True)
    ax.legend(fontsize=8, loc='lower left', frameon=True)

    swath_text = f"Swath width: {swath_width_m / 1000:.0f} km"
    ax.text(bounds[0], bounds[3] + 1, swath_text, fontsize=12, fontweight='bold', color='navy')

    plt.tight_layout()
    plt.show()

    # Return the satellites that passed with times
    return passing_sats
