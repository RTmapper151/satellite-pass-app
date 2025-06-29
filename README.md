# Earth Observation Satellite Pass Analyzer

## Goal  
The goal of this project is to identify Earth observation satellite passes over defined Areas of Interest (AOIs) using publicly available TLE data. It combines orbital modeling through Skyfield with spatial analysis tools to calculate and visualize when and where satellites are overhead. The final product is an interactive Streamlit app that supports multiple AOIs, outputs pass information, and maps each satellite’s ground track with context.

---

## Updates

### 2025-06-29  
Today is focused on organizing and documenting. I started by updating the README and digital journal to reflect how far the project has come. With the core analysis in good shape — reliably detecting satellite passes using bounding boxes and swath widths — the priority now shifts to refining the user interface. The Streamlit app is fully functional but needs more polish and better user flow. I’m outlining specific UI/UX tasks to work through next so I can move into a final cleanup phase.

### 2025-06-14  
Rebuilt the analysis loop to handle multiple AOIs from a folder instead of a single shapefile. This change significantly improved batch processing and made comparisons across regions more efficient. I added automatic CRS standardization (EPSG:4326) to avoid projection issues, which had caused minor problems earlier. Each AOI now receives a randomized identifier for output labeling, with the intention of manually renaming as needed. The output console now prints satellite names and timestamps clearly, and the ground track maps are drawing consistently across runs.

### 2025-06-11  
Refined the pass detection logic to move beyond simple point-in-polygon tests. I implemented bounding box checks using estimated swath widths to simulate the coverage area of each satellite. This approach produces much more realistic results, especially for wide-swath sensors, and reduces the likelihood of missing actual passes. I also filtered out nighttime passes to better align with optical satellite utility. This update marked a key improvement in how meaningful and accurate the pass detections are.

### 2025-06-05  
Successfully integrated TLE data from CelesTrak’s resource group into the analysis pipeline. Skyfield is now calculating full-day satellite ephemerides, and I’ve implemented basic visualizations of ground tracks overlaid on AOIs. With the timestamps visible and the orbits appearing correct, this laid the groundwork for a repeatable, day-by-day analysis model. At this stage, testing was limited to a single AOI and satellite, but the structure was ready to scale.

### 2025-05  
Initial development focused on testing feasibility. I set up the local Python environment with key libraries (`geopandas`, `matplotlib`, `shapely`, and `skyfield`) and built a prototype that mapped a satellite orbit over one AOI. This early proof-of-concept confirmed that spatial and orbital tools could be combined effectively, and it motivated the decision to develop a full-scale analysis and visualization platform with a user-facing interface.
