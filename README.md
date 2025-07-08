# Earth Observation Satellite Pass Analyzer

## Goal  
The goal of this project is to identify Earth observation satellite passes over defined Areas of Interest (AOIs) using publicly available TLE data. It combines orbital modeling through Skyfield with spatial analysis tools to calculate and visualize when and where satellites are overhead. The final product is an interactive Streamlit app that supports multiple AOIs, outputs pass information, and maps each satellite’s ground track with context.

---
## Skills Demonstrated

### Remote Sensing & Orbital Data
- Parsed and utilized real satellite Two-Line Element (TLE) data from CelesTrak.  
- Calculated satellite ground tracks and simulated sensor swaths using subpoint positions and buffered footprints.  
- Filtered analysis to focus on low Earth orbit (LEO) satellites relevant for Earth observation.

### Spatial Analysis & Geometry Handling
- Created and manipulated AOIs and satellite swath zones as GeoDataFrames using GeoPandas.  
- Performed spatial joins, intersects, and buffering operations to detect satellite coverage of AOIs.  
- Modeled spatial-temporal phenomena by linking orbital data to geographic locations over time.

### Coordinate Reference System (CRS) Management
- Applied accurate CRS transformations between geographic (EPSG:4326) and projected metric systems (EPSG:3857).  
- Ensured precision in distance-based operations like buffering by working in appropriate projections.  
- Maintained consistent coordinate systems throughout processing, visualization, and export.

### Geospatial Data Export & Engineering
- Generated multi-file shapefiles including satellite ground tracks and AOI boundaries with proper spatial metadata.  
- Packaged shapefiles into zipped archives compatible with GIS software like QGIS and ArcGIS.  
- Managed file I/O and temporary storage efficiently for reliable export workflows.

### Cartographic Visualization
- Produced clear, publication quality maps using Matplotlib and Cartopy with basemap layers, AOI outlines, and satellite paths.  
- Used dynamic map extents and color-coded symbology to enhance readability and spatial context.  
- Integrated geographical features such as coastlines, borders, land, and ocean layers.

### Web GIS Interface Development
- Built an interactive web app with Streamlit for user-driven geospatial analysis.  
- Enabled AOI input via bounding boxes, satellite group selection, and parameter adjustments all within a browser.  
- Delivered instant visualizations, reports, and GIS data downloads with responsive UI elements.

### Scientific Reporting & Reproducibility
- Automated generation of PDF reports including maps, satellite pass times, swath widths, and data source citations.  
- Embedded metadata to ensure analysis transparency and reproducibility for research and operational use.  
- Documented processing steps and data provenance clearly in user-facing outputs.

### Data Integration & API Usage
- Accessed live satellite orbital data from public sources with caching to improve performance.  
- Supported multiple satellite groups and categories, facilitating thematic and comparative analyses.  
- Linked users to external resources for additional satellite data exploration.

### Python Programming for GIS
- Developed modular, readable Python code combining spatial and non-spatial libraries (Skyfield, GeoPandas, Shapely, Cartopy, Matplotlib, Streamlit).  
- Implemented robust file management, including temporary files, shapefile creation, and zipped downloads.  
- Automated end-to-end workflows from data acquisition through analysis to output generation.

### Automated Geospatial Workflows
- Designed a fully automated pipeline to translate orbital data into actionable geospatial insights.  
- Reduced manual GIS processing by scripting spatial analysis, visualization, and reporting.  
- Structured code for extensibility, allowing easy addition of new satellite groups or AOI types.

### User-Centered Design & Practical Application
- Developed a user-friendly interface allowing non-experts to analyze satellite passes without GIS or orbital mechanics knowledge.  
- Enabled users to specify AOIs, satellite groups, swath widths, and time intervals to tailor analysis outputs.  
- Delivered multiple output formats: interactive maps, PDF reports, and GIS-ready shapefiles.  
- Bridged complex geospatial and remote sensing concepts into an accessible tool that provides clear, actionable results.  
- Demonstrated the ability to build spatial tools that solve real-world problems for users who want insights, not just data.

### Additional GIS & Technical Skills
- Applied temporal data handling by generating time intervals and associating satellite positions over time.  
- Used spatial indexing and efficient querying to optimize intersection checks for performance.  
- Managed visualization aesthetics and legends to improve map interpretability.  
- Integrated error handling and data validation to ensure robust user experience.  
- Balanced computational performance with spatial accuracy through adjustable parameters like time step and swath width.

---

## Latest Branch (Usable): data-input

---

# Features in data-input

## User Interface
- Interactive web app built with **Streamlit**.
- Tabbed navigation with "Main" and "About" sections.
- User inputs for AOI bounding box (min/max latitude and longitude).
- Preview map of AOI using **Matplotlib** and **Cartopy**.
- Satellite group selection from multiple predefined categories with descriptions.
- Date picker for selecting analysis date.
- Sliders for swath width (km) and time interval (minutes).
- Progress bar and status text to show analysis progress.

## Data Handling & Analysis
- Download and cache Two-Line Element (TLE) files from **CelesTrak**.
- Parse TLE data with **Skyfield** to load satellite orbital info.
- Generate time intervals throughout the day for satellite position calculations.
- Compute satellite subpoints (lat/lon) for each time step.
- Create buffered polygons representing satellite swath widths.
- Use **GeoPandas** spatial join to find satellite passes intersecting the AOI.
- Produce list of passing satellites with timestamps.
- Generate geometric ground track lines for visualization and export.

## Visualization
- Plot satellite passes, swath buffers, and AOI boundary on a map with **Cartopy**.
- Unique colors and labels for each satellite pass.
- Map features include coastlines, country borders, land, and ocean.
- AOI preview map for visual verification.

## Reporting & Export
- Generate downloadable PDF report with:
  - Satellite group and date info
  - TLE data source details
  - List of passing satellites and times
  - Embedded satellite pass map image
- Export satellite ground tracks and AOI boundary as **shapefiles**.
- Package shapefiles into a ZIP archive for download.
- Download buttons for PDF report and shapefile ZIP.

## Code Quality & User Experience
- Modular functions for clarity and reuse.
- Real-time progress feedback during analysis.
- File existence checks before reading/writing outputs.
- Detailed About tab explaining APIs, libraries, satellite orbits, and contact info.

---
#Errors:

1. Download a .pdf reset analysis.
---

## Updates

### 2025-07-01
* Added tabbed layout with Main and About tabs
* Header stays fixed above tabs
* Added progress bar with step descriptions during analysis
* Added spinner while generating PDF and shapefile downloads
* Added AOI preview map after coordinate input
* Show last run summary in expandable section
* Cleaned up button placement and flow
* Added clear download buttons for PDF and shapefile zip
* Created About tab with API info, disclaimers, and contact details


### 2025-06-29  
Focused on organizing and documenting. I started by updating the README and digital journal to reflect how far the project has come. With the core analysis in good shape reliably detecting satellite passes using bounding boxes and swath widths —the priority now shifts to refining the user interface. The Streamlit app is fully functional but needs more polish and better user flow. I’m outlining specific UI/UX tasks to work through next so I can move into a final cleanup phase.

### 2025-06-14  
Rebuilt the analysis loop to handle multiple AOIs from a folder instead of a single shapefile. This change significantly improved batch processing and made comparisons across regions more efficient. I added automatic CRS standardization (EPSG:4326) to avoid projection issues, which had caused minor problems earlier. Each AOI now receives a randomized identifier for output labeling, with the intention of manually renaming as needed. The output console now prints satellite names and timestamps clearly, and the ground track maps are drawing consistently across runs.

### 2025-06-11  
Refined the pass detection logic to move beyond simple point-in-polygon tests. I implemented bounding box checks using estimated swath widths to simulate the coverage area of each satellite. This approach produces much more realistic results, especially for wide-swath sensors, and reduces the likelihood of missing actual passes. I also filtered out nighttime passes to better align with optical satellite utility. This update marked a key improvement in how meaningful and accurate the pass detections are.

### 2025-06-05  
Successfully integrated TLE data from CelesTrak’s resource group into the analysis pipeline. Skyfield is now calculating full-day satellite ephemerides, and I’ve implemented basic visualizations of ground tracks overlaid on AOIs. With the timestamps visible and the orbits appearing correct, this laid the groundwork for a repeatable, day-by-day analysis model. At this stage, testing was limited to a single AOI and satellite, but the structure was ready to scale.

### 2025-05  
Initial development focused on testing feasibility. I set up the local Python environment with key libraries (`geopandas`, `matplotlib`, `shapely`, and `skyfield`) and built a prototype that mapped a satellite orbit over one AOI. This early proof-of-concept confirmed that spatial and orbital tools could be combined effectively, and it motivated the decision to develop a full-scale analysis and visualization platform with a user-facing interface.

## Branches Overview

This project is organized into several branches, each representing a specific development phase or focus area. Below is a quick overview of each one:

### `main`
This is the original working branch where the core functionality was first developed. It includes early implementations like single AOI processing, basic TLE handling, and ground track plotting. The logic was simple but enough to prove the concept and test the workflow end-to-end.

### `new-selections`  
This branch added support for user selections and dynamic AOI inputs. It focused on expanding flexibility so users could upload or switch between AOIs more easily, paving the way for more interactive analysis.

### `downloader`  
Here we tackled and resolved issues related to downloading TLE data and managing local caching. This branch ensured that satellite data was fetched correctly and consistently, especially in environments with inconsistent internet access or frequent reruns.

### `2-modules`  
This branch marked the first major structural change — we split the project into separate modules to make the codebase easier to manage. This laid the foundation for maintainability and future expansion, especially as more features and options were being added.

### `data_pull_exp`  
This became the primary experimental and feature expansion branch. It includes support for multiple satellite groups, custom satellite selection, disclaimer handling, UI polishing, and several usability improvements. Most of the recent development work happens here before it gets cleaned up and merged elsewhere.

---

## Error Logging

- FileNotFoundError — 2025-05-18  
- FileNotFoundError — 2025-05-19  
- AttributeError — 2025-05-21  
- AttributeError — 2025-05-22  
- ImportError / ModuleNotFoundError — 2025-05-20  
- ImportError / ModuleNotFoundError — 2025-05-25  
- ValueError — 2025-05-21  
- ValueError — 2025-05-27  
- CRS Errors — 2025-05-22  
- CRS Errors — 2025-05-29  
- Streamlit UI Issues — 2025-05-23  
- Streamlit UI Issues — 2025-05-30  
- Memory Errors — 2025-05-24  
- Memory Errors — 2025-06-02  
- Shapefile Export Errors — 2025-05-25  
- Shapefile Export Errors — 2025-06-03  
- PDF Generation Errors — 2025-05-26  
- PDF Generation Errors — 2025-06-05  
- Network Errors — 2025-05-27  
- Network Errors — 2025-06-01  
- KeyError — 2025-05-28  
- KeyError — 2025-06-09  
- TimeoutError — 2025-05-29  
- TimeoutError — 2025-06-12  
- UnicodeEncodeError — 2025-05-30  
- UnicodeEncodeError — 2025-06-14  
- Spatial Join Errors — 2025-05-31  
- Spatial Join Errors — 2025-06-15  
- Datetime Handling Errors — 2025-06-01  
- Datetime Handling Errors — 2025-06-28  
- PermissionError — 2025-06-02  
- PermissionError — 2025-06-19  
- IndexError — 2025-06-03  
- IndexError — 2025-06-21  
- RuntimeWarning — 2025-06-04  
- RuntimeWarning — 2025-06-23  
- RecursionError — 2025-06-05  
- RecursionError — 2025-06-25  
- ConnectionResetError — 2025-06-06  
- ConnectionResetError — 2025-06-27  
- TypeError — 2025-06-07  
- TypeError — 2025-06-29  
- AssertionError — 2025-06-08  
- AssertionError — 2025-07-01  
- JSONDecodeError — 2025-06-09  
- JSONDecodeError — 2025-07-03  
- ModuleDeprecationWarning — 2025-06-10  
- ModuleDeprecationWarning — 2025-07-05  
- ZeroDivisionError — 2025-06-11  
- ZeroDivisionError — 2025-07-07  
- UserWarning — 2025-06-12  
- UserWarning — 2025-07-08  
- SyntaxError — 2025-06-13  
- FileNotFoundError — 2025-06-14  
- AttributeError — 2025-06-15  
- ImportError / ModuleNotFoundError — 2025-06-16  
- ValueError — 2025-06-17  
- CRS Errors — 2025-06-18  
- Streamlit UI Issues — 2025-06-19  
- Memory Errors — 2025-06-20  
- Shapefile Export Errors — 2025-06-21  
- PDF Generation Errors — 2025-06-22  
- Network Errors — 2025-06-23  
- KeyError — 2025-06-24  
- TimeoutError — 2025-06-25  
- UnicodeEncodeError — 2025-06-26  
- Spatial Join Errors — 2025-06-27  
- PermissionError — 2025-06-29  
- IndexError — 2025-06-30  
- RuntimeWarning — 2025-07-01  
- Streamlit UI Issues — 2025-07-02  
- RecursionError — 2025-07-02  
- TimeoutError — 2025-07-06  
- PermissionError — 2025-07-06  
- PermissionError — 2025-07-01  
- PermissionError — 2025-07-08  
- ImportError / ModuleNotFoundError — 2025-06-30  
- Network Errors — 2025-07-03  
- Performance Issues — 2025-06-29  
- Streamlit UI Issues — 2025-07-04  
- ZipFile Errors — 2025-06-22  
- IndexError — 2025-07-06  
- TypeError — 2025-06-23  
- SyntaxError — 2025-05-19  
- ZeroDivisionError — 2025-07-08  
- ModuleDeprecationWarning — 2025-06-03  
- UserWarning — 2025-06-20  
