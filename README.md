# Satellite Pass Prediction App

A lightweight and accessible app for predicting when Earth Observation satellites pass over user-defined areas of interest (AOIs). Built with Python and Streamlit, this tool allows users to visualize satellite passes, export shapefiles and PDFs, and quickly assess satellite visibility windows over specific regions.

---

## Features

- Predict satellite passes over a custom AOI  
- Uses live TLE data from [CelesTrak](https://celestrak.org)  
- Supports bounding box input for AOIs  
- Interactive mapping with Cartopy  
- Export options: PDF summary report + shapefile of satellite ground tracks  
- Lightweight and deployable locally or via the web  

---

## Folder Structure

- satellite-pass-app/
- app.py # Main Streamlit app
- requirements.txt # Python dependencies
- README.md

---

## Installation

1. **Clone the repo:**

`git clone https://github.com/your-username/satellite-pass-app.git`
`cd satellite-pass-app`

2. **Set up a virtual environment (optional but recommended):**

`python -m venv venv`
`source venv/bin/activate` or `venv\Scripts\activate`

3. **Install dependencies:**

`pip install -r requirements.txt`

---

## How To Use

1. **Run the app:**

`streamlit run app.py`

2. **Enter bounding box coordinates for your area of interest (AOI).**

3. **Select parameters:**

4. **Click “Run Analysis”**

5.**View:**
- Interactive map of satellite passes
- Table of satellite names and pass times

6. **Export:**
- PDF summary report
- Ground track shapefile (ZIP)

___

## Output Files

- PDF Report: Contains AOI info, pass summary (sat name, pass time), and a map image.
- Shapefile: ZIP archive with polyline features representing each satellite’s ground track over your AOI.

___

## Requirements

- Python 3.9+
- Internet access (for TLE data from CelesTrak)
- streamlit
- geopandas
- matplotlib
- shapely
- skyfield
- numpy
- cartopy
- pandas
- fpdf
- pillow
- folium
- streamlit_folium

---

## Troubleshooting

1. Map not showing?
  - Check that you entered valid latitude/longitude coordinates.

2. App crashes on launch?
  - Ensure you installed all packages with `pip install -r requirements.txt`.

3. No passes found?
  - Try widening your bounding box or increasing the analysis interval.








