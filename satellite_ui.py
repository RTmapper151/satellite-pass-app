import streamlit as st
from datetime import date

st.title("Satellite AOI Passes")

# Inputs from user
min_lon = st.number_input("Min Longitude", value=123.0)
min_lat = st.number_input("Min Latitude", value=20.0)
max_lon = st.number_input("Max Longitude", value=125.0)
max_lat = st.number_input("Max Latitude", value=22.0)

selected_date = st.date_input("Select date", value=date(2025, 7, 30))

swath_km = st.slider("Swath width (km)", min_value=1, max_value=100, value=30)

interval_min = st.number_input("Interval between points (minutes)", min_value=1, max_value=60, value=10)

if st.button("Run Analysis"):
    st.write("Running satellite pass analysis...")
    # Here you will call your satellite analysis function
    # For now, just show the inputs back
    st.write(f"Bounding Box: ({min_lon}, {min_lat}, {max_lon}, {max_lat})")
    st.write(f"Date: {selected_date}")
    st.write(f"Swath Width: {swath_km} km")
    st.write(f"Interval: {interval_min} minutes")
