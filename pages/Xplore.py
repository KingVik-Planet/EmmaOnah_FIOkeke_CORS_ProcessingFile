#import all the Modules

import streamlit as st
import folium
import numpy as np
from scipy.interpolate import griddata
from pyproj import Proj, transform
from streamlit_folium import folium_static
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
import io
from geopy.distance import geodesic
import pandas as pd
import itertools
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import numpy as np


# Function to calculate distance between two coordinates using Haversine formula
def calculate_distance(coord1, coord2):
    return geodesic(coord1, coord2).kilometers


# Function to convert LLH to UTM
def llh_to_utm(lat, lon, height):
    utm_zone = int((lon + 180) / 6) + 1  # UTM zone calculation
    p = Proj(proj="utm", zone=utm_zone, ellps="WGS84", preserve_units=False)
    utm_x, utm_y = p(lon, lat)

    # Adjust UTM coordinates for zones 31 and 33 to make them positive
    if utm_zone == 31:
        utm_x -= 500000  # Add 500000 meters
    elif utm_zone == 33:
        utm_x += 500000  # Subtract 500000 meters

    return utm_x, utm_y, height


# Function to add UTM coordinates to station data
def add_utm_coordinates(station_data):
    for station, data in station_data.items():
        lat, lon, height = data["lat"], data["lon"], data["height"]
        utm_x, utm_y, utm_z = llh_to_utm(lat, lon, height)
        data["utm_x"] = utm_x
        data["utm_y"] = utm_y
        data["utm_z"] = utm_z




# Function to convert LLH to Cartesian coordinates
def llh_to_cartesian(lat, lon, height):
    # Placeholder conversion, replace with actual LLH to Cartesian conversion logic
    cart_x = lon * 1000  # Placeholder conversion
    cart_y = lat * 1000  # Placeholder conversion
    cart_z = height * 1000  # Placeholder conversion
    return cart_x, cart_y, cart_z

# Add Cartesian coordinates to station data
def add_cartesian_coordinates(station_data):
    for station, data in station_data.items():
        lat, lon, height = data["lat"], data["lon"], data["height"]
        cart_x, cart_y, cart_z = llh_to_cartesian(lat, lon, height)
        data["cart_x"] = cart_x
        data["cart_y"] = cart_y
        data["cart_z"] = cart_z



# Known station data
station_data = {
    "ABAK": {"lat": 6.315055601, "lon": 8.122842036, "height": 49.771},
    "ABIA": {"lat": 5.524274111, "lon": 7.520314611, "height": 159.203},
    "ASAB": {"lat": 6.19500475, "lon": 6.719156056, "height": 54.877},
    # "GEOS": {"lat": 6.330945781, "lon": 5.638304722, "height": 91.143},
    "OSUN": {"lat": 7.752762361, "lon": 4.52582788, "height": 326.628},
    "WARR": {"lat": 5.567309747, "lon": 5.8100107, "height": 16.382}
}

# Add UTM coordinates to station data
add_utm_coordinates(station_data)

# Extracting latitudes, longitudes, and heights
lats = [station_data[station]["lat"] for station in station_data]
lons = [station_data[station]["lon"] for station in station_data]
heights = [station_data[station]["height"] for station in station_data]



#University of Nigeria Nsukka Image
st.image("images/UNN.png", caption=" University of Nigeria", width=150)


# Create Streamlit app Title
from datetime import datetime

# Get the current date and time in the desired format
current_datetime = datetime.now().strftime('Date: %Y - %m - %d: Time: %H %M %S %f')

# Create Streamlit app Title with current date and time
st.markdown(f"""
    <div style='text-align: center; background-color: lightgreen; border-radius: 5px; padding: 1px;'>
        <h1 style='color: blue;'>Okeke_Onah GNSS CORS Processing Center Report</h1>
        <p style='color: blue; font-size: 1.5em;'>Generated at: {current_datetime}</p>
    </div>
""", unsafe_allow_html=True)



# User input for new station
st.markdown("""
    <div style='text-align: center; background-color: lightgreen; border-radius: 10px; padding: 10px;'>
        <h4 style='color: blue;'>Details of The Processed Station Data</h4>
    </div>
""", unsafe_allow_html=True)

st.markdown("<h4 style='color: blue;'>Name</h4>", unsafe_allow_html=True)
new_station_name = st.text_input("Station Name", key="name_input")
st.markdown("<h4 style='color: blue;'>Latitude</h4>", unsafe_allow_html=True)
new_station_lat = st.number_input("Station latitude", format="%.10f", key="latitude_input")
st.markdown("<h4 style='color: blue;'>Longitude</h4>", unsafe_allow_html=True)
new_station_lon = st.number_input("Station longitude", format="%.10f", key="longitude_input")
st.markdown("<h4 style='color: blue;'>Height</h4>", unsafe_allow_html=True)
new_station_height = st.number_input("Station height", format="%.10f", key="height_input")


# new_station_name = st.text_input("Name")
# new_station_lat = st.number_input("Latitude", format="%.10f")
# new_station_lon = st.number_input("Longitude", format="%.10f")
# new_station_height = st.number_input("Height", format="%.10f")

# New station data
new_station = {
    "name": new_station_name,
    "lat": new_station_lat,
    "lon": new_station_lon,
    "height": new_station_height
}


# Calculate UTM coordinates for the new station
new_station["utm_x"], new_station["utm_y"], new_station["utm_z"] = llh_to_utm(new_station["lat"], new_station["lon"], new_station["height"])



############
# Calculate Cartesian coordinates for the new station
new_station["cart_x"], new_station["cart_y"], new_station["cart_z"] = llh_to_cartesian(new_station["lat"], new_station["lon"], new_station["height"])

# Add Cartesian coordinates to station data
add_cartesian_coordinates(station_data)


# Interpolation methods
methods = ['linear']


# Perform interpolation and display results
for method in methods:
    interpolated_lat = griddata((lats, lons), lats, (new_station["lat"], new_station["lon"]), method=method)
    interpolated_lon = griddata((lats, lons), lons, (new_station["lat"], new_station["lon"]), method=method)
    interpolated_height = griddata((lats, lons), heights, (new_station["lat"], new_station["lon"]), method)

    # st.subheader(f"Interpolation Result for method: {method.capitalize()}")
    st.markdown(f"""
        <div style='text-align: center; background-color: lightgreen; border-radius: 5px; padding: 0.5px;'>
             <h5 style='color: blue;'>Note: </h5> 
        </div>
    """, unsafe_allow_html=True)
    st.write("""
        "All coordinates are based on the International GNSS Service (IGS) realization of the International Terrestrial 
    Reference Frame 2014 reference frame. The provided ITRF2014 coordinates correspond to a mean epoch of 
    the site observations, and all coordinates are referenced to the ground mark."
        """)
    st.write(f"Processed Latitude: {interpolated_lat:.10f},Adjusted Latitude: {interpolated_lat:.10f}, Difference: {new_station['lat'] - interpolated_lat:.10f}")
    st.write(f"Processed Longitude: {interpolated_lon:.10f},Adjusted Longitude: {interpolated_lon:.10f}, Difference: {new_station['lon'] - interpolated_lon:.10f}")
    st.write(f"Processed Height: {interpolated_height:.10f},Adjusted Height: {interpolated_height:.10f}, Difference: {new_station['height'] - interpolated_height:.10f}")

    # Extracting latitudes, longitudes, and heights
    lats = [station_data[station]["lat"] for station in station_data]
    lons = [station_data[station]["lon"] for station in station_data]
    heights = [station_data[station]["height"] for station in station_data]

    # Interpolate the new station data for latitude and longitude using nearest method
    interpolated_lat = griddata((lats, lons), lats, (new_station["lat"], new_station["lon"]), method="nearest")
    interpolated_lon = griddata((lats, lons), lons, (new_station["lat"], new_station["lon"]), method="nearest")

    # Plot known station data and the new interpolated point
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(lons, lats, c='blue', label='Known Stations')
    ax.scatter(interpolated_lon, interpolated_lat, c='red', label='Interpolated Point')



######################################################################################################################


    # Draw lines connecting the interpolated point to its nearest known station
    for station in station_data:
        ax.plot([station_data[station]["lon"], interpolated_lon], [station_data[station]["lat"], interpolated_lat],
                color='gray', linestyle='--', alpha=0.5)

    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # ax.set_title('Interpolation and Adjusted Using Nearest Method')
    # st.markdown("""
    #     <div style='text-align: center; background-color: lightgreen; border-radius: 5px; padding: 1px;'>
    #         <h4 style='color: blue;'>Connection Relation of the Rover and Known Stations</h4>
    #     </div>
    # """, unsafe_allow_html=True)
    ax.legend()
    ax.grid(True)

    # Display the plot using Streamlit
    # st.pyplot(fig)
#######################################################################################################################


    st.markdown("""
        <div style='text-align: center; background-color: lightgreen; border-radius: 5px; padding: 1px;'>
            <h4 style='color: blue;'>Connection Relation of the Rover and Known Stations</h4>
        </div>
    """, unsafe_allow_html=True)
    # Calculate the distances between the new station and each known station
    distances = [calculate_distance((new_station["lat"], new_station["lon"]), (lat, lon)) for lat, lon in
                 zip(lats, lons)]

    # Find the index of the station with the smallest distance
    min_distance_index = distances.index(min(distances))

    # Highlight the point with the smallest distance
    ax.scatter(lons[min_distance_index], lats[min_distance_index], c='red', label='Nearest Station')

    # Annotate the point with the smallest distance
    nearest_station_name = list(station_data.keys())[min_distance_index]
    ax.annotate(f'Nearest Station: {nearest_station_name} :: {min(distances):.2f} km',
                (lons[min_distance_index], lats[min_distance_index]),
                xytext=(lons[min_distance_index] + 0.05, lats[min_distance_index] + 0.05),
                arrowprops=dict(facecolor='black', arrowstyle='->'))

    # Display the plot using Streamlit
    st.pyplot(fig)





    # Display distances to known stations in a table
    dist_table_data = []
    for station, data in station_data.items():
        coord_new_station = (new_station["lat"], new_station["lon"])
        coord_known_station = (data["lat"], data["lon"])
        distance = calculate_distance(coord_new_station, coord_known_station)
        dist_table_data.append([station, f"{distance:.2f} km"])

    dist_table_data_with_headers = [["From New Station To:", "Distance (km)"]] + dist_table_data
    # st.subheader("Distances to Known Stations")
    st.markdown("""
        <div style='text-align: center; background-color: lightgreen; border-radius: 5px; padding: 1px;'>
            <h4 style='color: blue;'>Distances to Known Stations to the New Station</h4>
        </div>
    """, unsafe_allow_html=True)
    st.table(dist_table_data_with_headers)
    st.write("""
    "Cumulative distances, all expressed in kilometers, computed from the newly determined coordinates to their 
    respective known stations. These measurements represent crucial insights into spatial relationships, offering 
    comprehensive understanding and valuable context regarding the geographic distribution and 
    connectivity of the surveyed points. Such precise distance computation serve as fundamental metrics for analyzing 
    spatial dynamics and facilitating informed decision-making processes in various domains and check."
    
    """)


























    # Display adjustment (geographic) in a table
    adjustment_geo_table_data = []
    for station, data in station_data.items():
        lat_diff = interpolated_lat - data["lat"]
        lon_diff = interpolated_lon - data["lon"]
        height_diff = interpolated_height - data["height"]
        adjustment_geo_table_data.append([station, f"{new_station['lat'] - lat_diff:.10f}", f"{new_station['lon'] - lon_diff:.10f}", f"{new_station['height'] - height_diff:.10f}"])

    adjustment_geo_table_data_with_headers = [["Station", "Adjusted Latitude", "Adjusted Longitude", "Adjusted Height"]] + adjustment_geo_table_data
    # st.subheader("Adjustment (Geographic)")
    st.markdown("""
        <div style='text-align: center; background-color: lightgreen; border-radius: 5px; padding: 1px;'>
            <h4 style='color: blue;'>Adjustment (Geographic)</h4>
        </div>
    """, unsafe_allow_html=True)
    st.table(adjustment_geo_table_data_with_headers)
    st.write("""
     "In this section, geoid-ellipsoidal separations are computed utilizing a sophisticated spherical harmonic synthesis 
     of the global Earth Gravitational Model 2008 (EGM2008) geoid. This advanced technique ensures precise and accurate 
     determination of the separation between the geoid and the reference ellipsoid, providing essential insights into the
     Earth's gravitational field variations."
     """)

    # Display adjustment (UTM) in a table
    adjustment_utm_table_data = []
    for station, data in station_data.items():
        utm_x_diff = new_station["utm_x"] - station_data[station]["utm_x"]
        utm_y_diff = new_station["utm_y"] - station_data[station]["utm_y"]
        utm_z_diff = new_station["utm_z"] - station_data[station]["utm_z"]
        adjustment_utm_table_data.append([station, f"{new_station['utm_x'] + utm_x_diff:.10f}", f"{new_station['utm_y'] + utm_y_diff:.10f}", f"{new_station['utm_z'] + utm_z_diff:.10f}"])

    adjustment_utm_table_data_with_headers = [["Station", "Adjusted Easting", "Adjusted Northing", "Adjusted Height"]] + adjustment_utm_table_data
    # st.subheader("Adjustment (UTM)")
    st.markdown("""
        <div style='text-align: center; background-color: lightgreen; border-radius: 5px; padding: 1px;'>
            <h4 style='color: blue;'>Adjusted Universial Tranverse Mecator - UTM</h4>
        </div>
    """, unsafe_allow_html=True)

    st.table(adjustment_utm_table_data_with_headers)
########################################################################################
    st.markdown("""
        <div style='text-align: center; background-color: lightgreen; border-radius: 5px; padding: 1px;'>
            <h4 style='color: blue;'>Credit:</h4>
        </div>
    """, unsafe_allow_html=True)
    st.write("""
     Leaflet (https://leafletjs.com) \n
     © OpenStreetMap (https://www.openstreetmap.org/copyright) contributors \n
     Department of Geoinformatics and Surveying, University of Nigeria.
    """)
##########################################################################################
# Create a Folium map centered around the mean of latitudes and longitudes
m = folium.Map(location=[np.mean(list(map(lambda x: x["lat"], station_data.values()))),
                         np.mean(list(map(lambda x: x["lon"], station_data.values())))],
               zoom_start=6, control_scale=True)

# Add markers for each station
for station, data in station_data.items():
    folium.Marker(
        location=[data["lat"], data["lon"]],
        popup=f"{station} ({data['lat']}, {data['lon']}, {data['height']})",
        icon=folium.Icon(color="blue"),
    ).add_to(m)

# Add marker for the new station
folium.Marker(
    location=[new_station["lat"], new_station["lon"]],
    popup=f"{new_station['name']} ({new_station['lat']}, {new_station['lon']}, {new_station['height']})",
    icon=folium.Icon(color="green"),
).add_to(m)

# Display the map using Streamlit's `st` function
st.markdown("""
    <div style='text-align: center; background-color: lightgreen; border-radius: 5px; padding: 1px;'>
        <h4 style='color: blue;'>Map with Stations (Known Stations [Blue] & New Station [Green])</h4>
    </div>
""", unsafe_allow_html=True)
folium_static(m)




# Report content
report_content = """
This document is a report of the GPS data processing undertaken by the Okeke_Onah GNSS CORS Processing Center.
The Okeke_Onah GNSS CORS Processing Center compute precise coordinates in International Terrestrial
Reference Frame (ITRF) anywhere in Nigeria Basically with UTM Zone 32 North, Nigeria. 
The Service is designed to process only dual frequency GPS phase data.
An overview of the GPS processing strategy is included in this report.
Please direct any correspondence to our email: okekeonahcorsprocessingunn@gmail.com

Geoinformatics and Surveying Department
University of Nigeria, Enugu Campus
Enugu North Local Government Area of Enugu State, Nigeria
Website: https://unn.edu.ng

Home Page: https://okekeonahcorscoordinate.streamlit.app/
"""

# Create Streamlit app Title with current date and time
st.markdown(f"""
    <div style='text-align: center; background-color: lightgreen; border-radius: 5px; padding: 1px;'>
        <h1 style='color: blue;'>Okeke_Onah GNSS CORS Processing Center Report</h1>
        <p style='color: blue; font-size: 1em;'>{report_content}</p>
        <p style='color: blue; font-size: 1em;'>Generated at: {current_datetime}</p>
    </div>
""", unsafe_allow_html=True)