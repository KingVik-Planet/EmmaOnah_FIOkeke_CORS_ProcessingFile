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
    "GEOS": {"lat": 6.330945781, "lon": 5.638304722, "height": 91.143},
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



st.markdown("<h3 style='color: blue;'>Name</h3>", unsafe_allow_html=True)
new_station_name = st.text_input("Station Name", key="name_input")
st.markdown("<h3 style='color: blue;'>Latitude</h3>", unsafe_allow_html=True)
new_station_lat = st.number_input("Station latitude", format="%.10f", key="latitude_input")
st.markdown("<h3 style='color: blue;'>Longitude</h3>", unsafe_allow_html=True)
new_station_lon = st.number_input("Station longitude", format="%.10f", key="longitude_input")
st.markdown("<h3 style='color: blue;'>Height</h3>", unsafe_allow_html=True)
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


# #############################99999999999999999999999999999999966666666666666666666666666666666666666666666666666666666666666666666666666666666
# # Create PDF function
# def create_pdf(method, interpolated_lat, interpolated_lon, interpolated_height):
#     buffer = io.BytesIO()
#     c = canvas.Canvas(buffer, pagesize=letter)
#     c.drawString(100, 750, f"Interpolation and Adjustment Report ({method.capitalize()})")
#     y = 700
#
#     c.drawString(100, y, f"Interpolation Result for method: {method.capitalize()}")
#     y -= 20
#     c.drawString(100, y, f"Inputted Latitude: {new_station['lat']:.10f}, Interpolated Latitude: {interpolated_lat:.10f}, Difference: {new_station['lat'] - interpolated_lat:.10f}")
#     y -= 20
#     c.drawString(100, y, f"Inputted Longitude: {new_station['lon']:.10f}, Interpolated Longitude: {interpolated_lon:.10f}, Difference: {new_station['lon'] - interpolated_lon:.10f}")
#     y -= 20
#     c.drawString(100, y, f"Inputted Height: {new_station['height']:.10f}, Interpolated Height: {interpolated_height:.10f}, Difference: {new_station['height'] - interpolated_height:.10f}")
#
#     # Adjustment (Geographic)
#     y -= 40
#     c.drawString(100, y, "Adjustment (Geographic)")
#     y -= 20
#     for station, data in station_data.items():
#         lat_diff = interpolated_lat - data["lat"]
#         lon_diff = interpolated_lon - data["lon"]
#         height_diff = interpolated_height - data["height"]
#         c.drawString(100, y, f"{station} to {new_station['name']} - Lat: {new_station['lat'] - lat_diff:.10f}, Lon: {new_station['lon'] - lon_diff:.10f}, Height: {new_station['height'] - height_diff:.10f}")
#         y -= 20
#
#     # Adjustment (UTM)
#     y -= 40
#     c.drawString(100, y, "Adjustment (UTM)")
#     utm_x_diff = new_station["utm_x"] - station_data[station]["utm_x"]
#     utm_y_diff = new_station["utm_y"] - station_data[station]["utm_y"]
#     utm_z_diff = new_station["utm_z"] - station_data[station]["utm_z"]
#     for station, data in station_data.items():
#         c.drawString(100, y, f"{station} to {new_station['name']} - UTM_X: {new_station['utm_x'] + utm_x_diff:.10f}, UTM_Y: {new_station['utm_y'] + utm_y_diff:.10f}, UTM_Z: {new_station['utm_z'] + utm_z_diff:.10f}")
#         y -= 20
#
#     c.showPage()
#     c.save()
#     buffer.seek(0)
#     return buffer
# #55555555555555555555555555555555555555555555555555555555555555555555555555555555



# Perform interpolation and display results
for method in methods:
    interpolated_lat = griddata((lats, lons), lats, (new_station["lat"], new_station["lon"]), method=method)
    interpolated_lon = griddata((lats, lons), lons, (new_station["lat"], new_station["lon"]), method=method)
    interpolated_height = griddata((lats, lons), heights, (new_station["lat"], new_station["lon"]), method)

    # st.subheader(f"Interpolation Result for method: {method.capitalize()}")
    st.markdown(f"""
        <div style='text-align: center; background-color: lightgreen; border-radius: 5px; padding: 1px;'>
            <h4 style='color: blue;'>Adjusted Result for method: {method.capitalize()}</h4>
        </div>
    """, unsafe_allow_html=True)
    st.write("""
    
    "All coordinates are based on the International GNSS Service (IGS) realization of the International Terrestrial 
    Reference Frame 2014 (ITRF2014) reference frame. The provided ITRF2014 coordinates correspond to a mean epoch of 
    the site observations, and all coordinates are referenced to the ground mark."
    
    """)
    st.write(f"Processed Latitude: {interpolated_lat:.10f},Adjusted Latitude: {interpolated_lat:.10f}, Difference: {new_station['lat'] - interpolated_lat:.10f}")
    st.write(f"Processed Longitude: {interpolated_lon:.10f},Adjusted Longitude: {interpolated_lon:.10f}, Difference: {new_station['lon'] - interpolated_lon:.10f}")
    st.write(f"Processed Height: {interpolated_height:.10f},Adjusted Height: {interpolated_height:.10f}, Difference: {new_station['height'] - interpolated_height:.10f}")

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
    "Cumulative distances, all expressed in kilometers, have been meticulously computed from the newly determined 
    coordinates to their respective known stations. These measurements represent crucial insights into spatial 
    relationships, offering comprehensive understanding and valuable context regarding the geographic distribution and 
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
    # # Display adjustment (Cartesian) in a table
    # adjustment_cart_table_data = []
    # for station, data in station_data.items():
    #     cart_x_diff = new_station["cart_x"] - station_data[station]["cart_x"]
    #     cart_y_diff = new_station["cart_y"] - station_data[station]["cart_y"]
    #     cart_z_diff = new_station["cart_z"] - station_data[station]["cart_z"]
    #     adjustment_cart_table_data.append([station, f"{new_station['cart_x'] + cart_x_diff:.10f}", f"{new_station['cart_y'] + cart_y_diff:.10f}", f"{new_station['cart_z'] + cart_z_diff:.10f}"])
    #
    # adjustment_cart_table_data_with_headers = [["Station", "Adjusted Cart_X", "Adjusted Cart_Y", "Adjusted Cart_Z"]] + adjustment_cart_table_data
    # st.markdown("""
    #     <div style='text-align: center; background-color: lightgreen; border-radius: 5px; padding: 1px;'>
    #         <h4 style='color: blue;'>Adjusted Cartesian Coordinates</h4>
    #     </div>
    # """, unsafe_allow_html=True)
    # st.table(adjustment_cart_table_data_with_headers)
    #



# ###################################################
# # Function to calculate positional uncertainty
# def calculate_positional_uncertainty(observation_errors, satellite_orbit_errors, atmospheric_delays):
#     # Assuming observation_errors, satellite_orbit_errors, and atmospheric_delays are normally distributed with given standard deviations
#     # Compute combined standard deviation using root sum of squares (RSS)
#     combined_std_dev = np.sqrt(observation_errors ** 2 + satellite_orbit_errors ** 2 + atmospheric_delays ** 2)
#
#     # For 95% confidence level, calculate the 95% confidence interval (CI) using the z-score (1.96 for 95% CI)
#     ci_95 = 1.96 * combined_std_dev
#
#     return ci_95
#
#
# # Example values for standard deviations (in meters)
# observation_errors = 0.1  # Assuming 0.1 meters observation error
# satellite_orbit_errors = 0.2  # Assuming 0.2 meters satellite orbit error
# atmospheric_delays = 0.3  # Assuming 0.3 meters atmospheric delay error
#
# # Calculate positional uncertainty for each station
# uncertainties = {}
# for station, data in station_data.items():
#     positional_uncertainty = calculate_positional_uncertainty(observation_errors, satellite_orbit_errors,
#                                                               atmospheric_delays)
#     uncertainties[station] = positional_uncertainty
#
# # Display results in a table
# uncertainty_table_data = [[station, data["lon"], data["lat"], data["height"], uncertainties[station]] for station, data
#                           in station_data.items()]
# uncertainty_table = pd.DataFrame(uncertainty_table_data, columns=["Station", "Long(East)", "Lat(North)", "Height(m)",
#                                                                   "Positional Uncertainty (95% C.L.)"])
#
# # Display the uncertainty table
# st.markdown("""
#     <div style='text-align: center; background-color: lightgreen; border-radius: 5px; padding: 1px;'>
#         <h4 style='color: blue;'>Positional Uncertainty (95% C.L.)</h4>
#     </div>
# """, unsafe_allow_html=True)
# st.table(uncertainty_table)

# #####################################################
# import itertools
#
#
# # Function to calculate ambiguity resolution per baseline
# def calculate_ambiguity_resolution(station_data):
#     resolved_ambiguities = {}
#     baseline_lengths = {}
#     for pair in itertools.combinations(station_data.keys(), 2):
#         station1, station2 = pair
#         # Here you would calculate the resolved ambiguities for the pair
#         # You can use any method or algorithm you have for ambiguity resolution
#         # For this example, let's assume a random value between 0 and 10 for demonstration purposes
#         resolved_ambiguities[pair] = np.random.randint(0, 11)
#
#         # Calculate baseline length in kilometers using geodesic distance
#         coord1 = (station_data[station1]["lat"], station_data[station1]["lon"])
#         coord2 = (station_data[station2]["lat"], station_data[station2]["lon"])
#         baseline_lengths[pair] = calculate_distance(coord1, coord2)
#
#     return resolved_ambiguities, baseline_lengths
#
#
# # Calculate ambiguity resolution per baseline
# resolved_ambiguities, baseline_lengths = calculate_ambiguity_resolution(station_data)
#
# # Display results in a table
# ambiguity_resolution_table_data = []
# for pair in resolved_ambiguities.keys():
#     total_possible_ambiguities = baseline_lengths[pair] * 2  # Each baseline has two stations
#     percentage_resolved = (resolved_ambiguities[pair] / total_possible_ambiguities) * 100
#     ambiguity_resolution_table_data.append(
#         [f"{pair[0]} - {pair[1]}", resolved_ambiguities[pair], f"{percentage_resolved:.2f}%", baseline_lengths[pair]])
#
# ambiguity_resolution_table = pd.DataFrame(ambiguity_resolution_table_data,
#                                           columns=["Baseline", "Ambiguities Resolved", "Percentage Resolved",
#                                                    "Baseline Length (km)"])
#
# # Display the ambiguity resolution table
# st.markdown("""
#     <div style='text-align: center; background-color: lightgreen; border-radius: 5px; padding: 1px;'>
#         <h4 style='color: blue;'>Ambiguity Resolution per Baseline</h4>
#     </div>
# """, unsafe_allow_html=True)
# st.table(ambiguity_resolution_table)
#
# #########################################################
##ambiguity_resolution for II


# # Function to calculate ambiguity resolution per baseline
# def calculate_ambiguity_resolution(station_data):
#     resolved_ambiguities = {}
#     baseline_lengths = {}
#     for pair in itertools.combinations(station_data.keys(), 2):
#         station1, station2 = pair
#         # Here you would calculate the resolved ambiguities for the pair
#         # You can use any method or algorithm you have for ambiguity resolution
#         # For this example, let's assume a random value between 0 and 10 for demonstration purposes
#         resolved_ambiguities[pair] = np.random.randint(0, 11)
#
#         # Calculate baseline length in kilometers using geodesic distance
#         coord1 = (station_data[station1]["lat"], station_data[station1]["lon"])
#         coord2 = (station_data[station2]["lat"], station_data[station2]["lon"])
#         baseline_lengths[pair] = calculate_distance(coord1, coord2)
#
#     return resolved_ambiguities, baseline_lengths
#
#
# # Calculate ambiguity resolution per baseline
# resolved_ambiguities, baseline_lengths = calculate_ambiguity_resolution(station_data)
#
# # Display results in a table
# ambiguity_resolution_table_data = [[f"{pair[0]} - {pair[1]}", resolved_ambiguities[pair], baseline_lengths[pair]] for
#                                    pair in resolved_ambiguities.keys()]
# ambiguity_resolution_table = pd.DataFrame(ambiguity_resolution_table_data,
#                                           columns=["Baseline", "Ambiguities Resolved", "Baseline Length (km)"])
#
# # Display the ambiguity resolution table
# st.markdown("""
#     <div style='text-align: center; background-color: lightgreen; border-radius: 5px; padding: 1px;'>
#         <h4 style='color: blue;'>Ambiguity Resolution per Baseline</h4>
#     </div>
# """, unsafe_allow_html=True)
# st.table(ambiguity_resolution_table)

# #########################################################
#     # Generate PDF button
#     pdf_buffer = create_pdf(method, interpolated_lat, interpolated_lon, interpolated_height)
#     st.download_button(label=f"Download PDF ({method.capitalize()})", data=pdf_buffer, file_name=f"interpolation_report_{method}.pdf", mime="application/pdf")
#     st.write("\n")
# ##################################################

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