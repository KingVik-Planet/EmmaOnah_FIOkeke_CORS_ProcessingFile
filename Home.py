import streamlit as st
import pandas as pd

col1, col2 = st.columns(2)

with col1:
    st.image("images/UNN.png", caption = "")
    st.markdown("<span style='color: blue; font-weight: bold;'>University of Nigeria Enugu Campus</span> <i>- Enugu State Nigeria</i>",
                unsafe_allow_html=True)

    about = """
    <div style="font-size: 20px; border: 1px solid #ccc; border-radius: 5px; padding: 10px; background-color: lightgreen;">
        <div style="display: flex; border-radius: 1px; overflow: hidden;">
            <a href="https://www.researchgate.net/institution/University-of-Nigeria2" style="text-decoration: none; color: blue; display: inline-block; padding: 5px 5px;">Research Gate</a>|
            <a href="(https://www.unn.edu.ng/academics/faculties/environmental-studies/geo-informatics-and-surveying/" style="text-decoration: none; color: blue; display: inline-block; padding: 5px 5px;">Dept of GSV</a>|
            <a href="https://unn.edu.ng" style="text-decoration: none; color: blue; display: inline-block; padding: 5px 5px;">UNN</a>
        </div>
    </div>
    """

    st.markdown(about, unsafe_allow_html=True)

    st.image("images/GPS.jpg", width = 345)

with col2:
#     # st.title("Okeke_Onah GNSS CORS Processing Center")
#     st.markdown("<h1 style='text-align: center; color: blue;'>Okeke_Onah GNSS CORS Processing Center</h1>",
#                 unsafe_allow_html=True)
    st.markdown("""
        <div style='text-align: center; background-color: lightgreen; border-radius: 10px; padding: 10px;'>
            <h1 style='color: blue;'>Okeke_Onah GNSS CORS Processing Center</h1>
        </div>
    """, unsafe_allow_html=True)

    content = """
Okeke_Onah GNSS CORS Processing Center, your premier online static data processing platform 
Powered by [Dept of Geoinformatics and Surveying, University of Nigeria](https://www.unn.edu.ng/academics/faculties/environmental-studies/geo-informatics-and-surveying/).
 - Processing Time: The processing time typically ranges from 2 to 5 hours.
 - Processing Methods: Our platform employs three different methods to ensure accurate coordinate adjustments and processing.
- Adjustments and Plotting: Upon processing, our system automatically plots the points on the map and provide Output Note in PDF through the submitted email. 

"""
    st.info(content)

with col1:
    import streamlit as st
    import folium
    import numpy as np
    from streamlit_folium import folium_static


    def main():
        # Station data
        station_data = {
            "ABAK": {"lat": 6.315055601, "lon": 8.122842036, "height": 49.771},
            "ABIA": {"lat": 5.524274111, "lon": 7.520314611, "height": 159.203},
            "ASAB": {"lat": 6.19500475, "lon": 6.719156056, "height": 54.877},
            "GEOS": {"lat": 6.330945781, "lon": 5.638304722, "height": 91.143},
            "OSUN": {"lat": 7.752762361, "lon": 4.52582788, "height": 326.628},
            "WARR": {"lat": 5.567309747, "lon": 5.8100107, "height": 16.382},
            "FUTM": {"lat": 9.5348783, "lon": 6.4509933, "height": 0.000},
            "IBAD": {"lat": 7.4141183, "lon": 3.9136054, "height": 0.000},
            "SACR": {"lat": 6.5952066, "lon": 3.3526114, "height": 0.000}
        }


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

        # Display the map using Streamlit's `st` function
        folium_static(m)


    if __name__ == "__main__":
        main()

#