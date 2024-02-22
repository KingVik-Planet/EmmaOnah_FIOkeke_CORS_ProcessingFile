import streamlit as st

#Title
st.markdown("""
    <div style='text-align: center; background-color: lightgreen; border-radius: 10px; padding: 10px;'>
        <h1 style='color: blue;'>Meet Our Research Lab Project Team</h1>
    </div>
""", unsafe_allow_html=True)

subheader = ("""
Okeke_Onah GNSS Laboratory (OGL)is an online GPS data processing facility that's free to use. 
It uses localized CORS Network that covers some state in Nigeria. OGL is capable of working with data collected from any 
location within the network and provides differential corrections to any user within the network.
To utilize OGL, simply submit dual-frequency geodetic quality GPS RINEX data that's been observed in a 'static' mode to 
the GPS data processing system.
Once the data has been processed, an OGL report will be sent to you via email""")
st.info(subheader)

col1, col2 = st.columns(2)

with col1:

    st.image("images/FIO.jpg")
    st.markdown("<span style='color: blue; font-weight: bold;'>Professor Francis I. Okeke</span> <i> - Supervisor</i>",
                unsafe_allow_html=True)

    about = """
    <div style="font-size: 20px; border: 1px solid #ccc; border-radius: 5px; padding: 10px;">
        <div style="display: flex; border-radius: 1px; overflow: hidden;">
            <a href="https://www.researchgate.net/profile/Francis-Okeke-4" style="text-decoration: none; color: blue; display: inline-block; padding: 5px 5px;">Research Gate   </a>
            <a href="https://staffprofile.unn.edu.ng/profile/1854" style="text-decoration: none; color: blue; display: inline-block; padding: 5px 5px;">UNN-Staff Details</a>
        </div>
    </div>
    """

    st.markdown(about, unsafe_allow_html=True)


with col2:
    st.image("images/Onah.jpeg", width=290)
    st.markdown("<span style='color: blue; font-weight: bold;'>Surveyor Emmanuel U. Onah</span> <i>- PhD Student</i>",
                unsafe_allow_html=True)

    about = """
    <div style="font-size: 20px; border: 1px solid #ccc; border-radius: 5px; padding: 10px;">
        <div style="display: flex; border-radius: 1px; overflow: hidden;">
            <a href="https://www.researchgate.net/profile/Emmanuel-Udochukwu-Onah" style="text-decoration: none; color: blue; display: inline-block; padding: 5px 5px;">Research Gate   </a>
            <a href="https://linkedin.com/in/kingsley-chika-chukwu-235791154" style="text-decoration: none; color: blue; display: inline-block; padding: 5px 5px;">Staff Details</a>
        </div>
    </div>
    """
    st.markdown(about, unsafe_allow_html=True)