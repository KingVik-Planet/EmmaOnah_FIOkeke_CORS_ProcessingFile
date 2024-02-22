import streamlit as st
from contact_email import contact_email

# st.markdown("<h1 style='text-align: center; color: blue;'>Okeke_Onah GNSS CORS Processing Center</h1>",
#             unsafe_allow_html=True)
st.markdown("""
    <div style='text-align: center; background-color: lightgreen; border-radius: 10px; padding: 10px;'>
        <h1 style='color: blue;'>Okeke_Onah GNSS CORS Processing Center</h1>
    </div>
""", unsafe_allow_html=True)

st.markdown("<h2 style='text-align: center; color: blue;'>Upload the CORS Data Here:</h2>",
            unsafe_allow_html=True)
content = """Hello, Thank you for choosing us to handle your data. 
Please follow these simple steps for a seamless experience:

1. Enter Your Name: Provide your full name so we can address you properly.

2. Enter Your Email: Double-check your email for accuracy. Any errors might result in a delay or failure to receive feedback.

3. Short Message or Remarks: If you have any additional comments or remarks, feel free to include them here.

4. Upload Your Data: Select the file you wish to upload for processing.

5. Submit: Once you've completed the above steps, hit the submit button.

6. Confirmation Email: You will receive a confirmation email once your submission is successfully received. Please allow 3 to 15 minutes for this.

7. Feedback Delivery: Expect to receive feedback on your data within 4 to 5 hours of submission.

8. No Feedback?: In the rare event that you don't receive feedback after 5 hours, please don't hesitate to contact us.

9. 24-hour Response: If you encounter any issues or have questions, rest assured, we'll provide you with feedback within 24 hours.

Thank you once again for choosing us. Have a wonderful and productive day!
"""
st.info(content)

# with st.form(key="email_forms"):
#    # user_name = st.text_input("Which Name Should I Address You With?: (Compulsory) ") *
#     user_email = st.text_input("Enter your Email Address: (Compulsory) ")
#     #raw_message = st.text_area("Your Message: (Optional) ") *
#
#     # Add file uploader for attachment
#     attachment = st.file_uploader("Upload Attachment (Compulsory)")
#
#     subject = "New email"
#     if user_email:
#         subject = f"New email from {user_email}"
#
#     # Create the email message
#     message = f"From: {user_email}"#\n{raw_message}" *
#
#     button = st.form_submit_button("Submit")
#
#     if button:
#         # Send email with or without attachment
#         if attachment is not None:
#             attachment_content = attachment.getvalue()
#             attachment_name = attachment.name
#             contact_email(subject, message, attachment_content, attachment_name)
#             st.success("Your Message Was sent Successfully. We will get back to you soonest.")
#             # Reset input fields
#             user_name = ""
#             user_email = ""
#             #raw_message = "" *
#         else:
#             contact_email(subject, message)
#             st.success("Your Message Was sent Successfully. We will get back to you soonest.")
#             # Reset input fields
#             user_name = ""
#             user_email = ""
#             raw_message = ""
#
#         # Send a confirmation message to the user's email
#         if user_email and attachment is not None:
#             confirmation_msg = f"Thank you for submitting '{attachment_name}'. We will get back to you soonest."
#             contact_email("Submission Confirmation", confirmation_msg, attachment_content=None)
#
#
# # for * to track off






import datetime
import streamlit as st
from contact_email import contact_email

# Function to generate reference number
def generate_reference_number(file_name):
    current_time = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    return f"{file_name}_{current_time}"

with st.form(key="email_forms"):
    user_email = st.text_input("Enter your Email Address: (Compulsory)")

    attachment = st.file_uploader("Upload Attachment (Compulsory)")

    subject = "New email"
    if user_email:
        subject = f"New email from {user_email}"

    message = f"From: {user_email}"

    if not user_email:
        st.error("Please enter your email address.")
    if not attachment:
        st.error("Please upload an attachment.")

    submit_button_clicked = st.form_submit_button("Submit")

    if submit_button_clicked:
        if user_email and attachment:
            if attachment is not None:
                attachment_content = attachment.getvalue()
                attachment_name = attachment.name
                reference_number = generate_reference_number(attachment_name)
                confirmation_msg = f"Thank you for submitting '{attachment_name}'. Your reference number is '{attachment_name}_{datetime.datetime.now().strftime('%Y%m%d%H%M%S%f')}'. We will get back to you soonest."
                contact_email("Submission Confirmation", confirmation_msg, None, attachment_name, user_email)

                contact_email(subject, message, attachment_content, attachment_name, user_email)
                st.success("Your Message Was sent Successfully. We will get back to you soonest.")
            else:
                reference_number = generate_reference_number("NoAttachment")
                confirmation_msg = f"Thank you for submitting your message. Your reference number is '{reference_number}'. We will get back to you soonest."
                contact_email("Submission Confirmation", confirmation_msg, None, reference_number, user_email)

                contact_email(subject, message, None, None, user_email)
                st.success("Your Message Was sent Successfully. We will get back to you soonest.")
