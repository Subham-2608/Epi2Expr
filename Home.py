import streamlit as st
import datetime
import os

st.set_page_config(page_title="EpiStackXpress", initial_sidebar_state="expanded", layout="wide")


# Visitor counter — stored in a text file
COUNTER_FILE = "static/data/visitor_count.txt"

def get_count():
    if not os.path.exists(COUNTER_FILE):
        return 0
    with open(COUNTER_FILE, "r") as f:
        try:
            return int(f.read().strip())
        except:
            return 0

def increment_count():
    count = get_count() + 1
    with open(COUNTER_FILE, "w") as f:
        f.write(str(count))
    return count

# Only count once per session
if "visited" not in st.session_state:
    st.session_state.visited = True
    visitor_count = increment_count()
else:
    visitor_count = get_count()


# IST time

IST_OFFSET = datetime.timezone(datetime.timedelta(hours=5, minutes=30))
ist_now = datetime.datetime.now(IST_OFFSET)
ist_str = ist_now.strftime("%d %B %Y  |  %I:%M:%S %p  IST")


# Header

col1, col2, col3 = st.columns([1.5, 20, 2])
with col1:
    st.image("static/images/icarlogo.png", width=150)
with col2:
    st.markdown("<h1 style='text-align:center;'>EpiStackXpress: Epigenetic Stacking Ensemble for Gene eXpression</h1>", unsafe_allow_html=True)
with col3:
    st.image("static/images/iasri-logo.png", width=150)

st.markdown("---")
st.text("")


# Background + Workflow

# If image needs more space
col1_1, col2_1 = st.columns([1, 1.2])


with col1_1:
    st.header("Background")
    st.markdown(
        "<p style='text-align: justify;'>"
        "Gene expression is a fundamental biological process that determines when, where, "
        "and how much a gene is active in a cell. Epigenetic modifications, such as DNA "
        "methylation, histone modifications, and chromatin remodeling, regulate gene "
        "expression without altering the underlying DNA sequence. These epigenetic marks "
        "play a critical role across all living organisms, influencing development, differentiation, "
        "and responses to environmental stimuli including biotic and abiotic stresses. Accurate "
        "prediction of gene expression from epigenetic features has emerged as a key challenge "
        "in computational biology with implications for understanding disease mechanisms, "
        "stress adaptation and genome regulation across plant and animal systems."
        "</p>",
        unsafe_allow_html=True
    )
    st.header("About EpiStackXpress")
    st.markdown(
        "<p style='text-align: justify;'>"
        "EpiStackXpress is a machine learning based web server for predicting gene expression "
        "levels (log₂FPKM) in Oryza sativa and related species using epigenetic features. It  "
        "employs stacking ensemble strategy combining Support Vector Regression (SVR) and Extreme "
        "Gradient Boosting (XGBoost) with a Ridge regression meta-learner. The server supports "
        "prediction under two experimental conditions — control and treatment- "
        "enabling researchers to study epigenetic regulation of gene expression "
        "in rice and related plants under biotic as wll as abiotic stress."
        "</p>",
        unsafe_allow_html=True
    )

with col2_1:
    st.markdown("<div style='padding-top: 95px;'>", unsafe_allow_html=True)
    st.image("static/images/Workflow.png", use_container_width=True)
    st.markdown("</div>", unsafe_allow_html=True)

st.text("")
st.markdown("---")
v1, v2, v3 = st.columns([1, 1, 1])
with v1:
    st.markdown(f"🕐 **Current Time:** {ist_str}")
with v3:
    st.markdown(f"<p style='text-align:right;'>👥 <b>Total Visitors:</b> {visitor_count:,}</p>", unsafe_allow_html=True)
st.markdown("---")
st.markdown("<div style='background-color:#32CD32; text-align:center'><p style='color:white'>Copyright © 2026 ICAR-Indian Agricultural Statistics Research Institute, New Delhi-110012. All rights reserved.</p></div>", unsafe_allow_html=True)