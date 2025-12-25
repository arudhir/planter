"""
Planter Pipeline - Streamlit Web Interface

Main entry point for the Streamlit application.
"""
import streamlit as st
from pathlib import Path
from config import Config

# Page configuration
st.set_page_config(
    page_title="Planter Pipeline",
    page_icon="🌱",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .main > div {
        padding-top: 2rem;
    }
    .stButton > button {
        width: 100%;
    }
</style>
""", unsafe_allow_html=True)

# Logo and title
logo_path = Path(__file__).parent / "static/images/recombia-logo.png"
if logo_path.exists():
    st.image(str(logo_path), width=200)

st.title("🌱 Planter Pipeline")

st.markdown("""
Welcome to the Planter Pipeline web interface!

Use the sidebar to navigate between different tools:

- **🔍 Search**: Search protein sequences using MMSeqs2
- **⚙️ Pipeline**: Run the Planter pipeline on samples
- **💾 Database**: Manage the reference database and FASTA files
- **📊 Query**: Execute SQL queries against the database

Select a page from the sidebar to get started.
""")

# Show system information
with st.expander("System Information"):
    st.info(f"""
    **Database Path**: `{Config.DUCKDB_PATH}`
    **Reference FASTA**: `{Config.REPSEQ_FASTA}`
    **Temp Directory**: `/mnt/data4/tmp`
    """)
