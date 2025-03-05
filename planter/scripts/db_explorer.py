import duckdb
import pandas as pd
import streamlit as st

# connect to the database
db_path = "/mnt/data2/planter_outputs/planter.duckdb"
con = duckdb.connect(db_path)

st.title("Explore Your Sequences Database")

# fetch unique values for dropdowns
sample_ids = (
    con.execute("SELECT DISTINCT sample_id FROM sequences")
    .fetchdf()["sample_id"]
    .tolist()
)
go_terms = (
    con.execute("SELECT DISTINCT go_term FROM go_terms")
    .fetchdf()["go_term"]
    .dropna()
    .tolist()
)
ec_numbers = (
    con.execute("SELECT DISTINCT ec_number FROM ec_numbers")
    .fetchdf()["ec_number"]
    .dropna()
    .tolist()
)
cog_categories = (
    con.execute("SELECT DISTINCT cog_category FROM annotations")
    .fetchdf()["cog_category"]
    .dropna()
    .tolist()
)

# filter selection components
selected_sample_id = st.selectbox("Sample ID", options=["All"] + sample_ids)
selected_go_terms = st.multiselect("GO Terms", options=go_terms)
selected_ec_numbers = st.multiselect("EC Numbers", options=ec_numbers)
selected_cog_category = st.selectbox("COG Category", options=["All"] + cog_categories)

# build the SQL query dynamically based on selections
query = """
    SELECT s.seqhash_id, s.sequence, s.sample_id, s.assembly_date, 
           a.preferred_name, a.cog_category, g.go_term, e.ec_number 
    FROM sequences s
    LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
    LEFT JOIN go_terms g ON s.seqhash_id = g.seqhash_id
    LEFT JOIN ec_numbers e ON s.seqhash_id = e.seqhash_id
    WHERE 1=1
"""

# add filtering conditions
if selected_sample_id != "All":
    query += f" AND s.sample_id = '{selected_sample_id}'"

if selected_go_terms:
    go_terms_filter = "', '".join(selected_go_terms)
    query += f" AND g.go_term IN ('{go_terms_filter}')"

if selected_ec_numbers:
    ec_numbers_filter = "', '".join(selected_ec_numbers)
    query += f" AND e.ec_number IN ('{ec_numbers_filter}')"

if selected_cog_category != "All":
    query += f" AND a.cog_category = '{selected_cog_category}'"

# add pagination parameters
page_size = 500
total_count = con.execute(f"SELECT COUNT(*) FROM ({query})").fetchone()[0]
max_pages = max((total_count + page_size - 1) // page_size, 1)
page = st.slider("Page", 1, max_pages, 1)
offset = (page - 1) * page_size

# finalize query with LIMIT and OFFSET
query += f" LIMIT {page_size} OFFSET {offset}"

# execute the query and display results
try:
    df = con.execute(query).df()
    st.write(f"Showing page {page} of {max_pages} (total records: {total_count})")
    st.dataframe(df)
except Exception as e:
    st.error(f"Query failed: {e}")
