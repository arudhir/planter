"""
Database Query Page
"""
import streamlit as st
import duckdb
import pandas as pd
from pathlib import Path

st.set_page_config(page_title="Database Query", page_icon="📊", layout="wide")

st.title("📊 Database Query")

# Configuration
DUCKDB_PATH = "/mnt/data4/master.duckdb"
QUERY_DIR = Path("/home/ubuntu/planter/planter/database/queries/sql")

# Preset queries
PRESET_QUERIES = {
    "": "-- Select a preset query --",
    "database_summary": "Database Summary",
    "organism_summary": "Organism Summary",
    "sample_stats": "Sample Statistics",
    "cluster_stats": "Cluster Statistics",
    "go_term_summary": "GO Term Summary",
    "ec_number_summary": "EC Number Summary"
}

@st.cache_data
def load_preset_query(query_name: str) -> str:
    """Load a preset SQL query from file."""
    if not query_name:
        return ""

    query_file = QUERY_DIR / f"{query_name}.sql"
    if query_file.exists():
        with open(query_file, 'r') as f:
            return f.read()
    return f"-- Query file not found: {query_file}"

@st.cache_data(ttl=300)
def execute_query(query: str) -> tuple:
    """Execute a SQL query and return results."""
    # Safety checks
    unsafe_operations = ['CREATE', 'DROP', 'ALTER', 'INSERT', 'UPDATE', 'DELETE', 'TRUNCATE', 'VACUUM']
    query_upper = query.upper()

    for operation in unsafe_operations:
        if operation in query_upper:
            return None, f"Unsafe operation '{operation}' is not allowed"

    try:
        with duckdb.connect(DUCKDB_PATH, read_only=True) as conn:
            result = conn.execute(query).fetchdf()
        return result, None
    except Exception as e:
        return None, str(e)

# Query selection
col1, col2 = st.columns([3, 1])

with col1:
    selected_preset = st.selectbox(
        "Preset Queries",
        options=list(PRESET_QUERIES.keys()),
        format_func=lambda x: PRESET_QUERIES[x],
        help="Select a pre-made query or write your own below"
    )

with col2:
    st.write("") # Spacing
    st.write("") # Spacing
    if st.button("Clear Query", use_container_width=True):
        st.session_state.current_query = ""
        st.rerun()

# Load preset query if selected
if selected_preset and selected_preset != "":
    preset_query = load_preset_query(selected_preset)
    if 'current_query' not in st.session_state or st.session_state.get('last_preset') != selected_preset:
        st.session_state.current_query = preset_query
        st.session_state.last_preset = selected_preset

# SQL Editor
st.subheader("SQL Query Editor")

current_query = st.session_state.get('current_query', '')

query = st.text_area(
    "SQL Query",
    value=current_query,
    height=200,
    placeholder="SELECT * FROM sequences LIMIT 10;",
    help="Write your SQL query here. Read-only operations only."
)

# Update session state
if query != current_query:
    st.session_state.current_query = query

# Info about modifying preset queries
if selected_preset and selected_preset != "":
    st.info("💡 **Tip**: You can modify the query above. Replace example values with your own data.")

# Execute button
col1, col2, col3 = st.columns([1, 1, 2])

with col1:
    execute_btn = st.button("▶️ Run Query", type="primary", use_container_width=True)

with col2:
    if st.button("📋 Copy Query", use_container_width=True):
        st.code(query, language="sql")

# Execute query
if execute_btn:
    if not query or not query.strip():
        st.error("Please enter a SQL query")
    else:
        with st.spinner("Executing query..."):
            result_df, error = execute_query(query.strip())

        if error:
            st.error(f"Query failed: {error}")
        elif result_df is not None:
            # Store results in session state
            st.session_state.query_results = result_df
            st.session_state.query_text = query

# Display results
if 'query_results' in st.session_state:
    result_df = st.session_state.query_results

    st.success(f"✅ Query executed successfully! {len(result_df)} rows returned")

    # Show query that was executed
    with st.expander("📝 Executed Query"):
        st.code(st.session_state.query_text, language="sql")

    # Display results
    st.subheader("Query Results")

    # Pagination info
    total_rows = len(result_df)
    if total_rows > 1000:
        st.warning(f"⚠️ Large result set ({total_rows:,} rows). Showing first 1000 rows.")
        display_df = result_df.head(1000)
    else:
        display_df = result_df

    # Display dataframe
    st.dataframe(
        display_df,
        use_container_width=True,
        height=600
    )

    # Download button
    csv = result_df.to_csv(index=False)
    st.download_button(
        label=f"📥 Download Results as CSV ({total_rows:,} rows)",
        data=csv,
        file_name="query_results.csv",
        mime="text/csv",
        use_container_width=True
    )

    # Show statistics
    with st.expander("📈 Result Statistics"):
        st.write(f"**Rows**: {total_rows:,}")
        st.write(f"**Columns**: {len(result_df.columns)}")
        st.write(f"**Column Names**: {', '.join(result_df.columns.tolist())}")

        # Show data types
        st.write("**Data Types**:")
        st.dataframe(result_df.dtypes.to_frame(name='Type'), use_container_width=True)

# Help section
with st.expander("❓ Help & Examples"):
    st.markdown("""
    ### Available Tables
    - `sequences` - Protein sequences
    - `annotations` - Sequence annotations
    - `sra_metadata` - Sample metadata
    - `clusters` - Sequence clusters
    - `cluster_members` - Cluster membership
    - `go_terms` - GO term annotations
    - `ec_numbers` - EC number annotations
    - `kegg_info` - KEGG pathway information

    ### Example Queries

    **Get all samples:**
    ```sql
    SELECT DISTINCT sample_id, organism
    FROM sra_metadata
    LIMIT 100;
    ```

    **Search by organism:**
    ```sql
    SELECT s.seqhash_id, s.sample_id, m.organism, a.description
    FROM sequences s
    JOIN sra_metadata m ON s.sample_id = m.sample_id
    LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
    WHERE m.organism LIKE '%Rhodiola%'
    LIMIT 100;
    ```

    **Cluster statistics:**
    ```sql
    SELECT cluster_id, COUNT(*) as size
    FROM cluster_members
    GROUP BY cluster_id
    ORDER BY size DESC
    LIMIT 20;
    ```
    """)
