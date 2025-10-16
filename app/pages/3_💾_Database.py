"""
Database Management Page
"""
import streamlit as st
import subprocess
import os
import duckdb
from pathlib import Path

st.set_page_config(page_title="Database Management", page_icon="💾", layout="wide")

st.title("💾 Database Management")

# Configuration
DUCKDB_PATH = "/mnt/data4/master.duckdb"
REPSEQ_FASTA = "/mnt/data4/repseq.faa"
S3_BUCKET = "recombia.planter"
S3_DB_KEY = "master.duckdb"

# Database information
st.subheader("Database Information")

col1, col2, col3 = st.columns(3)

with col1:
    db_exists = os.path.exists(DUCKDB_PATH)
    if db_exists:
        db_size = os.path.getsize(DUCKDB_PATH) / (1024**3)  # GB
        st.metric("Database Status", "✅ Found", f"{db_size:.2f} GB")
    else:
        st.metric("Database Status", "❌ Not Found", "")

with col2:
    fasta_exists = os.path.exists(REPSEQ_FASTA)
    if fasta_exists:
        fasta_size = os.path.getsize(REPSEQ_FASTA) / (1024**2)  # MB
        st.metric("FASTA Status", "✅ Found", f"{fasta_size:.1f} MB")
    else:
        st.metric("FASTA Status", "❌ Not Found", "")

with col3:
    if db_exists:
        try:
            with duckdb.connect(DUCKDB_PATH, read_only=True) as conn:
                seq_count = conn.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
            st.metric("Total Sequences", f"{seq_count:,}", "")
        except:
            st.metric("Total Sequences", "Error", "")
    else:
        st.metric("Total Sequences", "-", "")

st.divider()

# Download Database from S3
st.subheader("📥 Download Database from S3")

st.info(f"**S3 Location**: s3://{S3_BUCKET}/{S3_DB_KEY}")

if st.button("Download Database from S3", type="primary", use_container_width=True):
    s3_path = f"s3://{S3_BUCKET}/{S3_DB_KEY}"

    # Create directory if needed
    os.makedirs(os.path.dirname(DUCKDB_PATH), exist_ok=True)

    with st.spinner(f"Downloading database from {s3_path}..."):
        try:
            # Run AWS CLI download
            cmd = ['aws', 's3', 'cp', s3_path, DUCKDB_PATH, '--no-progress']

            process = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600  # 10 minute timeout
            )

            if process.returncode == 0:
                st.success(f"✅ Database downloaded successfully to {DUCKDB_PATH}")
                st.rerun()  # Refresh to update metrics
            else:
                st.error(f"Download failed: {process.stderr}")

        except subprocess.TimeoutExpired:
            st.error("Download timed out after 10 minutes")
        except Exception as e:
            st.error(f"Error downloading database: {str(e)}")

st.divider()

# Create Reference FASTA
st.subheader("📄 Create Reference FASTA")

st.info("""
This will extract all representative sequences from the database and create a FASTA file for MMSeqs2 searching.
Includes both existing cluster representatives AND unclustered sequences from new samples.
""")

if not db_exists:
    st.warning("⚠️ Database not found. Please download it first.")
else:
    if st.button("Create Reference FASTA", type="primary", use_container_width=True):
        with st.spinner("Extracting representative sequences from database..."):
            try:
                # Create directory if needed
                os.makedirs(os.path.dirname(REPSEQ_FASTA), exist_ok=True)

                # Extract sequences using same logic as extract_representative_sequences()
                # Include BOTH existing representatives AND unclustered sequences
                with duckdb.connect(DUCKDB_PATH, read_only=True) as conn:
                    query = """
                    SELECT
                       '>' || seqhash_id || chr(10) || sequence
                    FROM
                        sequences
                    WHERE
                        repseq_id = seqhash_id           -- Existing representatives
                        OR repseq_id IS NULL;            -- Unclustered sequences (new samples)
                    """
                    result = conn.execute(query).fetchall()

                # Write to file
                with open(REPSEQ_FASTA, 'w') as f:
                    for row in result:
                        f.write(row[0] + '\n')

                st.success(f"✅ Created reference FASTA with {len(result):,} sequences")
                st.info(f"File saved to: {REPSEQ_FASTA}")
                st.rerun()  # Refresh to update metrics

            except Exception as e:
                st.error(f"Error creating FASTA: {str(e)}")

st.divider()

# Database Statistics
if db_exists:
    st.subheader("📊 Database Statistics")

    with st.spinner("Loading statistics..."):
        try:
            with duckdb.connect(DUCKDB_PATH, read_only=True) as conn:
                # Get table sizes
                tables = ['sequences', 'annotations', 'clusters', 'cluster_members',
                         'sra_metadata', 'go_terms', 'ec_numbers', 'kegg_info']

                stats = []
                for table in tables:
                    try:
                        count = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
                        stats.append({'Table': table, 'Row Count': f"{count:,}"})
                    except:
                        stats.append({'Table': table, 'Row Count': 'N/A'})

                # Display as dataframe
                import pandas as pd
                stats_df = pd.DataFrame(stats)
                st.dataframe(stats_df, use_container_width=True, hide_index=True)

                # Additional statistics
                col1, col2, col3 = st.columns(3)

                with col1:
                    rep_count = conn.execute("SELECT COUNT(*) FROM sequences WHERE is_representative = TRUE").fetchone()[0]
                    st.metric("Representative Sequences", f"{rep_count:,}")

                with col2:
                    cluster_count = conn.execute("SELECT COUNT(DISTINCT cluster_id) FROM cluster_members").fetchone()[0]
                    st.metric("Total Clusters", f"{cluster_count:,}")

                with col3:
                    organism_count = conn.execute("SELECT COUNT(DISTINCT organism) FROM sra_metadata").fetchone()[0]
                    st.metric("Unique Organisms", f"{organism_count:,}")

        except Exception as e:
            st.error(f"Error loading statistics: {str(e)}")

# Help section
with st.expander("❓ Help"):
    st.markdown("""
    ### Database Management Tasks

    **Download Database from S3**
    - Downloads the latest master.duckdb from S3
    - Requires AWS credentials to be configured
    - May take several minutes depending on database size

    **Create Reference FASTA**
    - Extracts all representative sequences from the database
    - Creates a FASTA file for MMSeqs2 searching
    - Only includes sequences where `is_representative = TRUE`
    - Required for the Search functionality to work

    ### File Locations
    - **Database**: `/mnt/data4/master.duckdb`
    - **Reference FASTA**: `/mnt/data4/repseq.faa`
    - **Temp Directory**: `/mnt/data4/tmp`
    """)
