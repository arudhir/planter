"""
MMSeqs2 Protein Sequence Search Page
"""
import streamlit as st
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.search import process_search_request
from config import Config

st.set_page_config(page_title="Search", page_icon="🔍", layout="wide")

st.title("🔍 MMSeqs2 Protein Sequence Search")

# Configuration
FASTA_PATH = Config.REPSEQ_FASTA
DUCKDB_PATH = Config.DUCKDB_PATH
EXAMPLE_FASTA_PATH = Config.EXAMPLE_FASTA

# Sequence input
st.subheader("Enter Protein Sequence")

# Initialize session state for sequence if not exists
if 'sequence' not in st.session_state:
    st.session_state.sequence = ""

col1, col2 = st.columns([4, 1])

with col1:
    sequence = st.text_area(
        "Protein Sequence (FASTA format or raw sequence)",
        value=st.session_state.sequence,
        height=200,
        placeholder="Enter your protein sequence here...",
        help="You can paste a FASTA format sequence (with >header) or just the raw sequence",
        key="sequence_input"
    )
    # Update session state when text area changes
    st.session_state.sequence = sequence

with col2:
    st.write("") # Spacing
    st.write("") # Spacing
    if st.button("Load Example", use_container_width=True):
        try:
            with open(EXAMPLE_FASTA_PATH, 'r') as f:
                example_seq = f.read().strip()
                st.session_state.sequence = example_seq
                st.rerun()
        except FileNotFoundError:
            st.error("Example file not found")

# Advanced parameters
with st.expander("⚙️ Advanced Search Parameters"):
    col1, col2 = st.columns(2)

    with col1:
        sensitivity = st.slider(
            "Sensitivity",
            min_value=1.0,
            max_value=7.5,
            value=4.0,
            step=0.5,
            help="Higher values = more sensitive but slower"
        )

        e_value = st.number_input(
            "E-value threshold",
            min_value=0.0,
            max_value=1.0,
            value=0.001,
            format="%.4f",
            help="Maximum E-value for matches"
        )

        coverage = st.number_input(
            "Coverage threshold",
            min_value=0.0,
            max_value=1.0,
            value=0.0,
            step=0.1,
            help="Minimum fraction of query covered"
        )

        max_seqs = st.number_input(
            "Maximum sequences",
            min_value=1,
            max_value=10000,
            value=300,
            step=10,
            help="Maximum hits to return"
        )

    with col2:
        alignment_mode = st.selectbox(
            "Alignment mode",
            options=[0, 1, 2, 3, 4],
            index=3,
            help="3 = backtrace (default)"
        )

        mask = st.selectbox(
            "Mask low complexity",
            options=[0, 1],
            index=1,
            help="1 = mask low complexity regions"
        )

        min_seq_id = st.number_input(
            "Minimum sequence identity",
            min_value=0.0,
            max_value=1.0,
            value=0.0,
            step=0.05,
            help="Minimum sequence identity for matches"
        )

# Search button
if st.button("🔍 Search", type="primary", use_container_width=True):
    if not sequence or not sequence.strip():
        st.error("Please enter a protein sequence")
    else:
        # Clean sequence (remove FASTA header if present)
        clean_sequence = sequence.strip()
        if clean_sequence.startswith('>'):
            lines = clean_sequence.split('\n')
            clean_sequence = ''.join(lines[1:])

        # Remove whitespace
        clean_sequence = ''.join(clean_sequence.split())

        # Validate sequence (basic check)
        if len(clean_sequence) < 10:
            st.error("Sequence too short (minimum 10 amino acids)")
        else:
            # Build search parameters
            search_params = {
                'sensitivity': sensitivity,
                'e_value': e_value,
                'coverage': coverage,
                'max_seqs': max_seqs,
                'alignment_mode': alignment_mode,
                'mask': mask,
                'min_seq_id': min_seq_id
            }

            # Run search
            with st.spinner(f"Searching {len(clean_sequence)} amino acids against database..."):
                results = process_search_request(
                    clean_sequence,
                    FASTA_PATH,
                    DUCKDB_PATH,
                    search_params
                )

            # Store results in session state
            st.session_state.search_results = results
            st.session_state.search_sequence = clean_sequence

# Display results
if 'search_results' in st.session_state:
    results = st.session_state.search_results

    if results.get('error'):
        st.error(f"Error: {results['error']}")
    elif results['data'].empty:
        st.warning("No results found")
    else:
        st.success(f"✅ Found {results['result_count']} matches!")

        # Show search parameters used
        with st.expander("Search Parameters Used"):
            params = results['params']
            cols = st.columns(4)
            cols[0].metric("Sensitivity", params['sensitivity'])
            cols[1].metric("E-value", f"{params['e_value']:.4f}")
            cols[2].metric("Coverage", params['coverage'])
            cols[3].metric("Max Seqs", params['max_seqs'])

        # Display results table
        st.subheader("Search Results")

        # Column selection (matching Flask app column order)
        display_columns = [
            'query', 'organism', 'sample_id', 'preferred_name', 'target', 'description', 'tseq',
            'cluster_members', 'cog_category',
            'pident', 'alnlen', 'mismatch', 'gapopen',
            'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits',
            'cluster_size'
        ]

        available_columns = [col for col in display_columns if col in results['data'].columns]

        # Format numeric columns
        df_display = results['data'][available_columns].copy()
        if 'pident' in df_display.columns:
            df_display['pident'] = df_display['pident'].round(2)
        if 'evalue' in df_display.columns:
            df_display['evalue'] = df_display['evalue'].apply(lambda x: f"{x:.2e}")
        if 'bits' in df_display.columns:
            df_display['bits'] = df_display['bits'].round(1)

        # Display dataframe with sorting
        st.dataframe(
            df_display,
            use_container_width=True,
            height=600
        )

        # Download button
        csv = results['data'].to_csv(index=False)
        st.download_button(
            label="📥 Download Results as CSV",
            data=csv,
            file_name="mmseqs2_search_results.csv",
            mime="text/csv",
            use_container_width=True
        )

        # Show sequence alignment for selected result
        with st.expander("🔬 View Sequence Alignments"):
            if 'tseq' in results['data'].columns:
                result_idx = st.selectbox(
                    "Select result to view alignment",
                    options=range(len(results['data'])),
                    format_func=lambda i: f"{i+1}. {results['data'].iloc[i]['target']} - {results['data'].iloc[i].get('organism', 'Unknown')}"
                )

                if result_idx is not None:
                    selected = results['data'].iloc[result_idx]
                    st.text(f"Query:  {st.session_state.search_sequence}")
                    st.text(f"Target: {selected.get('tseq', 'N/A')}")
                    st.text(f"Identity: {selected.get('pident', 0):.1f}% | E-value: {selected.get('evalue', 0):.2e}")
            else:
                st.info("Sequence alignments not available in results")
