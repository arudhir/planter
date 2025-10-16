"""
MMSeqs2 search utilities.

Contains all the backend logic for protein sequence searching.
"""
import subprocess
import tempfile
import os
import pandas as pd
import duckdb
from typing import Tuple, Dict, List


def run_mmseqs2_search(
    sequence: str,
    db_path: str,
    sensitivity: float = 4.0,
    e_value: float = 0.001,
    coverage: float = 0.0,
    max_seqs: int = 300,
    alignment_mode: int = 3,
    mask: int = 1,
    min_seq_id: float = 0.0
) -> Tuple[pd.DataFrame, str]:
    """
    Runs MMSeqs2 search and returns parsed results as a DataFrame.

    Args:
        sequence: The protein sequence to search
        db_path: Path to the reference database
        sensitivity: Search sensitivity (1.0-7.5, default 4.0)
        e_value: E-value threshold (default 0.001)
        coverage: Coverage threshold (0.0-1.0, default 0.0)
        max_seqs: Maximum number of hits per query (default 300)
        alignment_mode: Alignment detail level (0-4, default 3)
        mask: Low-complexity region masking (0 or 1, default 1)
        min_seq_id: Minimum sequence identity (0.0-1.0, default 0.0)

    Returns:
        Tuple of (results DataFrame, error message or None)
    """
    with tempfile.TemporaryDirectory(dir='/mnt/data4/tmp') as temp_dir:
        input_file = os.path.join(temp_dir, "input.fasta")
        output_file = os.path.join(temp_dir, "output.tsv")
        tmp_dir = os.path.join(temp_dir, "tmp")

        with open(input_file, 'w') as f:
            f.write(f">query\n{sequence}\n")

        mmseqs_command = [
            "mmseqs", "easy-search", input_file, db_path, output_file, tmp_dir,
            "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tseq",
            "-s", str(sensitivity),
            "-e", str(e_value),
            "-c", str(coverage),
            "--max-seqs", str(max_seqs),
            "--alignment-mode", str(alignment_mode),
            "--mask", str(mask),
            "--min-seq-id", str(min_seq_id),
            "-v", "0"  # Suppress verbose output
        ]

        try:
            process = subprocess.run(
                mmseqs_command,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
                timeout=60
            )

            if process.returncode != 0:
                return pd.DataFrame(), f"MMSeqs2 failed: {process.stderr}"

            if not os.path.exists(output_file):
                return pd.DataFrame(), "MMSeqs2 did not produce an output file"

            df = pd.read_csv(output_file, sep='\t', names=[
                'query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen',
                'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'tseq'
            ])

            return df, None

        except subprocess.TimeoutExpired:
            return pd.DataFrame(), "MMSeqs2 search timed out after 60 seconds"
        except Exception as e:
            return pd.DataFrame(), f"Error running MMSeqs2: {str(e)}"


def fetch_annotations_and_clusters(
    seqhash_ids: List[str],
    db_path: str
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Fetch annotations and cluster info from DuckDB.

    Args:
        seqhash_ids: List of sequence hash IDs to fetch
        db_path: Path to DuckDB database

    Returns:
        Tuple of (annotations DataFrame, clusters DataFrame)
    """
    if not seqhash_ids:
        return pd.DataFrame(), pd.DataFrame()

    try:
        with duckdb.connect(db_path) as db:
            # Fetch annotations
            annotations_query = f"""
                SELECT s.seqhash_id AS target, s.sample_id, s.length, m.organism,
                       a.description, a.cog_category, a.preferred_name
                FROM sequences s
                JOIN sra_metadata m ON s.sample_id = m.sample_id
                LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
                WHERE s.seqhash_id IN ({','.join([f"'{x}'" for x in seqhash_ids])})
            """
            annotations_df = db.execute(annotations_query).fetchdf()

            # Fetch clusters - for each hit sequence, get its cluster info including all members
            cluster_query = f"""
                SELECT
                    cm.seqhash_id AS target,
                    cm.cluster_id,
                    (SELECT GROUP_CONCAT(cm2.seqhash_id, ';')
                     FROM cluster_members cm2
                     WHERE cm2.cluster_id = cm.cluster_id) AS cluster_members,
                    (SELECT COUNT(*)
                     FROM cluster_members cm2
                     WHERE cm2.cluster_id = cm.cluster_id) AS cluster_size
                FROM cluster_members cm
                WHERE cm.seqhash_id IN ({','.join([f"'{x}'" for x in seqhash_ids])})
            """
            cluster_df = db.execute(cluster_query).fetchdf()

        return annotations_df, cluster_df

    except Exception as e:
        print(f"Database error: {e}")
        return pd.DataFrame(), pd.DataFrame()


def process_search_request(
    sequence: str,
    fasta_path: str,
    duckdb_path: str,
    search_params: Dict
) -> Dict:
    """
    Handles the entire search process and merges MMSeqs2 results with metadata.

    Args:
        sequence: Protein sequence to search
        fasta_path: Path to reference FASTA file
        duckdb_path: Path to DuckDB database
        search_params: Dictionary of search parameters

    Returns:
        Dictionary with headers and data
    """
    # Extract search parameters
    sensitivity = float(search_params.get('sensitivity', 4.0))
    e_value = float(search_params.get('e_value', 0.001))
    coverage = float(search_params.get('coverage', 0.0))
    max_seqs = int(search_params.get('max_seqs', 300))
    alignment_mode = int(search_params.get('alignment_mode', 3))
    mask = int(search_params.get('mask', 1))
    min_seq_id = float(search_params.get('min_seq_id', 0.0))

    # Run MMSeqs2 search
    mmseqs_df, error = run_mmseqs2_search(
        sequence, fasta_path,
        sensitivity=sensitivity,
        e_value=e_value,
        coverage=coverage,
        max_seqs=max_seqs,
        alignment_mode=alignment_mode,
        mask=mask,
        min_seq_id=min_seq_id
    )

    if error:
        return {'error': error, 'data': pd.DataFrame()}

    if mmseqs_df.empty:
        return {'error': 'No results found', 'data': pd.DataFrame()}

    # Get sequence IDs from results
    seqhash_ids = mmseqs_df['target'].unique().tolist()

    # Fetch annotations and clusters
    annotations_df, cluster_df = fetch_annotations_and_clusters(seqhash_ids, duckdb_path)

    # Merge annotations
    merged_df = mmseqs_df.merge(annotations_df, on='target', how='left')

    # Merge clusters
    merged_df = merged_df.merge(cluster_df, on='target', how='left')

    # Fill missing values
    merged_df.fillna({
        'organism': 'Unknown',
        'sample_id': 'Unknown',
        'description': 'No annotation',
        'cog_category': '',
        'preferred_name': '',
        'cluster_size': 1,
        'cluster_members': '',
        'cluster_id': ''
    }, inplace=True)

    return {
        'error': None,
        'data': merged_df,
        'params': search_params,
        'result_count': len(merged_df)
    }
