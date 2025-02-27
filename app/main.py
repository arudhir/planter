from flask import Flask, render_template, request, jsonify, current_app
import subprocess
import tempfile
import os
import logging
from pathlib import Path
import pandas as pd
import duckdb
from app.config import config
from planter.database.query_manager import DatabaseManager

def create_app(config_name='default'):
    """Creates and configures the Flask app."""
    app = Flask(__name__)
    app.config.from_object(config[config_name])
    logging.basicConfig(level=logging.DEBUG)

    @app.route('/', methods=['GET', 'POST'])
    def index():
        """Handles the main search page."""
        if request.method == 'POST':
            sequence = request.form['sequence']
            # import ipdb; ipdb.set_trace()
            results = process_search_request(sequence, app.config['REPSEQ_FASTA'], app.config['DUCKDB_PATH'])
            return render_template('results.html', results=results)
        return render_template('index.html')

    @app.route('/load_example', methods=['GET'])
    def load_example():
        """Loads an example sequence from a file."""
        try:
            with open(app.config['EXAMPLE_FASTA'], 'r') as f:
                example_sequence = f.read().strip()
            return jsonify({'sequence': example_sequence})
        except Exception as e:
            app.logger.error(f"Error loading example: {str(e)}")
            return jsonify({'error': 'Failed to load example sequence'}), 500

    def get_sequence_details(seqhash_ids):
        """Fetch annotations and expression for the given seqhash_ids from the database.
        
        Retrieves comprehensive information about sequences:
        - Annotations (GO terms, EC numbers, description, etc.)
        - Expression data (TPM values)
        - Gene to protein mappings
        
        Args:
            seqhash_ids: List of protein seqhash IDs to query
            
        Returns:
            Dictionary of sequence information indexed by seqhash_id
        """
        try:
            with DatabaseManager(app.config['DUCKDB_PATH']) as db:
                # First get the basic annotations using the canonical query
                annotations_df = db.queries.sequences.get_by_id(seqhash_ids[0])
                
                # Collect all sequence results
                results = {}
                
                for seqhash_id in seqhash_ids:
                    try:
                        # Get detailed information for this sequence
                        annotation = db.queries.sequences.get_by_id(seqhash_id)
                        
                        if not annotation.empty:
                            # Create a dictionary for this sequence and add it to results
                            # Expression data temporarily disabled
                            # try:
                            #     # Get expression data for the corresponding gene
                            #     expression_data = db.queries.sequences.get_annotation_with_expression(
                            #         sample_id=annotation.iloc[0].get('sample_id', None)
                            #     )
                            #     
                            #     # Filter expression data for this sequence
                            #     seq_expression = expression_data[expression_data['protein_id'] == seqhash_id]
                            #     
                            #     # Create a dictionary for this sequence
                            #     seq_dict = annotation.iloc[0].to_dict()
                            #     
                            #     # Add expression data if available
                            #     if not seq_expression.empty:
                            #         seq_dict.update({
                            #             'gene_id': seq_expression.iloc[0].get('gene_id', ''),
                            #             'tpm': seq_expression.iloc[0].get('tpm', 0),
                            #             'num_reads': seq_expression.iloc[0].get('num_reads', 0)
                            #         })
                            # except Exception as e:
                            #     # Log the error but continue with basic annotation
                            #     app.logger.error(f"Error processing expression data for {seqhash_id}: {e}")
                            
                            # Just use basic annotation data for now
                            seq_dict = annotation.iloc[0].to_dict()
                            results[seqhash_id] = seq_dict
                        else:
                            # Handle case where sequence is not found in database
                            app.logger.warning(f"Sequence {seqhash_id} not found in database")
                            results[seqhash_id] = {
                                'seqhash_id': seqhash_id,
                                'description': 'Sequence not found in database'
                            }
                    except Exception as e:
                        # Log any other errors during sequence processing
                        app.logger.error(f"Error processing sequence {seqhash_id}: {e}")
                        results[seqhash_id] = {
                            'seqhash_id': seqhash_id,
                            'description': f'Error: {str(e)}'
                        }
                
                app.logger.debug(f"Retrieved sequence details: {results}")
                return results
                
        except Exception as e:
            app.logger.error(f"Database error: {e}")
            return {}
    return app


import pandas as pd

def run_mmseqs2_search(sequence, db_path):
    """Runs MMSeqs2 search and returns parsed results as a DataFrame."""
    current_app.logger.debug(f"Running MMSeqs2 search with sequence: {sequence}")

    with tempfile.TemporaryDirectory() as temp_dir:
        input_file = os.path.join(temp_dir, "input.fasta")
        output_file = os.path.join(temp_dir, "output.tsv")
        tmp_dir = os.path.join(temp_dir, "tmp")

        with open(input_file, 'w') as f:
            f.write(f">query\n{sequence}\n")

        mmseqs_command = [
            "mmseqs", "easy-search", input_file, db_path, output_file, tmp_dir,
            "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tseq"
        ]

        process = subprocess.run(mmseqs_command, capture_output=True, text=True)
        if process.returncode != 0:
            current_app.logger.error(f"MMSeqs2 Error: {process.stderr}")
            return pd.DataFrame(), f"MMSeqs2 failed: {process.stderr}"

        if not os.path.exists(output_file):
            return pd.DataFrame(), "MMSeqs2 did not produce an output file"

        df = pd.read_csv(output_file, sep='\t', names=[
            'query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen', 
            'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'tseq'
        ])
        
    return df, None

def fetch_annotations_and_clusters(seqhash_ids, db_path):
    """Fetch annotations and cluster info from DuckDB."""
    if not seqhash_ids:
        return pd.DataFrame(), pd.DataFrame()  # return empty dfs if no input

    try:
        with duckdb.connect(db_path) as db:
            # fetch annotations
            annotations_query = f"""
                SELECT s.seqhash_id AS target, s.sample_id, m.organism, 
                       a.description, a.cog_category, a.preferred_name
                FROM sequences s
                JOIN sra_metadata m ON s.sample_id = m.sample_id
                LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
                WHERE s.seqhash_id IN ({','.join([f"'{x}'" for x in seqhash_ids])})
            """
            annotations_df = db.execute(annotations_query).fetchdf()

            # fetch clusters
            cluster_query = f"""
                SELECT ANY_VALUE(cm.seqhash_id) AS target, 
                    cm.cluster_id, 
                    GROUP_CONCAT(cm.seqhash_id, ';') AS cluster_members, 
                    COUNT(*) AS cluster_size
                FROM cluster_members cm
                WHERE cm.seqhash_id IN ({','.join([f"'{x}'" for x in seqhash_ids])})
                GROUP BY cm.cluster_id
            """
            cluster_df = db.execute(cluster_query).fetchdf()

        return annotations_df, cluster_df

    except Exception as e:
        current_app.logger.error(f"Database error: {e}")
        return pd.DataFrame(), pd.DataFrame()

import pandas as pd

def merge_results_with_metadata(parsed_results, annotations, cluster_info):
    """Merges MMSeqs2 results with annotations, organism/sample_id, and cluster data."""
    if not parsed_results:
        return []

    mmseqs_df = pd.DataFrame(parsed_results)  # Convert MMSeqs2 results into a DataFrame

    if annotations.empty:
        annotations = pd.DataFrame(columns=['target', 'organism', 'sample_id', 'description', 'cog_category', 'preferred_name'])
    if cluster_info.empty:
        cluster_info = pd.DataFrame(columns=['target', 'cluster_size', 'cluster_members'])

    # Merge annotations
    # import ipdb; ipdb.set_trace()
    merged_df = mmseqs_df.merge(annotations, on='target', how='left')

    # Merge cluster information
    merged_df = merged_df.merge(cluster_info, on='target', how='left')

    # Fill missing values
    merged_df.fillna({
        'organism': 'Unknown',
        'sample_id': 'Unknown',
        'description': 'No annotation',
        'cog_category': None,
        'preferred_name': None,
        'cluster_size': 1,
        'cluster_members': ''
    }, inplace=True)

    return merged_df.to_dict(orient='records')  # Convert back to list of dictionaries

def process_search_request(sequence, fasta_path, duckdb_path):
    """Handles the entire search process and merges MMSeqs2 results with metadata."""
    mmseqs_df, error = run_mmseqs2_search(sequence, fasta_path)
    if error:
        return {'headers': ['Error'], 'data': [{'Error': error}]}

    if mmseqs_df.empty:
        return {'headers': ['Error'], 'data': [{'Error': "No results found"}]}

    # Debugging breakpoint to inspect mmseqs_df
    # import ipdb; ipdb.set_trace()

    # Get sequence ids from mmseqs results
    seqhash_ids = mmseqs_df['target'].unique().tolist()

    # Fetch annotations and clusters from duckdb
    annotations_df, cluster_df = fetch_annotations_and_clusters(seqhash_ids, duckdb_path)

    # Debugging breakpoint before merge
    # ipdb.set_trace()

    # Merge annotations
    merged_df = mmseqs_df.merge(annotations_df, on='target', how='left')

    # Merge clusters
    merged_df = merged_df.merge(cluster_df, on='target', how='left').dropna()

    # Fill missing values
    # merged_df.fillna({
    #     'organism': 'Unknown',
    #     'sample_id': 'Unknown',
    #     'description': 'No annotation',
    #     'cog_category': None,
    #     'preferred_name': None,
    #     'cluster_size': 1,
    #     'cluster_members': ''
    # }, inplace=True)

    desired_headers = [
        'query', 'organism', 'sample_id', 'preferred_name', 'target',  'description', 'tseq',
        'cluster_members', 'cog_category',
        'pident', 'alnlen', 'mismatch', 'gapopen', 
        'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 
        'cluster_size'
    ]
    
    
    headers = [h for h in desired_headers if h in merged_df.columns]
    return {'headers': headers, 'data': merged_df.to_dict(orient='records')}


if __name__ == '__main__':
    app = create_app('development')
    app.run(host='0.0.0.0', port=8888, debug=True)
