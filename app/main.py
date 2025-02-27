from flask import Flask, render_template, request, jsonify
import subprocess
import tempfile
import os
import logging
from pathlib import Path
from app.config import config
from planter.database.query_manager import DatabaseManager

def create_app(config_name='default'):
    app = Flask(__name__)
    app.config.from_object(config[config_name])
    logging.basicConfig(level=logging.DEBUG)

    @app.route('/', methods=['GET', 'POST'])
    def index():
        if request.method == 'POST':
            sequence = request.form['sequence']
            results = run_mmseqs2_easy_search(sequence)
            return render_template('results.html', results=results)
        return render_template('index.html')

    @app.route('/load_example', methods=['GET'])
    def load_example():
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
                annotations_df = db.query_manager.sequence.get_by_id(seqhash_ids[0])
                
                # Collect all sequence results
                results = {}
                
                for seqhash_id in seqhash_ids:
                    # Get detailed information for this sequence
                    annotation = db.query_manager.sequence.get_by_id(seqhash_id)
                    
                    if not annotation.empty:
                        # Get expression data for the corresponding gene
                        expression_data = db.query_manager.sequence.get_annotation_with_expression(
                            sample_id=annotation.iloc[0].get('sample_id', None)
                        )
                        
                        # Filter expression data for this sequence
                        seq_expression = expression_data[expression_data['protein_id'] == seqhash_id]
                        
                        # Create a dictionary for this sequence
                        seq_dict = annotation.iloc[0].to_dict()
                        
                        # Add expression data if available
                        if not seq_expression.empty:
                            seq_dict.update({
                                'gene_id': seq_expression.iloc[0].get('gene_id', ''),
                                'tpm': seq_expression.iloc[0].get('tpm', 0),
                                'num_reads': seq_expression.iloc[0].get('num_reads', 0)
                            })
                        
                        results[seqhash_id] = seq_dict
                
                app.logger.debug(f"Retrieved sequence details: {results}")
                return results
                
        except Exception as e:
            app.logger.error(f"Database error: {e}")
            return {}

    def run_mmseqs2_easy_search(sequence):
        app.logger.debug(f"Starting MMSeqs2 search with sequence: {sequence}")
        
        with tempfile.TemporaryDirectory() as temp_dir:
            tmp_dir = Path(temp_dir)
            input_file = tmp_dir / "input.fasta"
            output_file = tmp_dir / "output.tsv"
            mmseqs_tmp = tmp_dir / "tmp"
            
            with open(input_file, 'w') as f:
                f.write(f">query\n{sequence}\n")
            
            try:
                mmseqs_command = [
                    "mmseqs", "easy-search", str(input_file), app.config['MMSEQS_DB'], 
                    str(output_file), str(mmseqs_tmp),
                    "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tseq"
                ]
                
                subprocess.run(mmseqs_command, capture_output=True, text=True, check=True)
                
                if output_file.exists():
                    with open(output_file, 'r') as f:
                        results = f.read().splitlines()
                    
                    # Parse MMSeqs2 results
                    headers = ['query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen', 
                             'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'tseq']
                    parsed_results = [dict(zip(headers, row.split('\t'))) for row in results]
                    
                    # Get target sequence IDs
                    seqhash_ids = [row['target'] for row in parsed_results]
                    app.logger.debug(f"Searching for seqhash_ids: {seqhash_ids}")
                    
                    # Get comprehensive sequence information using the canonical query
                    sequence_details = get_sequence_details(seqhash_ids)
                    
                    # Merge MMSeqs2 results with detailed sequence information
                    merged_results = []
                    for result in parsed_results:
                        target_id = result['target']
                        details = sequence_details.get(target_id, {})  # Get sequence details or default to empty dict
                        result.update(details)  # Add all sequence details to the result
                        merged_results.append(result)

                    # Define desired display order for headers
                    desired_headers = [
                        'query', 'target', 'pident', 'alnlen', 'evalue', 'bits',
                        'preferred_name', 'description', 'organism', 'tpm', 'gene_id', 
                        'cog_category', 'go_terms', 'ec_numbers', 'kegg_pathway', 'kegg_ko',
                        'qstart', 'qend', 'tstart', 'tend', 'mismatch', 'gapopen', 'tseq'
                    ]
                    
                    # Filter to headers that are actually present in the results
                    headers = [h for h in desired_headers if any(h in d for d in merged_results)]
                    return {'headers': headers, 'data': merged_results}
                else:
                    return {'headers': ['Error'], 'data': [{'Error': "MMSeqs2 did not produce output file"}]}
                    
            except Exception as e:
                app.logger.error(f"Error: {e}")
                return {'headers': ['Error'], 'data': [{'Error': str(e)}]}

    return app

if __name__ == '__main__':
    app = create_app('development')
    app.run(host='0.0.0.0', port=8888, debug=True)
