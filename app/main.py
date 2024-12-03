from flask import Flask, render_template, request, jsonify
import subprocess
import tempfile
import os
import logging
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

    def get_sequence_annotations(seqhash_ids):
        """Fetch annotations for given seqhash_ids from the database."""
        try:
            with DatabaseManager(app.config['DUCKDB_PATH']) as db:
                # Pass seqhash_ids as a single-element tuple containing a list
                annotations = db.query_manager.sequence_annotations(values=(seqhash_ids,))
                app.logger.debug(f"Raw annotations: {annotations}")
                if annotations.empty:
                    app.logger.warning(f"No annotations found for seqhash_ids: {seqhash_ids}")
                    return {}
                # Return annotations as a dictionary keyed by seqhash_id
                return annotations.set_index('seqhash_id').to_dict('index')
        except Exception as e:
            app.logger.error(f"Database error: {e}")
            return {}

    def run_mmseqs2_easy_search(sequence):
        app.logger.debug(f"Starting MMSeqs2 search with sequence: {sequence}")
        
        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = os.path.join(temp_dir, "input.fasta")
            output_file = os.path.join(temp_dir, "output.tsv")
            tmp_dir = os.path.join(temp_dir, "tmp")
            
            with open(input_file, 'w') as f:
                f.write(f">query\n{sequence}\n")
            
            try:
                mmseqs_command = [
                    "mmseqs", "easy-search", input_file, app.config['MMSEQS_DB'], output_file, tmp_dir,
                    "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tseq"
                ]
                
                process = subprocess.run(mmseqs_command, capture_output=True, text=True, check=True)
                
                if os.path.exists(output_file):
                    with open(output_file, 'r') as f:
                        results = f.read().splitlines()
                    
                    # Parse MMSeqs2 results
                    headers = ['query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen', 
                             'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'tseq']
                    parsed_results = [dict(zip(headers, row.split('\t'))) for row in results]
                    
                    # Get target sequence IDs
                    seqhash_ids = [row['target'] for row in parsed_results]
                    app.logger.debug(f"Searching for seqhash_ids: {seqhash_ids}")
                    annotations = get_sequence_annotations(seqhash_ids)
                    app.logger.debug(f"Retrieved annotations: {annotations}")
                                        
                    # Get annotations from database
                    annotations = get_sequence_annotations(seqhash_ids)
                    
                    # Merge MMSeqs2 results with annotations
                    merged_results = []
                    for result in parsed_results:
                        target_id = result['target']
                        annotation = annotations.get(target_id, {})  # Safely get annotations or default to an empty dict
                        result.update(annotation)  # Add annotation fields to the result
                        merged_results.append(result)

                    # Dynamically update headers based on merged results
                    headers = list({key for row in merged_results for key in row.keys()})
                    
                    # Update headers with annotation fields
                    desired_headers = [
                        'query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen', 
                        'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits',
                        'organism', 'description', 'preferred_name', 'cog_category', 
                        'go_terms', 'ec_numbers', 'kegg_ko', 'kegg_pathway', 'tseq'
                    ]
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
