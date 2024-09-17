from flask import Flask, render_template, request, jsonify
import subprocess
import tempfile
import os
import logging
from config import config

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
                    "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
                ]
                
                process = subprocess.run(mmseqs_command, capture_output=True, text=True, check=True)
                
                app.logger.debug(f"MMSeqs2 stdout: {process.stdout}")
                app.logger.debug(f"MMSeqs2 stderr: {process.stderr}")
                
                if os.path.exists(output_file):
                    with open(output_file, 'r') as f:
                        results = f.read().splitlines()
                    
                    # Parse results into list of dictionaries
                    headers = ['query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits']
                    parsed_results = [dict(zip(headers, row.split('\t'))) for row in results]
                    
                    return {'headers': headers, 'data': parsed_results}
                else:
                    return {'headers': ['Error'], 'data': [{'Error': "MMSeqs2 did not produce output file"}]}
            
            except subprocess.CalledProcessError as e:
                app.logger.error(f"MMSeqs2 failed: {e}")
                return {'headers': ['Error'], 'data': [{'Error': f"MMSeqs2 search failed. Return code: {e.returncode}, Error: {e.stderr}"}]}
            except Exception as e:
                app.logger.error(f"Unexpected error: {e}")
                return {'headers': ['Error'], 'data': [{'Error': f"An unexpected error occurred: {str(e)}"}]}

    return app

if __name__ == '__main__':
    app = create_app('development')
    app.run(host='0.0.0.0', port=8888, debug=True)
