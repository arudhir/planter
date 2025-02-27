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
                # Check if MMSEQS_DB is accessible
                if not os.path.exists(app.config['MMSEQS_DB']):
                    app.logger.error(f"MMSeqs database not found at: {app.config['MMSEQS_DB']}")
                    return {'headers': ['Error'], 'data': [{'Error': f"MMSeqs database not found or inaccessible. Please check configuration."}]}
                
                mmseqs_command = [
                    "mmseqs", "easy-search", str(input_file), app.config['MMSEQS_DB'], 
                    str(output_file), str(mmseqs_tmp),
                    "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tseq"
                ]
                
                try:
                    process_result = subprocess.run(mmseqs_command, capture_output=True, text=True, check=True)
                    app.logger.debug(f"MMSeqs2 stdout: {process_result.stdout}")
                    app.logger.debug(f"MMSeqs2 stderr: {process_result.stderr}")
                except subprocess.CalledProcessError as e:
                    app.logger.error(f"MMSeqs2 search failed: {e}")
                    app.logger.error(f"Command output: {e.stdout}")
                    app.logger.error(f"Command error: {e.stderr}")
                    return {'headers': ['Error'], 'data': [{'Error': f"MMSeqs2 search failed: {str(e)}"}]}
                
                if output_file.exists():
                    try:
                        with open(output_file, 'r') as f:
                            results = f.read().splitlines()
                        
                        if not results:
                            app.logger.warning("MMSeqs2 returned empty results")
                            return {'headers': ['Message'], 'data': [{'Message': "No similar sequences found in the database."}]}
                        
                        # Parse MMSeqs2 results
                        headers = ['query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen', 
                                 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'tseq']
                        parsed_results = []
                        
                        for row in results:
                            try:
                                fields = row.split('\t')
                                if len(fields) != len(headers):
                                    app.logger.warning(f"Skipping malformed result row: {row}")
                                    continue
                                parsed_results.append(dict(zip(headers, fields)))
                            except Exception as row_error:
                                app.logger.error(f"Error parsing MMSeqs2 result row: {row_error}")
                                continue
                        
                        if not parsed_results:
                            return {'headers': ['Message'], 'data': [{'Message': "Could not parse any search results."}]}
                        
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
                            'preferred_name', 'description', 'organism', 
                            # 'tpm', 'gene_id',  # Expression data temporarily disabled
                            'cog_category', 'go_terms', 'ec_numbers', 'kegg_pathway', 'kegg_ko',
                            'qstart', 'qend', 'tstart', 'tend', 'mismatch', 'gapopen', 'tseq'
                        ]
                        
                        # Filter to headers that are actually present in the results
                        headers = [h for h in desired_headers if any(h in d for d in merged_results)]
                        return {'headers': headers, 'data': merged_results}
                    except Exception as parse_error:
                        app.logger.error(f"Error parsing MMSeqs2 results: {parse_error}")
                        return {'headers': ['Error'], 'data': [{'Error': f"Error parsing search results: {str(parse_error)}"}]}
                else:
                    return {'headers': ['Error'], 'data': [{'Error': "MMSeqs2 did not produce output file"}]}
                    
            except Exception as e:
                app.logger.error(f"Error: {e}")
                return {'headers': ['Error'], 'data': [{'Error': str(e)}]}

    return app

if __name__ == '__main__':
    app = create_app('development')
    app.run(host='0.0.0.0', port=8888, debug=True)
