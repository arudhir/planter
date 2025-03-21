from flask import Flask, render_template, request, jsonify, current_app, redirect, url_for, flash
import subprocess
import tempfile
import os
import logging
import shutil
import json
import uuid
import threading
import re
from datetime import datetime
from pathlib import Path
import pandas as pd
import duckdb
from app.config import config
from planter.database.query_manager import DatabaseManager

# Global storage for pipeline jobs
# In a production app, this would be in a database
pipeline_jobs = {}

def create_app(config_name='default'):
    """Creates and configures the Flask app."""
    app = Flask(__name__)
    app.config.from_object(config[config_name])
    app.secret_key = os.environ.get('SECRET_KEY', 'dev-key-for-session-flashes')
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
            
    @app.route('/download_database', methods=['POST'])
    def download_database():
        """Downloads the master.duckdb file from S3."""
        try:
            s3_path = f"s3://{app.config['S3_BUCKET']}/{app.config['S3_DB_KEY']}"
            output_path = app.config['DUCKDB_PATH']
            
            # Create directory if it doesn't exist
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            
            # Run AWS CLI command to download the file with progress
            cmd = ['aws', 's3', 'cp', s3_path, output_path, '--no-progress']
            app.logger.info(f"Running command: {' '.join(cmd)}")
            
            # Actually start the download process
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            
            # Log that we've started the download
            app.logger.info(f"Started download process with PID: {process.pid}")
            
            return jsonify({
                'status': 'started', 
                'message': f"Download started from {s3_path}",
                'output_path': output_path,
                's3_path': s3_path
            })
            
        except Exception as e:
            app.logger.error(f"Download database error: {str(e)}")
            return jsonify({'status': 'error', 'message': str(e)}), 500
            
    @app.route('/download_progress', methods=['GET'])
    def download_progress():
        """Checks the download progress and status."""
        try:
            # Get source and destination from the query parameters
            s3_path = request.args.get('s3_path')
            output_path = request.args.get('output_path')
            
            if not s3_path or not output_path:
                return jsonify({'status': 'error', 'message': 'Missing parameters'}), 400
                
            # If output file exists, get its current size
            if os.path.exists(output_path):
                current_size = os.path.getsize(output_path)
                
                # Get the total size of the file on S3
                cmd = ['aws', 's3api', 'head-object', 
                      '--bucket', app.config['S3_BUCKET'], 
                      '--key', app.config['S3_DB_KEY']]
                      
                process = subprocess.run(cmd, capture_output=True, text=True)
                
                if process.returncode == 0:
                    import json
                    file_info = json.loads(process.stdout)
                    total_size = file_info.get('ContentLength', 0)
                    
                    # Calculate progress
                    if total_size > 0:
                        progress = min(100, int((current_size / total_size) * 100))
                        
                        # Check if download is complete
                        if current_size >= total_size:
                            return jsonify({
                                'status': 'success',
                                'progress': 100,
                                'message': f"Download complete. File saved to {output_path}",
                                'current_size': current_size,
                                'total_size': total_size,
                                'human_size': f"{current_size / (1024*1024):.2f} MB / {total_size / (1024*1024):.2f} MB"
                            })
                        else:
                            # Still downloading
                            return jsonify({
                                'status': 'downloading',
                                'progress': progress,
                                'message': f"Downloading: {progress}% complete",
                                'current_size': current_size,
                                'total_size': total_size,
                                'human_size': f"{current_size / (1024*1024):.2f} MB / {total_size / (1024*1024):.2f} MB"
                            })
                    else:
                        return jsonify({
                            'status': 'unknown',
                            'message': "Could not determine file size",
                            'current_size': current_size
                        })
                else:
                    # Couldn't get file info from S3, but file is downloading
                    return jsonify({
                        'status': 'downloading',
                        'message': f"Downloading (size unknown)",
                        'current_size': current_size
                    })
            else:
                # Start the download if it hasn't started yet
                cmd = ['aws', 's3', 'cp', s3_path, output_path]
                process = subprocess.Popen(cmd, 
                                          stdout=subprocess.PIPE, 
                                          stderr=subprocess.PIPE,
                                          text=True)
                
                return jsonify({
                    'status': 'started',
                    'message': "Download started",
                    'progress': 0
                })
                
        except Exception as e:
            app.logger.error(f"Download progress check error: {str(e)}")
            return jsonify({'status': 'error', 'message': str(e)}), 500
    
    @app.route('/create_refseq', methods=['POST'])
    def create_refseq():
        """Creates the reference FASTA file from the DuckDB database."""
        try:
            db_path = app.config['DUCKDB_PATH']
            output_path = app.config['REPSEQ_FASTA']
            
            # Verify database file exists
            if not os.path.exists(db_path):
                return jsonify({
                    'status': 'error',
                    'message': f"Database file not found at {db_path}. Please download it first."
                }), 404
            
            # Create directory if it doesn't exist
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            
            # Connect to DuckDB and extract representative sequences
            try:
                app.logger.info(f"Extracting representative sequences from {db_path} to {output_path}")
                with duckdb.connect(db_path) as conn:
                    query = """
                    SELECT
                       '>' || seqhash_id || chr(10) || sequence
                    FROM
                        sequences
                    WHERE
                        is_representative = TRUE;
                    """
                    result = conn.execute(query).fetchall()
                
                # Write the results to the output file
                with open(output_path, 'w') as f:
                    for row in result:
                        f.write(row[0] + '\n')
                
                return jsonify({
                    'status': 'success',
                    'message': f"Reference FASTA file created at {output_path}",
                    'path': output_path,
                    'count': len(result)
                })
                
            except Exception as e:
                app.logger.error(f"DuckDB query error: {str(e)}")
                return jsonify({'status': 'error', 'message': f"Database query error: {str(e)}"}), 500
            
        except Exception as e:
            app.logger.error(f"Create refseq error: {str(e)}")
            return jsonify({'status': 'error', 'message': str(e)}), 500
            
    @app.route('/get_preset_query', methods=['GET'])
    def get_preset_query():
        """Gets a preset SQL query template."""
        query_name = request.args.get('query_name')
        app.logger.debug(f"Requested query template: {query_name}")
        
        if not query_name:
            return jsonify({'status': 'error', 'message': 'No query name provided'}), 400
            
        # Load the SQL query from the database directory
        try:
            query_file = f"/home/ubuntu/planter/planter/database/queries/sql/{query_name}.sql"
            app.logger.debug(f"Looking for query file at: {query_file}")
            
            if not os.path.exists(query_file):
                app.logger.error(f"Query file not found: {query_file}")
                return jsonify({'status': 'error', 'message': f"Query '{query_name}' not found"}), 404
                
            with open(query_file, 'r') as f:
                query_sql = f.read()
                
            # Log placeholder check
            if query_name == 'search_samples_by_id':
                if 'SAMPLE_IDS_PLACEHOLDER' in query_sql:
                    app.logger.debug("Found SAMPLE_IDS_PLACEHOLDER in query template")
                else:
                    app.logger.error("SAMPLE_IDS_PLACEHOLDER NOT found in query template!")
            
            if query_name == 'search_sequence_by_seqhash':
                if 'SEQHASH_IDS_PLACEHOLDER' in query_sql:
                    app.logger.debug("Found SEQHASH_IDS_PLACEHOLDER in query template")
                else:
                    app.logger.error("SEQHASH_IDS_PLACEHOLDER NOT found in query template!")
                
            return jsonify({
                'status': 'success',
                'query': query_sql,
                'file_path': query_file
            })
        except Exception as e:
            app.logger.error(f"Get preset query error: {str(e)}")
            return jsonify({'status': 'error', 'message': str(e)}), 500
            
    @app.route('/execute_query', methods=['POST'])
    def execute_query():
        """Executes a SQL query against the DuckDB database."""
        try:
            # Get the query from the request
            data = request.json
            if not data or 'query' not in data:
                return jsonify({'status': 'error', 'message': 'No query provided'}), 400
                
            query = data['query'].strip()
            if not query:
                return jsonify({'status': 'error', 'message': 'Empty query'}), 400
                
            # Check if the database exists
            db_path = app.config['DUCKDB_PATH']
            if not os.path.exists(db_path):
                return jsonify({
                    'status': 'error',
                    'message': f"Database file not found at {db_path}. Please download it first."
                }), 404
                
            # Block unsafe operations
            unsafe_operations = ['CREATE', 'DROP', 'ALTER', 'INSERT', 'UPDATE', 'DELETE', 'TRUNCATE', 'VACUUM']
            for operation in unsafe_operations:
                if operation in query.upper():
                    return jsonify({
                        'status': 'error',
                        'message': f"Unsafe operation '{operation}' is not allowed"
                    }), 403
                    
            # Execute the query with a timeout
            app.logger.info(f"Executing query: {query}")
            
            try:
                # Count number of parameter placeholders (?)
                param_count = query.count('?')
                
                # Replace Jinja2-style parameters like {{ limit }} with actual values
                if '{{ limit }}' in query:
                    query = query.replace('{{ limit }}', '100')  # Default limit to 100
                
                # Create a list of NULL parameters if there are any placeholders
                params = [None] * param_count if param_count > 0 else []
                
                # Limited to 10 seconds execution time
                with duckdb.connect(db_path) as conn:
                    if params:
                        app.logger.info(f"Executing with {len(params)} NULL parameters")
                        result = conn.execute(query, params).fetchdf()
                    else:
                        result = conn.execute(query).fetchdf()
                    
                # Convert the result to a list of dictionaries
                results = result.to_dict(orient='records')
                columns = result.columns.tolist()
                
                # Apply a limit on the number of results to avoid overwhelming the UI
                max_results = 1000
                if len(results) > max_results:
                    results = results[:max_results]
                    truncated = True
                else:
                    truncated = False
                    
                response = {
                    'status': 'success',
                    'results': results,
                    'columns': columns,
                    'truncated': truncated,
                    'row_count': len(results),
                    'total_rows': len(result) if truncated else len(results)
                }
                
                return jsonify(response)
                
            except Exception as e:
                app.logger.error(f"Query execution error: {str(e)}")
                return jsonify({'status': 'error', 'message': f"Query execution failed: {str(e)}"}), 500
                
        except Exception as e:
            app.logger.error(f"Execute query error: {str(e)}")
            return jsonify({'status': 'error', 'message': str(e)}), 500
            
    @app.route('/run_pipeline', methods=['POST'])
    def run_pipeline():
        """Runs the Planter pipeline with the specified samples."""
        global pipeline_jobs  # Access the global variable
        
        try:
            # Parse JSON data from request
            data = request.json
            if not data:
                return jsonify({'status': 'error', 'message': 'No data provided'}), 400
                
            # Extract parameters
            cores = data.get('cores', 16)
            outdir = data.get('outdir', 'outputs')
            s3_bucket = data.get('s3_bucket', 'recombia.planter')
            dry_run = data.get('dry_run', False)
            background = data.get('background', False)
            
            # Parse samples
            samples_json = data.get('samples', '[]')
            try:
                # If samples_json is already a list, use it directly
                if isinstance(samples_json, list):
                    samples = samples_json
                else:
                    # Otherwise, parse it as JSON
                    samples = json.loads(samples_json)
            except json.JSONDecodeError:
                return jsonify({'status': 'error', 'message': 'Invalid samples format'}), 400
                
            if not samples:
                return jsonify({'status': 'error', 'message': 'No samples provided'}), 400
                
            # Generate a unique job ID
            job_id = str(uuid.uuid4())
            
            # Create job record
            pipeline_jobs[job_id] = {
                'id': job_id,
                'status': 'pending',
                'start_time': datetime.now().isoformat(),
                'cores': cores,
                'outdir': outdir,
                's3_bucket': s3_bucket,
                'samples': samples,
                'progress': 0,
                'logs': ['Pipeline job created.'],
                'process': None,
                'dry_run': dry_run,
                'background': background
            }
            
            # Start pipeline in a background thread if requested
            thread = threading.Thread(
                target=run_pipeline_job,
                args=(app, job_id, cores, outdir, s3_bucket, samples, pipeline_jobs, dry_run)
            )
            thread.daemon = True
            thread.start()
            
            message = 'Pipeline job started'
            if dry_run:
                message += ' (dry run)'
            if background:
                message += ' in background'
            
            return jsonify({
                'status': 'started',
                'job_id': job_id,
                'message': message,
                'background': background
            })
            
        except Exception as e:
            app.logger.error(f"Run pipeline error: {str(e)}")
            return jsonify({'status': 'error', 'message': str(e)}), 500
            
    @app.route('/pipeline_status/<job_id>', methods=['GET'])
    def pipeline_status(job_id):
        """Gets the status of a running pipeline job."""
        global pipeline_jobs  # Access the global variable
        
        try:
            if job_id not in pipeline_jobs:
                return jsonify({'status': 'error', 'message': 'Job not found'}), 404
                
            job = pipeline_jobs[job_id]
            
            # Check if the process has completed
            process = job.get('process')
            if process and job['status'] == 'running':
                # For background jobs, capture and process output when status is checked
                if job.get('background', False):
                    # Check if we need to process output
                    if not job.get('output_processed', False):
                        # Read any available stdout
                        progress_pattern = re.compile(r'(\d+)% complete')
                        
                        # Non-blocking read from stdout
                        while process.stdout.readable():
                            line = process.stdout.readline()
                            if not line:
                                break
                                
                            line = line.strip()
                            if line:
                                job['logs'].append(line)
                                
                                # Check for progress information
                                progress_match = progress_pattern.search(line)
                                if progress_match:
                                    progress = int(progress_match.group(1))
                                    job['progress'] = progress
                                    job['message'] = f"Running: {progress}% complete"
                
                # Check if process has completed
                if process.poll() is not None:
                    # Process has finished
                    returncode = process.returncode
                    
                    # Process any remaining output 
                    if process.stdout:
                        for line in process.stdout.readlines():
                            if line.strip():
                                job['logs'].append(line.strip())
                    
                    # Process any stderr
                    if process.stderr:
                        stderr = process.stderr.read().decode('utf-8') if process.stderr else ''
                        if stderr:
                            for line in stderr.splitlines():
                                if line.strip():
                                    job['logs'].append(f"ERROR: {line.strip()}")
                    
                    if returncode == 0:
                        job['status'] = 'completed'
                        job['progress'] = 100
                        job['message'] = "Pipeline completed successfully"
                        if job.get('dry_run', False):
                            job['message'] += " (dry run)"
                        job['logs'].append('Pipeline completed successfully.')
                    else:
                        job['status'] = 'failed'
                        job['logs'].append(f'Pipeline failed with return code {returncode}.')
                    
                    # Mark output as fully processed
                    job['output_processed'] = True
            
            # Get the latest logs
            logs_since_last_check = []
            if job.get('last_log_index', 0) < len(job['logs']):
                logs_since_last_check = job['logs'][job.get('last_log_index', 0):]
                job['last_log_index'] = len(job['logs'])
            
            return jsonify({
                'status': job['status'],
                'progress': job['progress'],
                'message': job.get('message', ''),
                'logs': logs_since_last_check,
                'background': job.get('background', False),
                'dry_run': job.get('dry_run', False)
            })
            
        except Exception as e:
            app.logger.error(f"Pipeline status error: {str(e)}")
            return jsonify({'status': 'error', 'message': str(e)}), 500

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

def run_pipeline_job(app, job_id, cores, outdir, s3_bucket, samples, jobs_dict, dry_run=False):
    """Background function to run the Snakemake pipeline as a subprocess."""
    with app.app_context():
        try:
            # Update job status
            job = jobs_dict[job_id]
            job['status'] = 'running'
            job['progress'] = 5
            job['logs'].append(f"Starting pipeline with {len(samples)} samples...")
            
            # Format the samples list for the command
            samples_json = json.dumps(samples)
            
            # Create the command to run
            cmd = [
                'docker-compose', 'run', '--rm', 'planter', 
                'snakemake', '--snakefile', 'planter/workflow/Snakefile',
                '--cores', str(cores),
                '--rerun-incomplete'
            ]
            
            # Add dry run flag if requested
            if dry_run:
                cmd.append('--dry-run')
                job['logs'].append("Running in dry-run mode (will not execute tasks)")
            
            # Add config parameters
            cmd.extend([
                '--config', f"outdir={outdir}", f"s3_bucket={s3_bucket}", f"samples={samples_json}"
            ])
            
            # Log the command being run
            cmd_str = ' '.join(cmd)
            app.logger.info(f"Running pipeline command: {cmd_str}")
            job['logs'].append(f"Executing: {cmd_str}")
            
            # Start the process
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                bufsize=1  # Line buffered
            )
            
            # Store the process in the job record
            job['process'] = process
            
            # Set up process handling based on background mode
            if job.get('background', False):
                # For background jobs, don't capture output immediately
                # We'll read it when the user checks the status
                job['logs'].append("Running in background mode")
                job['message'] = "Pipeline running in background"
                return
            
            # For foreground jobs, process output in real-time
            # Read output in real-time to update progress
            progress_pattern = re.compile(r'(\d+)% complete')
            
            # Process stdout
            for line in iter(process.stdout.readline, ''):
                line = line.strip()
                if line:
                    job['logs'].append(line)
                    
                    # Check for progress information
                    progress_match = progress_pattern.search(line)
                    if progress_match:
                        progress = int(progress_match.group(1))
                        job['progress'] = progress
                        job['message'] = f"Running: {progress}% complete"
            
            # Wait for the process to complete
            process.wait()
            
            # Process any remaining stderr
            stderr = process.stderr.read()
            if stderr:
                for line in stderr.splitlines():
                    if line.strip():
                        job['logs'].append(f"ERROR: {line.strip()}")
            
            # Update final status
            if process.returncode == 0:
                job['status'] = 'completed'
                job['progress'] = 100
                job['message'] = "Pipeline completed successfully"
                if dry_run:
                    job['message'] += " (dry run)"
                job['logs'].append("Pipeline execution completed successfully.")
            else:
                job['status'] = 'failed'
                job['message'] = f"Pipeline failed with exit code {process.returncode}"
                job['logs'].append(f"Pipeline execution failed with exit code {process.returncode}")
            
        except Exception as e:
            app.logger.error(f"Pipeline job error: {str(e)}")
            if job_id in jobs_dict:
                jobs_dict[job_id]['status'] = 'failed'
                jobs_dict[job_id]['message'] = f"Error: {str(e)}"
                jobs_dict[job_id]['logs'].append(f"ERROR: {str(e)}")


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
