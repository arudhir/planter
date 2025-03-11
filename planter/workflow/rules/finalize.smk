import os
import logging
from pathlib import Path
import duckdb
from typing import List, Union
import pandas as pd
import time
import shutil
from datetime import datetime

from planter.database.utils.s3 import create_zip_archive, upload_to_s3
from planter.database.schema.schema_version import get_db_schema_version, ensure_compatibility
from planter.database.utils.duckdb_utils import (
    extract_representative_sequences, 
    create_duckdb,
    merge_duckdbs,
    validate_duckdb_schema,
    update_clusters
)
# Setup module logger
logger = logging.getLogger(__name__)

rule get_qc_stats:
    input:
        fastp = rules.fastp_raw.output.json,
        salmon_metadata = rules.quant.output.stats,
        quantsf = rules.quant.output.formatted_tsv,
        transcripts = rules.rename_headers.output.fasta,
        eggnog = rules.eggnog.output.annotations,
    output:
        qc_stats = Path(config['outdir']) / '{sample}/{sample}_stats.json',
    threads: workflow.cores        
    run:
        shell(
            './planter/scripts/get_qc_stats.py '
            '--sample {wildcards.sample} '
            '--fastp {input.fastp} '
            '--salmon_metadata {input.salmon_metadata} '
            '--transcripts {input.transcripts} '
            '--quantsf {input.quantsf} '
            '--eggnog {input.eggnog} '
            '--output_file {output.qc_stats}'
        )

rule create_duckdb:
    input:
        analyze_eggnog = expand(rules.analyze_eggnog.output, sample=config['samples']),
        quant = expand(rules.quant.output, sample=config['samples']),
    output:
        duckdb = Path(config['outdir']) / '{sample}/{sample}.duckdb'
    params:
        outdir = lambda wildcards: Path(config['outdir'])
    run:
        logger.info(f"Creating DuckDB for sample {wildcards.sample}")
        create_duckdb(sample_id=wildcards.sample, outdir=params.outdir, duckdb_out=output.duckdb)

rule fetch_master_db:
    output:
        master_db = temp(Path(config['outdir']) / 'master.duckdb.cpy'),
    run:
        shell(
            'aws s3 cp s3://recombia.planter/master.duckdb {output.master_db}'
        )

rule mmseqs_clustering:
    input:
        master_db = expand(rules.fetch_master_db.output),
        proteins = expand(rules.transdecoder.output.longest_orfs_pep, sample=config['samples'])
    output:
        concat_proteins = temp(Path(config['outdir']) / 'tmp/concat_proteins.pep'),
        repseq_path = temp(Path(config['outdir']) / 'tmp/repseq.faa'),
        cluster_file = Path(config['outdir']) / 'tmp/newClusterDB.tsv',
        done = temp(Path(config['outdir']) / 'tmp/mmseqs_clustering_done.txt')
    params:
        tmp_dir = Path(config["outdir"]) / "tmp"
    log:
        path = Path(config['outdir']) / 'logs' / 'mmseqs_clustering.log'
    run:
        # # Set up logging to both file and console
        # log_file = str(log.path)
        # os.makedirs(os.path.dirname(log_file), exist_ok=True)
        # logging.basicConfig(
        #     level=logging.INFO,
        #     format='%(asctime)s - %(levelname)s - %(message)s',
        #     handlers=[
        #         logging.FileHandler(log_file),
        #         logging.StreamHandler()
        #     ]
        # )
        # logger = logging.getLogger('mmseqs_clustering')
        
        logger.info(f"Starting MMSeqs clustering at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info(f"Processing samples: {config['samples']}")
        
        try:
            # 1. Setup temp directory
            logger.info(f"Setting up temporary directory at {params.tmp_dir}")
            if os.path.exists(params.tmp_dir):
                shutil.rmtree(params.tmp_dir)
            os.makedirs(params.tmp_dir, exist_ok=True)

            # 2. Concatenate all proteins
            logger.info("Concatenating protein sequences")
            with open(str(output.concat_proteins), 'wb') as concat_file:
                for protein_file in input.proteins:
                    with open(protein_file, 'rb') as pf:
                        shutil.copyfileobj(pf, concat_file)
            
            # 3. Extract representative sequences
            logger.info(f"Extracting representative sequences from {input.master_db}")
            repseq_path = extract_representative_sequences(
                db_path=input.master_db,
                output_path=output.repseq_path
            )
            
            # See how many sequences are in the repseq file - capture the output properly
            repseq_count = subprocess.check_output(f"grep -c '>' {output.repseq_path}", shell=True).decode('utf-8').strip()
            logger.info(f"Number of sequences in repseq file: {repseq_count}")

            # 4. Run MMSeqs2 clustering update - use exactly the same command format as original
            logger.info("Running MMSeqs2 clustering update")
            start_time = time.time()
            shell(
                f"""
                python ./planter/scripts/mmseqs_cluster_update.py -i {output.repseq_path} -o {params.tmp_dir} {output.concat_proteins}
                """
            )
            clustering_time = time.time() - start_time
            logger.info(f"MMSeqs2 clustering completed in {clustering_time:.2f} seconds")
            
            # Ensure the cluster file exists
            if not os.path.exists(output.cluster_file):
                raise FileNotFoundError(f"Expected cluster file not found at {output.cluster_file} after MMSeqs2 clustering")
            
            # 5. Mark as done
            logger.info("MMSeqs clustering completed successfully")
            shell(f"touch {output.done}")
            
        except Exception as e:
            logger.error(f"Error in MMSeqs clustering: {str(e)}", exc_info=True)
            raise

rule update_database:
    input:
        duckdb = expand(rules.create_duckdb.output, sample=config['samples']),
        master_db = expand(rules.fetch_master_db.output),
        cluster_file = rules.mmseqs_clustering.output.cluster_file
    output:
        done = temp(Path(config['outdir']) / 'tmp/update_database_done.txt')
    params:
        canonical_db = "s3://recombia.planter/master.duckdb",
        schema_path = Path('/usr/src/planter/planter/database/schema/migrations/004_add_gene_protein_map.sql'),
        # Target schema version (set to None to use latest available)
        target_schema_version = None
    log:
        path = Path(config['outdir']) / 'logs' / 'update_database.log'
    run:
        # Set up logging to both file and console
        log_file = str(log.path)
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        logger = logging.getLogger('update_database')
        
        logger.info(f"Starting database update at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        
        try:
            # Step 1: Merge all sample DuckDBs into master
            logger.info("Merging sample databases into master database")
            start_time = time.time()  # Move to decorator probably
            
            # Merge databases with schema compatibility
            master_db_path = merge_duckdbs(
                duckdb_paths=input.duckdb,
                master_db_path=input.master_db,
                schema_sql_path=params.schema_path,
                upgrade_schema=True,
                target_schema_version=params.target_schema_version
            )
            merge_time = time.time() - start_time
            logger.info(f"Database merge completed in {merge_time:.2f} seconds")
            
            # Step 2: Update with cluster information
            logger.info("Updating master database with cluster information")
            cluster_file = input.cluster_file
            if not os.path.exists(cluster_file):
                raise FileNotFoundError(f"Cluster file not found at {cluster_file}")
                
            start_time = time.time()
            # Use the update_clusters function with logging
            update_clusters(
                db_path=input.master_db, 
                tsv_path=cluster_file,
                backup_first=True  # Adds .backup to the end of the file, so master.duckdb.cpy --> master.duckdb.cpy.backup
            )
            update_time = time.time() - start_time
            logger.info(f"Cluster information update completed in {update_time:.2f} seconds")
            
            # Step 3: Validate the database
            logger.info("Validating database...")
            try:
                # Get basic statistics
                con = duckdb.connect(str(output.updated_master))
                sequence_count = con.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
                cluster_count = con.execute("SELECT COUNT(*) FROM clusters").fetchone()[0] 
                member_count = con.execute("SELECT COUNT(*) FROM cluster_members").fetchone()[0]
                sample_count = con.execute("SELECT COUNT(DISTINCT sample_id) FROM sra_metadata").fetchone()[0]
                
                logger.info(f"Database summary: {sequence_count} sequences, {cluster_count} clusters, {member_count} cluster members, {sample_count} samples")
                con.close()
            except Exception as e:
                logger.error(f"Error validating database: {str(e)}")
            
            
            # Step 4: Upload to S3 if needed
            logger.info(f"Uploading updated database to S3 at {params.canonical_db}")
            shell(f"aws s3 cp {input.master_db}.backup {params.canonical_db}")
            
            # Mark as done
            logger.info("Database update completed successfully")
            shell(f"touch {output.done}")
            
        except Exception as e:
            logger.error(f"Error updating database: {str(e)}", exc_info=True)
            raise



rule upload_to_s3:
    input:
        analyze_eggnog = expand(rules.analyze_eggnog.output, sample=config['samples']),
        quant = expand(rules.quant.output, sample=config['samples']),
        duckdb = expand(rules.create_duckdb.output, sample=config['samples']),
    output:
        done = expand(Path(config['outdir']) / '{sample}/{sample}_s3_upload.done', sample=config['samples'])
    run:
        samples = config['samples']
        if not isinstance(samples, list):
            samples = [samples]
        
        for sample in samples:
            print('Finalizing sample: ', sample)
            output_dir = Path(config['outdir']) / sample
            bucket = config['s3_bucket']

            success = upload_to_s3(output_dir, sample, bucket)
            # Only "touch" the .done file if everything was either skipped or successfully uploaded
            if success:
                (Path(config['outdir']) / f'{sample}/{sample}_s3_upload.done').touch()
            else:
                print(f"Encountered errors uploading {sample}; not creating .done file.")
