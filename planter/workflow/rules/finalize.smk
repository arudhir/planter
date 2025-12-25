import os
import logging
import subprocess
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
    extract_representative_sequences_v2,
    apply_clustering_schema,
    get_clustering_stats,
    create_duckdb,
    merge_duckdbs,
    validate_duckdb_schema,
    update_clusters
)
from planter.clustering import ClusteringManager, ClusteringParams
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
        master_db = Path(config['outdir']) / 'master.duckdb.cpy',
    run:
        shell(
            'aws s3 cp s3://recombia.planter/master.duckdb {output.master_db}'
        )

rule mmseqs_clustering:
    """
    Simplified clustering using the new immutable architecture.

    Instead of complex incremental updates, we:
    1. Merge all sample DBs into master first
    2. Export all sequences
    3. Run full clustering
    4. Store as a new immutable clustering version

    This is cleaner, more reliable, and easier to debug.
    """
    input:
        master_db = expand(rules.fetch_master_db.output),
        duckdb = expand(rules.create_duckdb.output, sample=config['samples']),
        proteins = expand(rules.transdecoder.output.longest_orfs_pep, sample=config['samples'])
    output:
        done = temp(Path(config['outdir']) / 'tmp/clustering_done.txt')
    params:
        tmp_dir = Path(config["outdir"]) / "tmp",
        schema_path = Path('/usr/src/planter/planter/database/schema/migrations/004_add_gene_protein_map.sql')
    log:
        path = Path(config['outdir']) / 'logs' / 'clustering.log'
    run:
        logger.info(f"Starting clustering at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info(f"Processing samples: {config['samples']}")

        try:
            # 1. Setup temp directory
            logger.info(f"Setting up temporary directory at {params.tmp_dir}")
            os.makedirs(params.tmp_dir, exist_ok=True)

            # 2. First merge all sample databases into master
            # This ensures all sequences are in the database before clustering
            logger.info("Merging sample databases into master database")
            start_time = time.time()

            master_db_path = merge_duckdbs(
                duckdb_paths=input.duckdb,
                master_db_path=str(input.master_db[0]),
                schema_sql_path=params.schema_path,
                upgrade_schema=True
            )
            merge_time = time.time() - start_time
            logger.info(f"Database merge completed in {merge_time:.2f} seconds")

            # 3. Apply the new immutable clustering schema
            logger.info("Applying immutable clustering schema")
            apply_clustering_schema(master_db_path)

            # 4. Run full clustering using the ClusteringManager
            logger.info("Running full sequence clustering")
            clustering_manager = ClusteringManager(master_db_path)

            clustering_params = ClusteringParams(
                min_seq_id=0.3,
                coverage=0.8,
                threads=workflow.cores
            )

            start_time = time.time()
            run_id, result = clustering_manager.run_full_clustering(
                params=clustering_params,
                work_dir=params.tmp_dir,
                notes=f"Samples: {', '.join(config['samples'])}"
            )
            clustering_time = time.time() - start_time

            if result.success:
                logger.info(f"Clustering completed in {clustering_time:.2f} seconds")
                logger.info(f"Run ID: {run_id}")
                logger.info(f"Sequences: {result.sequence_count}, Clusters: {result.cluster_count}")
            else:
                raise RuntimeError(f"Clustering failed: {result.error_message}")

            # 5. Mark as done
            shell(f"touch {output.done}")
            logger.info("Clustering completed successfully")

        except Exception as e:
            logger.error(f"Error in clustering: {str(e)}", exc_info=True)
            raise

rule update_database:
    """
    Validate and upload the updated database to S3.

    With the new architecture, merge and clustering are done in mmseqs_clustering,
    so this rule just validates and uploads.
    """
    input:
        master_db = rules.fetch_master_db.output,
        clustering_done = rules.mmseqs_clustering.output.done
    output:
        done = temp(Path(config['outdir']) / 'tmp/update_database_done.txt')
    params:
        canonical_db = "s3://recombia.planter/master.duckdb"
    log:
        path = Path(config['outdir']) / 'logs' / 'update_database.log'
    run:
        # Set up logging
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

        logger.info(f"Starting database validation and upload at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

        master_db_path = str(input.master_db[0])

        try:
            # Step 1: Validate the database
            logger.info("Validating database...")

            con = duckdb.connect(master_db_path)
            try:
                # Basic counts
                sequence_count = con.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
                sample_count = con.execute("SELECT COUNT(DISTINCT sample_id) FROM sra_metadata").fetchone()[0]

                # Check clustering state using new architecture
                stats = get_clustering_stats(master_db_path)

                logger.info(f"Database summary:")
                logger.info(f"  - Sequences: {sequence_count}")
                logger.info(f"  - Samples: {sample_count}")
                logger.info(f"  - Clustering architecture: {stats.get('architecture', 'unknown')}")
                logger.info(f"  - Clusters: {stats.get('cluster_count', 0)}")

                if stats.get('current_run_id'):
                    logger.info(f"  - Current clustering run: {stats['current_run_id']}")
                    logger.info(f"  - Avg cluster size: {stats.get('avg_cluster_size', 0):.2f}")

            finally:
                con.close()

            # Step 2: Upload to S3
            # FIXED: Upload the actual database, not the backup
            logger.info(f"Uploading updated database to S3 at {params.canonical_db}")
            shell(f"aws s3 cp {master_db_path} {params.canonical_db}")

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
