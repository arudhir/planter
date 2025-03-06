import os
import logging
from pathlib import Path
import duckdb
from typing import List, Union
import pandas as pd

from planter.database.utils.s3 import create_zip_archive, upload_to_s3
from planter.database.utils.duckdb_utils import merge_duckdbs, update_duckdb_with_cluster_info
from planter.database.schema.schema_version import get_db_schema_version, ensure_compatibility

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

storage:
    provider = "s3"




rule create_duckdb:
    input:
        analyze_eggnog = expand(rules.analyze_eggnog.output, sample=config['samples']),
        quant = expand(rules.quant.output, sample=config['samples']),
    output:
        duckdb = Path(config['outdir']) / '{sample}/{sample}.duckdb'
    params:
        outdir = lambda wildcards: Path(config['outdir'])
    run:
        shell(
            'python ./planter/scripts/create_duckdb.py '
            ' --sample_id {wildcards.sample} '
            ' --outdir {params.outdir} '
            ' --duckdb_out {output.duckdb}'
        )

rule update_database:
    input:
        master_db = storage.s3("s3://recombia.planter/master-database.duckdb"),
        proteins = expand(rules.transdecoder.output.longest_orfs_pep, sample=config['samples']),
        duckdb = expand(rules.create_duckdb.output, sample=config['samples']),
    output:
        updated_master = Path(config['outdir']) / 'updated_master.duckdb',
        concat_proteins = temp(Path(config['outdir']) / 'concat_proteins.pep'),
        done = temp(Path(config['outdir']) / 'update_database_done.txt')
    params:
        canonical_db = "s3://recombia.planter/master.duckdb",
        tmp_dir = Path(config["outdir"]) / "tmp",
        schema_path = Path('/usr/src/planter/planter/database/schema/migrations/001_initial_schema.sql'),
        # Target schema version (set to None to use latest available)
        target_schema_version = None
    log:
        Path(config['outdir']) / 'logs' / 'update_database.log'
    run:
        import logging
        import time
        import shutil
        from datetime import datetime
        
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
        logger.info(f"Processing samples: {config['samples']}")
        
        try:
            # Create a local copy of the master database for processing
            local_master = str(output.updated_master)
            logger.info(f"Creating local copy of master database at {local_master}")
            shutil.copy(str(input.master_db), local_master)
            
            # Check and ensure schema compatibility
            schema_version, was_upgraded = ensure_compatibility(
                local_master,
                required_version=params.target_schema_version
            )
            if was_upgraded:
                logger.info(f"Database schema was automatically upgraded to version {schema_version}")
            else:
                logger.info(f"Database schema version: {schema_version} (no upgrade needed)")
            
            # 1. Setup temp directory
            logger.info(f"Setting up temporary directory at {params.tmp_dir}")
            if os.path.exists(params.tmp_dir):
                shutil.rmtree(params.tmp_dir)
            os.makedirs(params.tmp_dir, exist_ok=True)
            
            # 2. Concatenate all proteins
            logger.info("Concatenating protein sequences")
            with open(output.concat_proteins, 'wb') as concat_file:
                for protein_file in input.proteins:
                    with open(protein_file, 'rb') as pf:
                        shutil.copyfileobj(pf, concat_file)
            
            # 3. Extract representative sequences from master DB
            logger.info("Extracting representative sequences from master database")
            repseq_path = params.tmp_dir / "repseq.faa"
            shell(
                f"""
                duckdb {local_master} -noheader -list -c "
                    SELECT 
                        '>' || seqhash_id || chr(10) || sequence
                    FROM 
                        sequences
                    WHERE 
                        repseq_id = seqhash_id;
                " > {repseq_path}
                """
            )
            
            # 4. Run MMSeqs2 clustering update
            logger.info("Running MMSeqs2 clustering update")
            start_time = time.time()
            shell(
                f"""
                python ./planter/scripts/mmseqs_cluster_update.py -i {repseq_path} -o {params.tmp_dir} {output.concat_proteins}
                """
            )
            clustering_time = time.time() - start_time
            logger.info(f"MMSeqs2 clustering completed in {clustering_time:.2f} seconds")
            
            # 5. Merge sample DuckDBs into master DB with schema compatibility
            logger.info("Merging sample databases into master database")
            start_time = time.time()
            
            # Check schema versions of sample databases
            sample_schema_versions = {}
            for sample_db in input.duckdb:
                sample_version = get_db_schema_version(sample_db)
                sample_schema_versions[str(sample_db)] = sample_version
                
            logger.info(f"Sample database schema versions: {sample_schema_versions}")
            
            # Merge databases with schema compatibility
            master_db_path = merge_duckdbs(
                duckdb_paths=input.duckdb,
                master_db_path=local_master,
                schema_sql_path=params.schema_path,
                upgrade_schema=True,
                target_schema_version=params.target_schema_version
            )
            merge_time = time.time() - start_time
            logger.info(f"Database merge completed in {merge_time:.2f} seconds")
            
            # 6. Update the master DB with new cluster information
            logger.info("Updating master database with new cluster information")
            cluster_file = params.tmp_dir / "newClusterDB.tsv"
            if not os.path.exists(cluster_file):
                raise FileNotFoundError(f"Cluster file not found at {cluster_file}")
                
            start_time = time.time()
            # Use schema compatibility features
            update_duckdb_with_cluster_info(
                local_master, 
                cluster_file,
                upgrade_schema=True
            )
            update_time = time.time() - start_time
            logger.info(f"Cluster information update completed in {update_time:.2f} seconds")
            
            # 7. Verify database integrity
            logger.info("Verifying database integrity")
            con = duckdb.connect(local_master)
            try:
                # Get schema version for reporting
                db_version = get_db_schema_version(local_master)
                logger.info(f"Final database schema version: {db_version}")
                
                # Run basic integrity checks
                tables = con.execute("SHOW TABLES").fetchall()
                logger.info(f"Found {len(tables)} tables in database")
                
                # Check for basic stats
                sequence_count = con.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
                cluster_count = con.execute("SELECT COUNT(*) FROM clusters").fetchone()[0]
                sample_count = con.execute("SELECT COUNT(DISTINCT sample_id) FROM sra_metadata").fetchone()[0]
                
                logger.info(f"Database summary: {sequence_count} sequences, {cluster_count} clusters, {sample_count} samples")
                
                # Verify that all sequences have a repseq_id
                null_repseqs = con.execute("SELECT COUNT(*) FROM sequences WHERE repseq_id IS NULL").fetchone()[0]
                if null_repseqs > 0:
                    logger.warning(f"Warning: {null_repseqs} sequences have NULL repseq_id")
                
                # Schema-specific checks
                if db_version >= 2:
                    # Check representative flag consistency in v2+ schema
                    inconsistent_reps = con.execute("""
                        SELECT COUNT(*) FROM sequences 
                        WHERE (is_representative = TRUE AND seqhash_id != repseq_id)
                        OR (is_representative = FALSE AND seqhash_id = repseq_id)
                    """).fetchone()[0]
                    
                    if inconsistent_reps > 0:
                        logger.warning(f"Warning: {inconsistent_reps} sequences have inconsistent representative flags")
            finally:
                con.close()
            
            # 8. Upload to S3
            logger.info(f"Uploading updated database to S3 at {params.canonical_db}")
            shell(f"aws s3 cp {local_master} {params.canonical_db}")
            
            # 9. Mark as done
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
