import os
from pathlib import Path
import shutil
import hashlib
import boto3
import botocore
import logging
import ipdb
import urllib3
import duckdb
from typing import List, Union

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
boto3.set_stream_logger(name='botocore', level=logging.INFO)

def create_zip_archive(output_dir):
    zip_archive = shutil.make_archive(
        base_name=f'{output_dir}', 
        root_dir=output_dir,
        format='zip'
    )
    return zip_archive


def upload_to_s3(output_dir, sample, bucket):
    """
    Upload all files from output_dir to `bucket` under an S3 prefix = sample.
    If a file already exists in S3, skip it (do not fail).
    Return True if all files were either skipped or uploaded successfully;
    return False if an unrecoverable error occurs.
    """
    s3 = boto3.client('s3')
    output_path = Path(output_dir)
    success = True  # We will set this to False only if we hit an actual error.
    
    for root, dirs, files in os.walk(output_dir):
        for filename in files:
            local_path = Path(root) / filename
            s3_key = str(Path(sample) / local_path.relative_to(output_path))
            
            try:
                # Check if the object already exists:
                s3.head_object(Bucket=bucket, Key=s3_key)
                # If no exception, the object exists => skip upload
                print(f"Skipping {local_path}: {s3_key} already exists in {bucket}.")
            except botocore.exceptions.ClientError as e:
                # If it's a 404 error => S3 object not found => proceed to upload
                if e.response['Error']['Code'] == "404":
                    print(f"Uploading {local_path} to s3://{bucket}/{s3_key}")
                    try:
                        s3.upload_file(str(local_path), bucket, s3_key)
                    except Exception as upload_error:
                        print(f"Error uploading {local_path} to {s3_key}: {upload_error}")
                        success = False
                else:
                    # Some other error (e.g., permissions) => fail this run
                    print(f"Error checking S3 object {s3_key}: {e}")
                    success = False
    
    return success


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

# We need to wait until all of the proteins are done being made, then concatenate them and update the repseq.
# This is a bit of a pain because we need to wait until all of the samples are done.
# We can do this by adding a rule that waits for all of the samples to be done, then concatenates the proteins and updates the repseq.
# rule concatenate_peptides:
#     input:
#         peptides = expand(rules.transdecoder.output.longest_orfs_pep, sample=config['samples']),
#     output:
#         all_peptides = Path(config['outdir']) / 'all_peptides.fasta',
#     run:
#         shell('cat {input.peptides} > {output.all_peptides}')

storage:
    provider = "s3",

def merge_duckdbs(
    duckdb_paths: List[Union[str, Path]],
    master_db_path: Union[str, Path],
    schema_sql_path: Union[str, Path]
) -> None:
    """
    Merge multiple DuckDB databases into a master DuckDB.
    
    Parameters:
      duckdb_paths (List[Union[str, Path]]): List of paths to source DuckDB files.
      master_db_path (Union[str, Path]): Path to the master (merged) DuckDB.
      schema_sql_path (Union[str, Path]): Path to the SQL file defining the schema.
    
    The function:
      - Creates (or opens) the master database.
      - Executes the schema SQL to create tables if they don't exist.
      - Iterates through each source database, attaches it,
        and inserts data into the master tables in dependency order.
      - Uses INSERT OR IGNORE to avoid duplicate key errors.
      - Detaches each source database after merging.
    """
    
    master_db_path = str(master_db_path)
    schema_sql_path = Path(schema_sql_path)
    
    # Read the schema SQL
    schema_sql = schema_sql_path.read_text()
    
    with duckdb.connect(master_db_path) as master_conn:
        # Set up the schema in the master database
        master_conn.execute(schema_sql)
        
        # Process each source DuckDB
        for i, source_db in enumerate(duckdb_paths):
            alias = f"db{i}"
            source_db_str = str(source_db)
            print(f"Attaching {source_db_str} as {alias}...")
            master_conn.execute(f"ATTACH '{source_db_str}' AS {alias};")
            
            # Insert data in dependency order
            master_conn.execute(f"""
                INSERT OR IGNORE INTO sra_metadata
                SELECT * FROM {alias}.sra_metadata;
            """)
            master_conn.execute(f"""
                INSERT OR IGNORE INTO sequences
                SELECT * FROM {alias}.sequences;
            """)
            master_conn.execute(f"""
                INSERT OR IGNORE INTO annotations
                SELECT * FROM {alias}.annotations;
            """)
            master_conn.execute(f"""
                INSERT OR IGNORE INTO go_terms
                SELECT * FROM {alias}.go_terms;
            """)
            master_conn.execute(f"""
                INSERT OR IGNORE INTO ec_numbers
                SELECT * FROM {alias}.ec_numbers;
            """)
            master_conn.execute(f"""
                INSERT OR IGNORE INTO kegg_info
                SELECT * FROM {alias}.kegg_info;
            """)
            master_conn.execute(f"""
                INSERT OR IGNORE INTO clusters
                SELECT * FROM {alias}.clusters;
            """)
            master_conn.execute(f"""
                INSERT OR IGNORE INTO cluster_members
                SELECT * FROM {alias}.cluster_members;
            """)
            
            # Optionally, merge schema_version if needed:
            # master_conn.execute(f"""
            #     INSERT OR IGNORE INTO schema_version
            #     SELECT * FROM {alias}.schema_version;
            # """)
            
            master_conn.execute(f"DETACH {alias};")
            print(f"Finished merging {source_db_str}\n")
        
        # Optional commit; DuckDB auto-commits by default.
        master_conn.commit()
    
    print("All databases have been merged into:", master_db_path)
    return master_db_path

def update_duckdb_with_cluster_info(master_db_path: str, tsv_path: str, inplace=False):
    """
    Update the 'sequences' table in the master DuckDB using a cluster mapping TSV.
    
    Parameters:
        master_db_path (str): Path to the master DuckDB file.
        tsv_path (str): Path to the TSV file containing cluster mappings.
                        The TSV should have two columns (no header):
                          - Column 0: repseq_id (the representative seqhash)
                          - Column 1: seqhash_id (the member seqhash)
                          
    The function updates the 'repseq_id' column in the 'sequences' table such that
    for every row where sequences.seqhash_id matches the TSV's member, repseq_id is set
    to the TSV's repseq_id. It also sets is_representative to TRUE for sequences that are representatives.
    """
    print(f"Updating master DB with cluster info from {tsv_path}")
    # Connect to the master DuckDB.
    con = duckdb.connect(master_db_path)
    
    # Create a temporary table by reading the TSV file.
    # We assume the TSV has no header and uses tab as the separator.
    con.execute(f"""
    CREATE TEMPORARY TABLE new_clusters AS 
      SELECT * FROM read_csv_auto('{tsv_path}', header=False, sep='\t', names=['repseq_id', 'seqhash_id'])
      AS (repseq_id VARCHAR, seqhash_id VARCHAR);
    """)
    
    # # Optionally, update the is_representative flag.
    # # For every sequence that appears as a repseq_id in the new_clusters table, mark it as representative.
    con.execute("""
    UPDATE sequences
    SET is_representative = TRUE
    WHERE seqhash_id IN (
        SELECT DISTINCT repseq_id FROM new_clusters
    );
    """)
    
    con.close()
    print("Master DB updated using cluster mapping.")



# rule batch_cluster_update:
#     """
#     Batch cluster update:
    
#     Concatenate all samples' peptide files, run mmseqs_cluster_update.py
#     to update the canonical repseq FASTA, upload it to S3, and create a
#     single done flag.
#     """
#     input:
#         # The current canonical repseq downloaded from S3.
#         old_repseq = storage.s3("s3://recombia.planter/repseq.faa"),
#         # Peptide files for all samples in the batch.
#         pep = expand("{outdir}/{sample}/transdecoder/{sample}.pep",
#                      outdir=config["outdir"],
#                      sample=config["samples"])
#     output:
#         # Final updated repseq file.
#         final_repseq = Path(config["outdir"]) / "cluster" / "repseq.faa",
#         # The new cluster TSV file (if your script produces it)
#         cluster_tsv  = Path(config["outdir"]) / "cluster" / "newClusterDB.tsv",
#         # A single done flag file.
#         done = Path(config["outdir"]) / "cluster" / "repseq_update_done.txt"
#     params:
#         # Use a dedicated temporary subdirectory inside the cluster folder.
#         cluster_dir = Path(config["outdir"]) / "cluster",
#         tmp_dir = Path(config["outdir"]) / "cluster" / "tmp",
#         # The canonical repseq S3 URI.
#         repseq_s3 = "s3://recombia.planter/repseq.faa"
#     threads: workflow.cores
#     run:
#         shell(
#             r"""
#             set -e  # Exit immediately if any command fails.
            
#             # Ensure the cluster output directory and temporary subdirectory exist.
#             mkdir -p {params.cluster_dir}
#             rm -rf {params.tmp_dir}
#             mkdir -p {params.tmp_dir}
            
#             # Concatenate all peptide files into a temporary aggregated file.
#             cat {input.pep} > {params.tmp_dir}/all_new.pep
            
#             # Run the mmseqs cluster update script.
#             # This script reads the current repseq (-i {input.old_repseq}),
#             # uses the aggregated peptides from {params.tmp_dir}/all_new.pep,
#             # and writes its updated repseq to {params.tmp_dir}/newRepSeqDB.fasta.
#             python ./planter/scripts/mmseqs_cluster_update.py -i {input.old_repseq} -o {params.tmp_dir} {params.tmp_dir}/all_new.pep
            
#             # Copy the updated repseq produced by mmseqs_cluster_update.py to the final output location.
#             cp {params.tmp_dir}/newRepSeqDB.fasta {output.final_repseq}
            
#             # Optionally, copy the cluster TSV file if your script produces it.
#             if [ -f {params.tmp_dir}/newClusterDB.tsv ]; then
#                 cp {params.tmp_dir}/newClusterDB.tsv {output.cluster_tsv}
#             fi
            
#             # Upload the final repseq to S3.
#             aws s3 cp {output.final_repseq} {params.repseq_s3}
            
#             # Clean up the temporary directory.
#             rm -rf {params.tmp_dir}
            
#             # Touch the done flag.
#             touch {output.done}
#             """
#         )

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
        proteins = expand(rules.transdecoder.output.longest_orfs_pep, sample=config['samples']),  # Redundant I suppose
        duckdb = expand(rules.create_duckdb.output, sample=config['samples']),
    output:
        updated_master = Path(config['outdir']) / 'updated_master.duckdb',
        concat_proteins = temp(Path(config['outdir']) / 'concat_proteins.pep'),
        # newClusterDB = temp(Path(config['outdir']) / 'newClusterDB.tsv'),
        done = temp(Path(config['outdir']) / 'update_database_done.txt')
    params:
        canonical_db = "s3://recombia.planter/master.duckdb",
        tmp_dir = Path(config["outdir"]) / "tmp"
    run:
        # import ipdb; ipdb.set_trace()
        shell(
                """
                rm -rf {params.tmp_dir}
                mkdir -p {params.tmp_dir}
                # Concatenate all of the proteins into a single file.
                cat {input.proteins} > {output.concat_proteins}

                # Extract repseq.faa from master.duckdb
                duckdb {input.master_db} -noheader -list -c "
                    SELECT 
                        '>' || seqhash_id || chr(10) || sequence
                    FROM 
                        sequences
                    WHERE 
                        repseq_id = seqhash_id;
                    " > {params.tmp_dir}/repseq.faa
                """
        )
        shell("""
            # Run mmseqs_cluster_update.py
            python ./planter/scripts/mmseqs_cluster_update.py -i {params.tmp_dir}/repseq.faa -o {params.tmp_dir} {output.concat_proteins}
            """
        )
        # Merge duckdbs
        master_db_path = merge_duckdbs(
            duckdb_paths=input.duckdb,
            master_db_path=input.master_db,
            schema_sql_path=Path('/usr/src/planter/planter/database/schema/migrations/001_initial_schema.sql')  # TODO: Move to config
        )
        # Update the concat_master.duckdb with the newClusterDB.tsv info
        update_duckdb_with_cluster_info(input.master_db, params.tmp_dir / "newClusterDB.tsv")
        print('Done updating master DB with cluster info.')

        # Copy the newly updated master DB to the output file
        shell("cp {input.master_db} {output.updated_master}")
        # Upload the updated master DB to S3
        shell(f"aws s3 cp {output.updated_master} {params.canonical_db}")
        
        # Touch the done flag.
        shell("touch {output.done}")

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
