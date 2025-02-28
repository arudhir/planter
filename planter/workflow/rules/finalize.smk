import os
from pathlib import Path
import duckdb
from typing import List, Union
import pandas as pd

from planter.database.utils.s3 import create_zip_archive, upload_to_s3
from planter.database.utils.duckdb_utils import merge_duckdbs, update_duckdb_with_cluster_info


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
    provider = "s3",




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
