import os
from pathlib import Path
import shutil
import hashlib
import boto3
import botocore
import logging
import ipdb
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
rule concatenate_peptides:
    input:
        peptides = expand(rules.transdecoder.output.longest_orfs_pep, sample=config['samples']),
    output:
        all_peptides = Path(config['outdir']) / 'all_peptides.fasta',
    run:
        shell('cat {input.peptides} > {output.all_peptides}')

storage:
    provider = "s3",
rule cluster_update:
    input:
        old_reps = storage.s3("s3://recombia.planter/repseq.faa"),
        pep = rules.concatenate_peptides.output.all_peptides
    output:
        new_reps = Path(config['outdir']) / 'cluster/new_sequences.fasta',
        removed_reps = Path(config['outdir']) / 'cluster/removed_sequences.fasta',
        updated_repseq = Path(config['outdir']) / 'cluster/newRepSeqDB.fasta',
        cluster_tsv = Path(config['outdir']) / 'cluster/newClusterDB.tsv',
    params:
        outdir = lambda wildcards: Path(config['outdir']) / 'cluster'
    shell:
        """
        ./planter/scripts/mmseqs_cluster_update.py \
            -i {input.old_reps} \
            -p {input.pep} \
            -o {params.outdir} \
            --output-new {output.new_reps} \
            --output-removed {output.removed_reps}
        """


# This doesn't handle the case where several samples are run in parallel.
# We need to track per-sample updates or do it in batch.
rule update_repseq:
    input:
        updated_repseq = rules.cluster_update.output.updated_repseq
    output:
        done_flag = Path(config['outdir']) / 'cluster/repseq_update_done.txt'  # ✅ Only track per-sample updates
    params:
        repseq_faa = storage.s3("s3://recombia.planter/repseq.faa")  # ✅ Pass the S3 path as a param
    shell:
        """
        aws s3 cp {input.updated_repseq} {params.repseq_faa}
        touch {output.done_flag}
        """

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
