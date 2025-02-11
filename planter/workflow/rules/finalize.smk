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
# rule concatenate_peptides:
#     input:
#         peptides = expand(rules.transdecoder.output.longest_orfs_pep, sample=config['samples']),
#     output:
#         all_peptides = Path(config['outdir']) / 'all_peptides.fasta',
#     run:
#         shell('cat {input.peptides} > {output.all_peptides}')

storage:
    provider = "s3",
storage:
    provider = "s3",

rule batch_cluster_update:
    """
    Batch cluster update:
    
    Concatenate all samples' peptide files, run mmseqs_cluster_update.py
    to update the canonical repseq FASTA, upload it to S3, and create a
    single done flag.
    """
    input:
        # The current canonical repseq downloaded from S3.
        old_repseq = storage.s3("s3://recombia.planter/repseq.faa"),
        # Peptide files for all samples in the batch.
        pep = expand("{outdir}/{sample}/transdecoder/{sample}.pep",
                     outdir=config["outdir"],
                     sample=config["samples"])
    output:
        # Final updated repseq file.
        final_repseq = Path(config["outdir"]) / "cluster" / "repseq.faa",
        # The new cluster TSV file (if your script produces it)
        cluster_tsv  = Path(config["outdir"]) / "cluster" / "newClusterDB.tsv",
        # A single done flag file.
        done = Path(config["outdir"]) / "cluster" / "repseq_update_done.txt"
    params:
        # Use a dedicated temporary subdirectory inside the cluster folder.
        cluster_dir = Path(config["outdir"]) / "cluster",
        tmp_dir = Path(config["outdir"]) / "cluster" / "tmp",
        # The canonical repseq S3 URI.
        repseq_s3 = "s3://recombia.planter/repseq.faa"
    threads: workflow.cores
    run:
        shell(
            r"""
            set -e  # Exit immediately if any command fails.
            
            # Ensure the cluster output directory and temporary subdirectory exist.
            mkdir -p {params.cluster_dir}
            rm -rf {params.tmp_dir}
            mkdir -p {params.tmp_dir}
            
            # Concatenate all peptide files into a temporary aggregated file.
            cat {input.pep} > {params.tmp_dir}/all_new.pep
            
            # Run the mmseqs cluster update script.
            # This script reads the current repseq (-i {input.old_repseq}),
            # uses the aggregated peptides from {params.tmp_dir}/all_new.pep,
            # and writes its updated repseq to {params.tmp_dir}/newRepSeqDB.fasta.
            python ./planter/scripts/mmseqs_cluster_update.py -i {input.old_repseq} -o {params.tmp_dir} {params.tmp_dir}/all_new.pep
            
            # Copy the updated repseq produced by mmseqs_cluster_update.py to the final output location.
            cp {params.tmp_dir}/newRepSeqDB.fasta {output.final_repseq}
            
            # Optionally, copy the cluster TSV file if your script produces it.
            if [ -f {params.tmp_dir}/newClusterDB.tsv ]; then
                cp {params.tmp_dir}/newClusterDB.tsv {output.cluster_tsv}
            fi
            
            # Upload the final repseq to S3.
            aws s3 cp {output.final_repseq} {params.repseq_s3}
            
            # Clean up the temporary directory.
            rm -rf {params.tmp_dir}
            
            # Touch the done flag.
            touch {output.done}
            """
        )
# rule cluster_update:
#     input:
#         old_reps = storage.s3("s3://recombia.planter/repseq.faa"),
#         pep = rules.transdecoder.output.longest_orfs_pep
#     output:
#         new_repseq = Path(config['outdir']) / '{sample}/cluster/newRepSeqDB.fasta',
#         cluster_tsv = Path(config['outdir']) / '{sample}/cluster/newClusterDB.tsv',
#     params:
#         outdir = lambda wildcards: Path(config['outdir']) / f'{wildcards.sample}/cluster'
#     threads: workflow.cores
#     run:
#         shell(
#             """
#             ./planter/scripts/mmseqs_cluster_update.py \
#                 -i {input.old_reps} \
#                 -o {params.outdir} \
#                 {input.pep}
#             """
#         )

# rule update_repseq:
#     input:
#         new_repseq = rules.cluster_update.output.new_repseq
#     output:
#         done_flag = Path(config['outdir']) / '{sample}/cluster/repseq_update_done.txt'
#     params:
#         repseq_faa = "s3://recombia.planter/repseq.faa"  # Use the S3 URI directly
#     resources:
#         repseq_update=1
#     shell:
#         """
#         aws s3 cp {input.new_repseq} {params.repseq_faa}
#         touch {output.done_flag}
#         """

# # This doesn't handle the case where several samples are run in parallel.
# # We need to track per-sample updates or do it in batch.
# rule update_repseq:
#     input:
#         updated_repseq = rules.cluster_update.output.updated_repseq
#     output:
#         done_flag = Path(config['outdir']) / 'cluster/repseq_update_done.txt'  # ✅ Only track per-sample updates
#     params:
#         repseq_faa = storage.s3("s3://recombia.planter/repseq.faa")  # ✅ Pass the S3 path as a param
#     shell:
#         """
#         aws s3 cp {input.updated_repseq} {params.repseq_faa}
#         touch {output.done_flag}
#         """

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
