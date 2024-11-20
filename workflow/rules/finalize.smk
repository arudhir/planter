from pathlib import Path
import shutil
import hashlib
import boto3
import logging
boto3.set_stream_logger(name='botocore', level=logging.INFO)

def create_zip_archive(output_dir):
    zip_archive = shutil.make_archive(
        base_name=f'{output_dir}', 
        root_dir=output_dir,
        format='zip'
    )
    return zip_archive


def upload_to_s3(zip_archive):
    s3 = boto3.client('s3')
    zip_filename = os.path.basename(zip_archive)  # only the filename, no local path
    print(f"Uploading {zip_filename} to S3 bucket: {config['s3_bucket']}")
    try:
        s3.upload_file(zip_archive, config['s3_bucket'], zip_filename)  # using just the filename in S3
        print(f"Successfully uploaded {zip_filename} to {config['s3_bucket']}")
        success = True
    except Exception as e:
        print(f"Failed to upload {zip_filename}: {e}")
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
            'scripts/get_qc_stats.py '
            '--sample {wildcards.sample} '
            '--fastp {input.fastp} '
            '--salmon_metadata {input.salmon_metadata} '
            '--transcripts {input.transcripts} '
            '--quantsf {input.quantsf} '
            '--eggnog {input.eggnog} '
            '--output_file {output.qc_stats}'
        )


rule finalize:
    input:
        analyze_eggnog = expand(rules.analyze_eggnog.output, sample=config['samples']),
        quant = expand(rules.quant.output, sample=config['samples']),
    output:
        zip_archive = expand(Path(config['outdir']) / '{sample}.zip', sample=config['samples']),
        done = temp(expand(Path(config['outdir']) / '{sample}/{sample}_s3_upload.done', sample=config['samples']))
    run:
        samples = config['samples']
        if type(config['samples']) != list:
            samples = [config['samples']]
        for sample in samples:
            print('Finalizing sample: ', sample)
            output_dir = Path(config['outdir']) / sample
            zip_archive = create_zip_archive(output_dir)
            success = upload_to_s3(zip_archive)
            if success:
                (Path(config['outdir']) / f'{sample}/{sample}_s3_upload.done').touch()
