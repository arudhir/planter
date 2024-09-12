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
    except Exception as e:
        print(f"Failed to upload {zip_filename}: {e}")

rule finalize:
    input:
        eggnog = expand(rules.eggnog.output, sample=config['samples']),
        output_dirs = expand(Path(config['outdir']) / '{sample}', sample=config['samples'])
    output:
        zip_archive = expand(Path(config['outdir']) / '{sample}.zip', sample=config['samples']),
    run:
        for output_dir in input.output_dirs:
            zip_archive = create_zip_archive(output_dir)
            upload_to_s3(zip_archive)
