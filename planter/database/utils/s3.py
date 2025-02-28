#!/usr/bin/env python3
"""
S3 utility functions for the planter workflow.
"""
import os
from pathlib import Path
import boto3
import botocore
import logging
from typing import Union

# Disable insecure request warnings
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
boto3.set_stream_logger(name='botocore', level=logging.INFO)


def create_zip_archive(output_dir: Union[str, Path]) -> str:
    """
    Create a zip archive of the given directory.
    
    Args:
        output_dir: Directory to archive
        
    Returns:
        Path to the created zip archive
    """
    import shutil
    zip_archive = shutil.make_archive(
        base_name=f'{output_dir}', 
        root_dir=output_dir,
        format='zip'
    )
    return zip_archive


def upload_to_s3(output_dir: Union[str, Path], sample: str, bucket: str) -> bool:
    """
    Upload all files from output_dir to `bucket` under an S3 prefix = sample.
    If a file already exists in S3, skip it (do not fail).
    
    Args:
        output_dir: Directory containing files to upload
        sample: Sample ID to use as S3 prefix
        bucket: S3 bucket name
        
    Returns:
        True if all files were either skipped or uploaded successfully;
        False if an unrecoverable error occurs.
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