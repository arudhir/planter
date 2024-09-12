import shutil
import hashlib

def create_zip_archive(output_dir):
    logger.info('Creating zip archive')
    shutil.make_archive(
        base_name=f'outputs', 
        root_dir=output_dir,
        format='zip'
    )

def md5sum(filename):
    logger.info('Creating md5sum')
    with open(filename, 'rb') as f:
        md5sum = hashlib.md5(f.read()).hexdigest()
    return md5sum

rule finalize:
    input:
        output_dir = config['outdir']
    output:
        zip_archive = f'{output_dir}.zip'
        md5sum = f'{output_dir}.md5sum'
    shell:
        'create_zip_archive({input.output_dir})'
        'md5sum({input.output_dir})'