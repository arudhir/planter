import os
from pathlib import Path

class Config:
    DUCKDB_PATH = os.environ.get('DUCKDB_PATH') or '/mnt/data4/master.duckdb'

    # Toggle: set to False to search ALL sequences, True for representatives only
    USE_REPRESENTATIVES_ONLY = os.environ.get('USE_REPRESENTATIVES_ONLY', 'false').lower() == 'true'

    # Generate FASTA from database
    if os.environ.get('REPSEQ_FASTA'):
        REPSEQ_FASTA = os.environ.get('REPSEQ_FASTA')
    else:
        from planter.database.utils.duckdb_utils import get_fasta_path
        REPSEQ_FASTA = str(get_fasta_path(DUCKDB_PATH, representatives_only=USE_REPRESENTATIVES_ONLY))

    EXAMPLE_FASTA = os.environ.get('EXAMPLE_FASTA') or '/home/ubuntu/planter/tests/test_enzymes.faa'
    S3_BUCKET = os.environ.get('S3_BUCKET') or 'recombia.planter'
    S3_DB_KEY = os.environ.get('S3_DB_KEY') or 'master.duckdb'
    REPSEQ_OUTPUT_DIR = os.environ.get('REPSEQ_OUTPUT_DIR') or '/mnt/data4'
    DEBUG = False

class DevelopmentConfig(Config):
    DEBUG = True

class ProductionConfig(Config):
    # Production-specific configs
    pass

config = {
    'development': DevelopmentConfig,
    'production': ProductionConfig,
    'default': DevelopmentConfig
}

