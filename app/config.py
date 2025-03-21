import os

class Config:
    REPSEQ_FASTA = os.environ.get('REPSEQ_FASTA') or '/mnt/data4/repseq.faa'
    EXAMPLE_FASTA = os.environ.get('EXAMPLE_FASTA') or '/home/ubuntu/planter/tests/test_enzymes.faa'
    DUCKDB_PATH = os.environ.get('DUCKDB_PATH') or '/mnt/data4/master.duckdb'
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

