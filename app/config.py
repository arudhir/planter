import os

class Config:
    REPSEQ_FASTA = os.environ.get('REPSEQ_FASTA') or '/mnt/data3/repseq.faa'
    EXAMPLE_FASTA = os.environ.get('EXAMPLE_FASTA') or '/home/ubuntu/planter/tests/test_enzymes.faa'
    DUCKDB_PATH = os.environ.get('DUCKDB_PATH') or '/mnt/data3/master.duckdb'
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

