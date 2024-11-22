import os

class Config:
    MMSEQS_DB = os.environ.get('MMSEQS_DB') or '/mnt/data2/planter_outputs/repseq/update_SRR8859648/newRepSeqDB.fasta'
    EXAMPLE_FASTA = os.environ.get('EXAMPLE_FASTA') or '/home/ubuntu/planter/tests/test_enzymes.faa'
    DUCKDB_PATH = os.environ.get('DUCKDB_PATH') or '/mnt/data2/planter_outputs/planter2.duckdb'
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

