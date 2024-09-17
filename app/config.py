import os

class Config:
    MMSEQS_DB = os.environ.get('MMSEQS_DB') or '/Users/aru/Development/planter/tests/seqdb.fa'
    EXAMPLE_FASTA = os.environ.get('EXAMPLE_FASTA') or '/Users/aru/Development/planter/tests/test_enzymes.faa'
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

