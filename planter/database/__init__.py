from .schema.schema import SchemaManager
import logging
import duckdb

class SequenceDBBuilder:
    def __init__(self, db_path: str):
        self.db_path = db_path
        self.con = None
        self.schema_manager = None
        self.logger = logging.getLogger(__name__)

    def __enter__(self):
        self.con = duckdb.connect(self.db_path)
        self.schema_manager = SchemaManager(self.con)
        return self

    def init_database(self):
        """Initialize database with clean schema"""
        self.clean_database()
        self.schema_manager.init_database()
