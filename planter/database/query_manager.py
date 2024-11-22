import os
import duckdb
from planter.database.queries import BaseQueryManager

class DatabaseManager:
    """Manages database connection and queries."""
    def __init__(self, db_path: str):
        sql_dir = os.path.join(os.path.dirname(__file__), 'queries', 'sql')  # directory for SQL files
        self.con = duckdb.connect(db_path)
        self.query_manager = BaseQueryManager(self.con, sql_dir)

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()

    def close(self):
        """Close the database connection."""
        self.con.close()
