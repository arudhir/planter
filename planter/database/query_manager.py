import os
import duckdb
from planter.database.queries import BaseQueryManager
from planter.database.queries.sequence_queries import SequenceQueries
from planter.database.queries.cluster_queries import ClusterQueries
from planter.database.queries.organism_queries import OrganismQueries
from planter.database.queries.sample_queries import SampleQueries

class QueryManager:
    """Manages access to various query components."""
    def __init__(self, con, sql_dir: str):
        self.con = con
        self.sql_dir = sql_dir
        
        # Initialize query components
        self.sequences = SequenceQueries(con, sql_dir)
        self.clusters = ClusterQueries(con, sql_dir)
        self.organisms = OrganismQueries(con, sql_dir)
        self.samples = SampleQueries(con, sql_dir)

class DatabaseManager:
    """Manages database connection and queries."""
    def __init__(self, db_path: str):
        sql_dir = os.path.join(os.path.dirname(__file__), 'queries', 'sql')  # directory for SQL files
        self.con = duckdb.connect(db_path)
        self.queries = QueryManager(self.con, sql_dir)

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()

    def close(self):
        """Close the database connection."""
        self.con.close()
