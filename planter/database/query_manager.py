import duckdb
from planter.database.queries.organism_queries import OrganismQueries
from planter.database.queries.sample_queries import SampleQueries
from planter.database.queries.sequence_queries import SequenceQueries

class DatabaseManager:
    """Handles querying the sequence database."""
    def __init__(self, db_path: str):
        self.con = duckdb.connect(db_path)
        self.organism_queries = OrganismQueries(self.con)
        self.sample_queries = SampleQueries(self.con)
        self.sequence_queries = SequenceQueries(self.con)

    def close(self):
        self.con.close()
