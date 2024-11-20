from . import BaseQueryManager

class SequenceQueries(BaseQueryManager):
    def get_by_id(self, seqhash_id: str):
        """Retrieve complete information for a specific sequence."""
        return self.execute_query("sequence_by_id", [seqhash_id])
