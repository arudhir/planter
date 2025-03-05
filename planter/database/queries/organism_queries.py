from . import BaseQueryManager


class OrganismQueries(BaseQueryManager):
    def get_summary(self):
        """Get summary of sequence counts per organism."""
        return self._execute_query("organism_summary")
