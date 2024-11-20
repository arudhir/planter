from . import BaseQueryManager

class OrganismQueries(BaseQueryManager):
    def get_summary(self):
        """Get summary of sequence counts per organism."""
        return self.execute_query("organism_summary")
