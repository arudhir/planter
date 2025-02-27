from . import BaseQueryManager

class SampleQueries(BaseQueryManager):
    def get_metadata(self, sample_id: str):
        """Get complete sample information including metadata."""
        return self._execute_query("sample_metadata", values=[sample_id])
