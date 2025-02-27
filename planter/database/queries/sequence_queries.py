from . import BaseQueryManager

class SequenceQueries(BaseQueryManager):
    def get_by_id(self, seqhash_id: str):
        """Retrieve complete information for a specific sequence."""
        return self._execute_query("sequence_by_id", values=[seqhash_id])
        
    def get_expression_summary(self, sample_id=None):
        """Get a summary of expression data across the database or for a specific sample."""
        # Always provide values, using NULL for default
        return self._execute_query("expression_summary", values=[sample_id, sample_id, sample_id])
    
    def get_expression_distribution(self, sample_id=None):
        """Get expression level distribution (low, medium, high, very high TPM)."""
        # Always provide the sample_id parameter, using NULL for default
        return self._execute_query("expression_distribution", values=[sample_id])
    
    def get_top_expressed_sequences(self, limit=10, sample_id=None):
        """Get the top expressed sequences by TPM."""
        params = {'limit': limit}
        # Always provide the sample_id parameter, using NULL for default
        return self._execute_query("top_expressed_sequences", params=params, values=[sample_id])
            
    def get_expression_for_sequence(self, seqhash_id):
        """Get expression data for a specific sequence across all samples."""
        # Pass the seqhash_id three times for the WITH clause parameter substitution
        return self._execute_query("sequence_expression", values=[seqhash_id, seqhash_id, seqhash_id])
