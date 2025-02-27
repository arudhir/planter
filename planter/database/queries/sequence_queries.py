from . import BaseQueryManager
from typing import List, Optional, Dict, Any
import os

class SequenceQueries:
    def __init__(self, con, sql_dir: str):
        self._base_query_manager = BaseQueryManager(con, sql_dir)
        self.con = con
        self.sql_dir = sql_dir
    
    def _execute_query(self, query_name, params=None, values=None):
        return self._base_query_manager._execute_query(query_name, params, values)
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
        
    def get_annotation_with_expression(self, sample_id=None, limit=50):
        """Get annotations together with expression levels.
        
        This query joins across annotations, gene-protein mapping, and expression tables
        to provide a comprehensive view of proteins with their expression data.
        
        Args:
            sample_id: Optional filter for a specific sample
            limit: Maximum number of results to return (default: 50)
            
        Returns:
            DataFrame with protein annotations and expression levels
        """
        params = {'limit': limit}
        return self._execute_query("annotation_with_expression", params=params, values=[sample_id])
    
    def search_sequences(self, 
                         sample_ids: Optional[List[str]] = None,
                         min_length: Optional[int] = None,
                         max_length: Optional[int] = None,
                         description: Optional[str] = None,
                         organism: Optional[str] = None,
                         cog_categories: Optional[List[str]] = None,
                         go_terms: Optional[List[str]] = None,
                         limit: int = 100) -> "pd.DataFrame":
        """Search for sequences based on various criteria.
        
        Args:
            sample_ids: List of sample IDs to filter by
            min_length: Minimum sequence length
            max_length: Maximum sequence length
            description: Description text to search for (uses LIKE)
            organism: Organism name to filter by (uses LIKE)
            cog_categories: List of COG categories to filter by
            go_terms: List of GO terms to filter by
            limit: Maximum number of results to return
            
        Returns:
            DataFrame with matching sequences and their annotations
        """
        # Build query parameters
        params: Dict[str, Any] = {'limit': limit}
        
        # Add conditional clauses
        if sample_ids:
            params['sample_id_condition'] = sample_ids
            
        if min_length:
            params['min_length_condition'] = f"s.length >= {min_length}"
            
        if max_length:
            params['max_length_condition'] = f"s.length <= {max_length}"
            
        if description:
            params['description_condition'] = f"a.description LIKE '%{description}%'"
            
        if organism:
            params['organism_condition'] = f"m.organism LIKE '%{organism}%'"
            
        if cog_categories:
            params['cog_category_condition'] = cog_categories
            
        if go_terms:
            params['go_term_condition'] = go_terms
            
        # Read the SQL query directly
        sql_path = os.path.join(self._base_query_manager.sql_dir, "search_sequences.sql")
        with open(sql_path, 'r') as file:
            query_template = file.read()
        
        # Render the template with parameters
        from jinja2 import Template
        template = Template(query_template)
        query = template.render(**params)
        
        # Execute the query directly
        return self._base_query_manager.con.execute(query).fetchdf()
