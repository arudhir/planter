from . import BaseQueryManager

class ClusterQueries(BaseQueryManager):
    def get_cluster_info(self, cluster_id: str):
        """Get detailed info for a specific cluster."""
        return self._execute_query("cluster_info", values=[cluster_id])
    
    def get_cluster_stats(self):
        """Get summary statistics for all clusters."""
        return self._execute_query("cluster_stats")