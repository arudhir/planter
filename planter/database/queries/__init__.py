import os

class BaseQueryManager:
    def __init__(self, con):
        self.con = con

    def _load_query(self, query_name: str) -> str:
        """Load SQL query from file."""
        path = os.path.join(os.path.dirname(__file__), 'sql', f"{query_name}.sql")
        with open(path, 'r') as f:
            return f.read()

    def execute_query(self, query_name: str, params: list = []):
        """Execute a query by name with parameters."""
        query = self._load_query(query_name)
        return self.con.execute(query, params).fetchdf()
