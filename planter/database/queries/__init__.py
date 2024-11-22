import os
from functools import partial

class BaseQueryManager:
    def __init__(self, con, sql_dir: str):
        self.con = con
        self.sql_dir = sql_dir
        self._load_queries()

    def _load_queries(self):
        """Dynamically load all `.sql` files as callable methods."""
        for file in os.listdir(self.sql_dir):
            if file.endswith('.sql'):
                query_name = file[:-4]  # remove ".sql" extension
                setattr(self, query_name, partial(self._execute_query, query_name))

    def _execute_query(self, query_name: str, *params):
        """Load and execute a query by name with optional parameters."""
        sql_path = os.path.join(self.sql_dir, f"{query_name}.sql")
        with open(sql_path, 'r') as file:
            query = file.read()
        return self.con.execute(query, params).fetchdf()
