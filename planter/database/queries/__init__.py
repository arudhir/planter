import os
from functools import partial

from jinja2 import Template


class BaseQueryManager:
    def __init__(self, con, sql_dir: str):
        self.con = con
        self.sql_dir = sql_dir
        self._load_queries()

    def _load_queries(self):
        """Dynamically load all `.sql` files as callable methods."""
        for file in os.listdir(self.sql_dir):
            if file.endswith(".sql"):
                query_name = file[:-4]  # remove ".sql" extension
                setattr(self, query_name, partial(self._execute_query, query_name))

    def _execute_query(self, query_name: str, params: dict = None, values: list = None):
        sql_path = os.path.join(self.sql_dir, f"{query_name}.sql")
        with open(sql_path, "r") as file:
            query_template = file.read()

        if params:
            if not isinstance(params, dict):
                raise TypeError(
                    f"'params' must be a dictionary, got {type(params)} instead."
                )
            template = Template(query_template)
            query = template.render(**params)
        else:
            query = query_template

        return (
            self.con.execute(query, values).fetchdf()
            if values
            else self.con.execute(query).fetchdf()
        )
