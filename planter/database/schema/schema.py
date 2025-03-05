import logging
from pathlib import Path
from typing import List


class SchemaManager:
    """Manages database schema initialization and migrations"""

    def __init__(self, connection):
        self.conn = connection
        self.logger = logging.getLogger(__name__)
        self.migrations_dir = Path(__file__).parent / "migrations"

    def _read_migration(self, filename: str) -> str:
        """Read SQL from migration file"""
        with open(self.migrations_dir / filename) as f:
            return f.read()

    def _get_migrations(self) -> List[Path]:
        """Get sorted list of migration files"""
        return sorted(self.migrations_dir.glob("*.sql"))

    def init_database(self):
        """Initialize database with latest schema"""
        try:
            for migration_file in self._get_migrations():
                self.logger.info(f"Applying migration: {migration_file.name}")
                migration_sql = self._read_migration(migration_file.name)
                self.conn.execute(migration_sql)

        except Exception as e:
            self.logger.error(f"Schema initialization failed: {str(e)}")
            raise
