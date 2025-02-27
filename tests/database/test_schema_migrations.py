#!/usr/bin/env python3
"""
Tests for database schema migrations.
"""
import os
import tempfile
import unittest
import duckdb
from pathlib import Path

from planter.database.schema.schema import SchemaManager


class TestSchemaMigrations(unittest.TestCase):
    """Test cases for database schema migrations."""

    def setUp(self):
        """Set up test environment."""
        # Create a temporary directory for the database
        self.temp_dir = tempfile.mkdtemp()
        self.db_path = os.path.join(self.temp_dir, "test_schema.duckdb")
        self.conn = duckdb.connect(self.db_path)

    def tearDown(self):
        """Clean up after tests."""
        self.conn.close()
        os.unlink(self.db_path)
        os.rmdir(self.temp_dir)

    def test_migrations_apply_in_order(self):
        """Test that migrations are applied in the correct order."""
        # Initialize schema
        schema_manager = SchemaManager(self.conn)
        
        # Manually apply migrations for testing
        for migration_file in schema_manager._get_migrations():
            migration_sql = schema_manager._read_migration(migration_file.name)
            self.conn.execute(migration_sql)
            self.conn.execute("""
                INSERT INTO schema_version (version, migration_name)
                VALUES (
                    (SELECT COALESCE(MAX(version), 0) + 1 FROM schema_version),
                    ?
                )
            """, [migration_file.name])
        
        # Check that schema_version table was created
        result = self.conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='schema_version'"
        ).fetchall()
        self.assertEqual(len(result), 1)
        
        # Check that migrations were applied in order
        migrations = self.conn.execute(
            "SELECT migration_name FROM schema_version ORDER BY version"
        ).fetchall()
        
        # Get expected migration files from the project
        migrations_dir = Path(__file__).parent.parent.parent / 'planter' / 'database' / 'schema' / 'migrations'
        expected_migrations = sorted([f.name for f in migrations_dir.glob('*.sql')])
        
        # Verify that all migrations were applied in the correct order
        applied_migrations = [m[0] for m in migrations]
        self.assertEqual(applied_migrations, expected_migrations)

    def test_required_tables_exist(self):
        """Test that all required tables are created by migrations."""
        # Initialize schema
        schema_manager = SchemaManager(self.conn)
        
        # Manually apply migrations for testing
        for migration_file in schema_manager._get_migrations():
            migration_sql = schema_manager._read_migration(migration_file.name)
            self.conn.execute(migration_sql)
        
        # Get all tables
        tables = self.conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        ).fetchall()
        table_names = [t[0] for t in tables]
        
        # Check for required tables
        required_tables = [
            'schema_version',
            'sra_metadata',
            'sequences',
            'annotations',
            'go_terms',
            'ec_numbers',
            'kegg_info',
            'clusters',
            'cluster_members',
            'expression'  # New table from our migration
        ]
        
        for table in required_tables:
            self.assertIn(table, table_names)

    def test_expression_table_schema(self):
        """Test that the expression table has the correct schema."""
        # Initialize schema
        schema_manager = SchemaManager(self.conn)
        
        # Manually apply migrations for testing
        for migration_file in schema_manager._get_migrations():
            migration_sql = schema_manager._read_migration(migration_file.name)
            self.conn.execute(migration_sql)
        
        # Get expression table schema using PRAGMA
        columns = self.conn.execute("PRAGMA table_info(expression)").fetchall()
        
        # Extract column names and types
        column_dict = {col[1]: col[2].upper() for col in columns}
        
        # Check required columns
        self.assertIn('seqhash_id', column_dict)
        self.assertIn('sample_id', column_dict)
        self.assertIn('tpm', column_dict)
        self.assertIn('num_reads', column_dict)
        self.assertIn('effective_length', column_dict)
        
        # Check data types (DuckDB uses VARCHAR, DOUBLE)
        self.assertEqual(column_dict['tpm'], 'DOUBLE')
        self.assertEqual(column_dict['num_reads'], 'DOUBLE')
        self.assertEqual(column_dict['effective_length'], 'DOUBLE')

    def test_table_structure(self):
        """Verify that the expression table has the correct structure with foreign keys"""
        # Initialize schema
        schema_manager = SchemaManager(self.conn)
        
        # Manually apply migrations for testing
        for migration_file in schema_manager._get_migrations():
            migration_sql = schema_manager._read_migration(migration_file.name)
            self.conn.execute(migration_sql)
        
        # Check that the 'expression' table exists
        result = self.conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='expression'"
        ).fetchall()
        self.assertEqual(len(result), 1)
        
        # Get table SQL creation statement
        table_sql = self.conn.execute("""
            SELECT sql FROM sqlite_master 
            WHERE type='table' AND name='expression'
        """).fetchone()[0]
        
        # Verify the SQL contains the foreign key constraints
        self.assertIn("FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id)", table_sql)
        self.assertIn("FOREIGN KEY (sample_id) REFERENCES sra_metadata(sample_id)", table_sql)


if __name__ == '__main__':
    unittest.main()