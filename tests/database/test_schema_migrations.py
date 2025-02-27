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
        schema_manager.init_database()
        
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
        schema_manager.init_database()
        
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
        schema_manager.init_database()
        
        # Check expression table columns
        columns = self.conn.execute("""
            SELECT column_name, data_type 
            FROM information_schema.columns 
            WHERE table_name='expression'
        """).fetchall()
        
        # Convert to dict for easier checking
        column_dict = {col[0]: col[1] for col in columns}
        
        # Check required columns
        self.assertIn('seqhash_id', column_dict)
        self.assertIn('sample_id', column_dict)
        self.assertIn('tpm', column_dict)
        self.assertIn('num_reads', column_dict)
        self.assertIn('effective_length', column_dict)
        
        # Check data types
        self.assertEqual(column_dict['tpm'].upper(), 'DOUBLE')
        self.assertEqual(column_dict['num_reads'].upper(), 'DOUBLE')
        self.assertEqual(column_dict['effective_length'].upper(), 'DOUBLE')

    def test_foreign_key_constraints(self):
        """Test that foreign key constraints are properly defined."""
        # Initialize schema
        schema_manager = SchemaManager(self.conn)
        schema_manager.init_database()
        
        # Check foreign keys for expression table
        fk_query = """
        SELECT 
            fk.constraint_table as child_table,
            fk.column_name as child_column,
            fk.foreign_table as parent_table,
            fk.foreign_column as parent_column
        FROM information_schema.foreign_keys fk
        WHERE fk.constraint_table = 'expression'
        """
        
        foreign_keys = self.conn.execute(fk_query).fetchall()
        
        # Should have two foreign keys: seqhash_id -> sequences.seqhash_id and sample_id -> sra_metadata.sample_id
        self.assertEqual(len(foreign_keys), 2)
        
        # Collect FKs in an easier format to check
        fk_dict = {}
        for fk in foreign_keys:
            child_col = fk[1]
            parent_table = fk[2]
            parent_col = fk[3]
            fk_dict[child_col] = (parent_table, parent_col)
        
        # Check specific FKs
        self.assertIn('seqhash_id', fk_dict)
        self.assertEqual(fk_dict['seqhash_id'][0], 'sequences')
        self.assertEqual(fk_dict['seqhash_id'][1], 'seqhash_id')
        
        self.assertIn('sample_id', fk_dict)
        self.assertEqual(fk_dict['sample_id'][0], 'sra_metadata')
        self.assertEqual(fk_dict['sample_id'][1], 'sample_id')


if __name__ == '__main__':
    unittest.main()