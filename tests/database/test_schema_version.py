#!/usr/bin/env python3
"""
Tests for schema versioning and compatibility.
"""
import os
import shutil
import tempfile
import unittest
from pathlib import Path

import duckdb

from planter.database.schema.schema_version import (SCHEMA_VERSIONS,
                                                    ensure_compatibility,
                                                    get_db_schema_version,
                                                    upgrade_schema)


class TestSchemaVersion(unittest.TestCase):
    """Test cases for the schema versioning system."""

    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.db_v1_path = os.path.join(self.temp_dir, "db_v1.duckdb")
        self.db_v2_path = os.path.join(self.temp_dir, "db_v2.duckdb")

        # Create a v1 schema database
        self._create_v1_database(self.db_v1_path)

        # Create a v2 schema database
        self._create_v2_database(self.db_v2_path)

    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.temp_dir)

    def _create_v1_database(self, db_path):
        """Create a database with v1 schema."""
        con = duckdb.connect(db_path)

        # Create v1 schema
        con.execute(
            """
            CREATE TABLE sequences (
                seqhash_id VARCHAR PRIMARY KEY,
                sequence VARCHAR,
                sample_id VARCHAR,
                assembly_date TIMESTAMP,
                repseq_id VARCHAR,
                length INTEGER
            )
        """
        )

        con.execute(
            """
            CREATE TABLE clusters (
                cluster_id VARCHAR PRIMARY KEY,
                representative_seqhash_id VARCHAR,
                size INTEGER
            )
        """
        )

        con.execute(
            """
            CREATE TABLE cluster_members (
                seqhash_id VARCHAR PRIMARY KEY,
                cluster_id VARCHAR
            )
        """
        )

        # Insert some test data
        con.execute(
            """
            INSERT INTO sequences VALUES
            ('seq1', 'ACTG', 'sample1', NULL, 'seq1', 4),
            ('seq2', 'TGCA', 'sample1', NULL, 'seq1', 4)
        """
        )

        con.close()

    def _create_v2_database(self, db_path):
        """Create a database with v2 schema."""
        con = duckdb.connect(db_path)

        # Create v2 schema
        con.execute(
            """
            CREATE TABLE sequences (
                seqhash_id VARCHAR PRIMARY KEY,
                sequence VARCHAR,
                sample_id VARCHAR,
                assembly_date TIMESTAMP,
                is_representative BOOLEAN DEFAULT FALSE,
                repseq_id VARCHAR,
                length INTEGER
            )
        """
        )

        con.execute(
            """
            CREATE TABLE clusters (
                cluster_id VARCHAR PRIMARY KEY,
                representative_seqhash_id VARCHAR,
                size INTEGER
            )
        """
        )

        con.execute(
            """
            CREATE TABLE cluster_members (
                seqhash_id VARCHAR PRIMARY KEY,
                cluster_id VARCHAR
            )
        """
        )

        # Insert some test data
        con.execute(
            """
            INSERT INTO sequences VALUES
            ('seq1', 'ACTG', 'sample1', NULL, TRUE, 'seq1', 4),
            ('seq2', 'TGCA', 'sample1', NULL, FALSE, 'seq1', 4)
        """
        )

        # Create schema_version table
        con.execute(
            """
            CREATE TABLE schema_version (
                version INTEGER PRIMARY KEY,
                migration_name VARCHAR NOT NULL,
                applied_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """
        )

        con.execute(
            """
            INSERT INTO schema_version (version, migration_name)
            VALUES (2, 'upgrade_to_v2')
        """
        )

        con.close()

    def test_get_db_schema_version(self):
        """Test detecting database schema version."""
        # Test v1 schema detection
        v1_version = get_db_schema_version(self.db_v1_path)
        self.assertEqual(v1_version, 1, "Should detect v1 schema")

        # Test v2 schema detection
        v2_version = get_db_schema_version(self.db_v2_path)
        self.assertEqual(v2_version, 2, "Should detect v2 schema")

    def test_upgrade_schema(self):
        """Test upgrading database schema."""
        # Create a copy of the v1 database to upgrade
        upgrade_db_path = os.path.join(self.temp_dir, "upgrade_db.duckdb")
        shutil.copy(self.db_v1_path, upgrade_db_path)

        # Upgrade to v2
        new_version = upgrade_schema(upgrade_db_path, target_version=2)
        self.assertEqual(new_version, 2, "Schema should be upgraded to v2")

        # Verify the upgrade
        con = duckdb.connect(upgrade_db_path)

        # Check is_representative column was added
        has_rep_column = con.execute(
            """
            SELECT COUNT(*) FROM pragma_table_info('sequences') 
            WHERE name = 'is_representative'
        """
        ).fetchone()[0]
        self.assertTrue(has_rep_column, "is_representative column should be added")

        # Check schema_version table was created
        has_version_table = con.execute(
            """
            SELECT COUNT(*) FROM sqlite_master 
            WHERE type='table' AND name='schema_version'
        """
        ).fetchone()[0]
        self.assertTrue(has_version_table, "schema_version table should be created")

        # Check version in schema_version table
        version = con.execute("SELECT MAX(version) FROM schema_version").fetchone()[0]
        self.assertEqual(version, 2, "schema_version should be 2")

        con.close()

    def test_ensure_compatibility(self):
        """Test ensuring schema compatibility."""
        # Create a copy of the v1 database
        compat_db_path = os.path.join(self.temp_dir, "compat_db.duckdb")
        shutil.copy(self.db_v1_path, compat_db_path)

        # Test ensuring compatibility to v2
        version, was_upgraded = ensure_compatibility(compat_db_path, required_version=2)
        self.assertEqual(version, 2, "Should report version 2")
        self.assertTrue(was_upgraded, "Should report that upgrade was performed")

        # Verify the database was actually upgraded
        con = duckdb.connect(compat_db_path)
        has_rep_column = con.execute(
            """
            SELECT COUNT(*) FROM pragma_table_info('sequences') 
            WHERE name = 'is_representative'
        """
        ).fetchone()[0]
        self.assertTrue(has_rep_column, "Database should be upgraded to v2 schema")
        con.close()

        # Test ensuring compatibility when already at required version
        version, was_upgraded = ensure_compatibility(compat_db_path, required_version=2)
        self.assertEqual(version, 2, "Should report version 2")
        self.assertFalse(was_upgraded, "Should report that no upgrade was needed")

        # Test ensuring compatibility to latest version
        version, was_upgraded = ensure_compatibility(compat_db_path)
        self.assertEqual(
            version, max(SCHEMA_VERSIONS.keys()), "Should report latest version"
        )
        self.assertFalse(was_upgraded, "Should report that no upgrade was needed")

    def test_compatibility_with_newer_version(self):
        """Test backward compatibility with newer schema version."""
        # Create a database with future schema v3 (simulated)
        future_db_path = os.path.join(self.temp_dir, "future_db.duckdb")
        con = duckdb.connect(future_db_path)

        # Create v2 schema with an extra column to simulate v3
        con.execute(
            """
            CREATE TABLE sequences (
                seqhash_id VARCHAR PRIMARY KEY,
                sequence VARCHAR,
                sample_id VARCHAR,
                assembly_date TIMESTAMP,
                is_representative BOOLEAN DEFAULT FALSE,
                repseq_id VARCHAR,
                length INTEGER,
                future_column VARCHAR
            )
        """
        )

        # Insert test data
        con.execute(
            """
            INSERT INTO sequences VALUES
            ('seq1', 'ACTG', 'sample1', NULL, TRUE, 'seq1', 4, 'future_data')
        """
        )

        # Create schema_version table with future version
        con.execute(
            """
            CREATE TABLE schema_version (
                version INTEGER PRIMARY KEY,
                migration_name VARCHAR NOT NULL,
                applied_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """
        )

        con.execute(
            """
            INSERT INTO schema_version (version, migration_name)
            VALUES (3, 'upgrade_to_v3')
        """
        )

        con.close()

        # Test that our code correctly identifies this as a future version
        version = get_db_schema_version(future_db_path)
        self.assertEqual(version, 3, "Should detect future schema v3")

        # Test that code doesn't try to downgrade the schema
        version, was_upgraded = ensure_compatibility(future_db_path, required_version=2)
        self.assertEqual(version, 3, "Should report version 3")
        self.assertFalse(was_upgraded, "Should not downgrade schema")


if __name__ == "__main__":
    unittest.main()
