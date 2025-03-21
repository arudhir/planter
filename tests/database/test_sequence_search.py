#!/usr/bin/env python3
"""
Tests for the sequence search functionality in both the database and Flask app.
"""
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import duckdb
import pandas as pd
from flask import Flask

from planter.database.query_manager import QueryManager
# Import local mock app instead of the actual app package
from tests.database.app_mock import create_app


class TestSequenceSearch(unittest.TestCase):
    """Test cases for sequence search functionality."""

    def setUp(self):
        """Set up a test database and Flask app."""
        # Create a temporary database
        self.temp_dir = tempfile.mkdtemp()
        self.db_path = Path(self.temp_dir) / "test_search.duckdb"

        # Create a connection to the database
        self.con = duckdb.connect(str(self.db_path))

        # Set up the database schema
        self.sample_id = "SRR12068547"
        self._create_test_database()

        # Create a test Flask app
        self.app = create_app("testing")
        self.app.config["TESTING"] = True
        self.app.config["DUCKDB_PATH"] = str(self.db_path)

        # Use the existing test_enzymes.faa file
        test_faa_path = str(Path(__file__).parent.parent / "test_enzymes.faa")
        if os.path.exists(test_faa_path):
            self.app.config["EXAMPLE_FASTA"] = test_faa_path
        else:
            # If the file doesn't exist, create a temporary one
            temp_faa = Path(self.temp_dir) / "test_enzymes.faa"
            with open(temp_faa, "w") as f:
                f.write(">test_enzyme\nACTG\n")
            self.app.config["EXAMPLE_FASTA"] = str(temp_faa)

        # Create a test client
        self.client = self.app.test_client()

    def tearDown(self):
        """Clean up after tests."""
        self.con.close()
        os.remove(self.db_path)
        os.rmdir(self.temp_dir)

    def _create_test_database(self):
        """Create a test database with schema and sample data."""
        # Basic schema for testing
        self.con.execute(
            "CREATE TABLE sra_metadata (sample_id VARCHAR PRIMARY KEY, organism VARCHAR);"
        )
        self.con.execute(
            "INSERT INTO sra_metadata VALUES ('SRR12068547', 'Test Organism');"
        )

        self.con.execute(
            """
        CREATE TABLE sequences (
            seqhash_id VARCHAR PRIMARY KEY,
            sequence VARCHAR NOT NULL,
            sample_id VARCHAR NOT NULL,
            length INTEGER,
            is_representative BOOLEAN DEFAULT FALSE
        );
        """
        )

        self.con.execute(
            """
        INSERT INTO sequences VALUES
        ('test_seqhash_1', 'ACTG', 'SRR12068547', 4, TRUE),
        ('test_seqhash_2', 'GGCC', 'SRR12068547', 4, FALSE);
        """
        )

        self.con.execute(
            """
        CREATE TABLE annotations (
            seqhash_id VARCHAR,
            description VARCHAR,
            preferred_name VARCHAR,
            sample_id VARCHAR,
            PRIMARY KEY (seqhash_id, sample_id)
        );
        """
        )

        self.con.execute(
            """
        INSERT INTO annotations VALUES
        ('test_seqhash_1', 'Test enzyme 1', 'enzyme1', 'SRR12068547'),
        ('test_seqhash_2', 'Test enzyme 2', 'enzyme2', 'SRR12068547');
        """
        )

    def test_load_example_error_handling(self):
        """Test that the load_example endpoint properly handles errors."""
        # First, test with a valid example file path
        response = self.client.get("/load_example")
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIn("sequence", data)

        # Now test with an invalid file path
        self.app.config["EXAMPLE_FASTA"] = "/path/that/does/not/exist.faa"
        response = self.client.get("/load_example")
        self.assertEqual(response.status_code, 500)

        # Verify we get a proper error message
        data = json.loads(response.data)
        self.assertIn("error", data)
        self.assertEqual(data["error"], "Failed to load example sequence")

        # Reset to valid path for other tests
        test_faa_path = str(Path(__file__).parent.parent / "test_enzymes.faa")
        if os.path.exists(test_faa_path):
            self.app.config["EXAMPLE_FASTA"] = test_faa_path
        else:
            temp_faa = Path(self.temp_dir) / "test_enzymes.faa"
            self.app.config["EXAMPLE_FASTA"] = str(temp_faa)


if __name__ == "__main__":
    unittest.main()
