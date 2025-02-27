#!/usr/bin/env python3
"""
Tests for the database builder functionality.
"""
import os
import json
import tempfile
import shutil
from pathlib import Path
import unittest
import duckdb
import pandas as pd
from unittest.mock import patch, MagicMock

from planter.database.builder import SequenceDBBuilder, SamplePaths


class TestSequenceDBBuilder(unittest.TestCase):
    """Test cases for the SequenceDBBuilder class."""

    def setUp(self):
        """Set up test environment."""
        # Create a temporary directory
        self.temp_dir = tempfile.mkdtemp()
        self.db_path = os.path.join(self.temp_dir, "test.duckdb")
        self.output_dir = Path(self.temp_dir) / "output"
        self.output_dir.mkdir(exist_ok=True)
        
        # Create sample directory structure
        self.sample_id = "TEST123"
        self.sample_dir = self.output_dir / self.sample_id
        self.sample_dir.mkdir(exist_ok=True)
        
        # Create transdecoder directory
        self.transdecoder_dir = self.sample_dir / "transdecoder"
        self.transdecoder_dir.mkdir(exist_ok=True)
        
        # Create eggnog directory
        self.eggnog_dir = self.sample_dir / "eggnog"
        self.eggnog_dir.mkdir(exist_ok=True)
        
        # Create quants directory
        self.quants_dir = self.sample_dir / "quants"
        self.quants_dir.mkdir(exist_ok=True)
        
        # Create test files
        self.create_test_files()

    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.temp_dir)

    def create_test_files(self):
        """Create test files for the database builder."""
        # Create a sample pep file
        with open(self.transdecoder_dir / f"{self.sample_id}.pep", "w") as f:
            f.write(">v1_DLS_1234567890abcdef1234567890abcdef1234567890abcdef1234567890abcdef\n")
            f.write("MEPKSLYVGDLHGAFYSIWKLVKSENDNRYLLLVDLNEKIAEMIFNLHGKLDVLSQLPQKK\n")
            f.write(">v1_DLS_0987654321fedcba0987654321fedcba0987654321fedcba0987654321fedcba\n")
            f.write("MARNVLNAIDVLSRLETHLNGLIRLAVDKMDLSEVITSLPKMRRSFSENLNQLTKRVQEL\n")
        
        # Create a sample eggnog annotations file
        with open(self.eggnog_dir / f"{self.sample_id}.emapper.annotations", "w") as f:
            f.write("# eggNOG-mapper - orthology assignment, gene annotation and functional interpretation\n")
            f.write("# emapper version: 2.1.6\n")
            f.write("v1_DLS_1234567890abcdef1234567890abcdef1234567890abcdef1234567890abcdef\t-\t1e-10\t100.0\t-\t-\t-\tHypothetical protein\t-\tGO:0008150,GO:0003674,GO:0005575\tec:1.1.1.1\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n")
            f.write("v1_DLS_0987654321fedcba0987654321fedcba0987654321fedcba0987654321fedcba\t-\t1e-5\t80.0\t-\t-\t-\tAnother protein\t-\tGO:0016491\tec:2.2.2.2\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n")
        
        # Create a sample quant.json file
        quant_data = [
            {
                "Name": "v1_DLS_1234567890abcdef1234567890abcdef1234567890abcdef1234567890abcdef",
                "Length": 60,
                "EffectiveLength": 58.12,
                "TPM": 10.5,
                "NumReads": 100.0,
                "sample": self.sample_id
            },
            {
                "Name": "v1_DLS_0987654321fedcba0987654321fedcba0987654321fedcba0987654321fedcba",
                "Length": 60,
                "EffectiveLength": 58.12,
                "TPM": 5.2,
                "NumReads": 50.0,
                "sample": self.sample_id
            }
        ]
        with open(self.quants_dir / f"{self.sample_id}.quant.json", "w") as f:
            json.dump(quant_data, f)

    @patch('planter.database.builder.get_sra_info')
    def test_database_initialization(self, mock_get_sra_info):
        """Test that the database is properly initialized."""
        # Mock SRA info response
        mock_get_sra_info.return_value = {
            'organism': 'Test Organism',
            'study_title': 'Test Study',
            'study_abstract': 'Abstract',
            'bioproject': 'BIOPROJECT123',
            'biosample': 'BIOSAMPLE123',
            'library_strategy': 'RNA-Seq',
            'library_source': 'TRANSCRIPTOMIC',
            'library_selection': 'cDNA',
            'library_layout': 'PAIRED',
            'instrument': 'ILLUMINA',
            'run_spots': '1000000',
            'run_bases': '100000000',
            'run_published': '2023-01-01'
        }
        
        # Initialize database builder
        with SequenceDBBuilder(self.db_path, self.output_dir) as builder:
            # Build the database
            builder.build_database([self.sample_id])
            
            # Check that the required tables exist
            tables = builder.con.execute("""
                SELECT name FROM sqlite_master 
                WHERE type='table'
            """).fetchall()
            
            table_names = [t[0] for t in tables]
            self.assertIn('sequences', table_names)
            self.assertIn('annotations', table_names)
            self.assertIn('go_terms', table_names)
            self.assertIn('ec_numbers', table_names)
            self.assertIn('expression', table_names)
            self.assertIn('sra_metadata', table_names)
            
            # Verify sequences were loaded
            sequences = builder.con.execute("SELECT * FROM sequences").fetchall()
            self.assertEqual(len(sequences), 2)
            
            # Verify annotations were loaded
            annotations = builder.con.execute("SELECT * FROM annotations").fetchall()
            self.assertEqual(len(annotations), 2)
            
            # Verify GO terms were loaded
            go_terms = builder.con.execute("SELECT * FROM go_terms").fetchall()
            self.assertEqual(len(go_terms), 4)  # 3 for first sequence, 1 for second
            
            # Verify EC numbers were loaded
            ec_numbers = builder.con.execute("SELECT * FROM ec_numbers").fetchall()
            self.assertEqual(len(ec_numbers), 2)  # 1 for each sequence
            
            # Verify expression data was loaded
            expression = builder.con.execute("SELECT * FROM expression").fetchall()
            self.assertEqual(len(expression), 2)  # 1 for each sequence
            
            # Check expression values for first sequence
            expr = builder.con.execute("""
                SELECT tpm, num_reads, effective_length 
                FROM expression 
                WHERE seqhash_id = ?
            """, ["v1_DLS_1234567890abcdef1234567890abcdef1234567890abcdef1234567890abcdef"]).fetchone()
            
            self.assertAlmostEqual(expr[0], 10.5)  # TPM
            self.assertAlmostEqual(expr[1], 100.0)  # NumReads
            self.assertAlmostEqual(expr[2], 58.12)  # EffectiveLength

    @patch('planter.database.builder.get_sra_info')
    def test_database_summary(self, mock_get_sra_info):
        """Test database summary functionality."""
        # Mock SRA info response
        mock_get_sra_info.return_value = {
            'organism': 'Test Organism',
            'study_title': 'Test Study',
        }
        
        # Initialize database builder and build database
        with SequenceDBBuilder(self.db_path, self.output_dir) as builder:
            builder.build_database([self.sample_id])
            
            # Get database summary
            summary = builder.get_database_summary()
            
            # Check summary statistics
            self.assertEqual(summary['total_sequences'].iloc[0], 2)
            self.assertEqual(summary['total_samples'].iloc[0], 1)
            self.assertEqual(summary['annotated_sequences'].iloc[0], 2)
            self.assertEqual(summary['sequences_with_go'].iloc[0], 2)
            self.assertEqual(summary['sequences_with_ec'].iloc[0], 2)
            self.assertEqual(summary['sequences_with_expression'].iloc[0], 2)


if __name__ == '__main__':
    unittest.main()