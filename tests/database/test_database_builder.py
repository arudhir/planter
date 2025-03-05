#!/usr/bin/env python3
"""
Tests for the database builder functionality.
"""
import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import duckdb
import pandas as pd

from planter.database.builder import SamplePaths, SequenceDBBuilder


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
            f.write(
                ">v1_DLS_1234567890abcdef1234567890abcdef1234567890abcdef1234567890abcdef.p1\n"
            )
            f.write("MEPKSLYVGDLHGAFYSIWKLVKSENDNRYLLLVDLNEKIAEMIFNLHGKLDVLSQLPQKK\n")
            f.write(
                ">v1_DLS_0987654321fedcba0987654321fedcba0987654321fedcba0987654321fedcba.p1\n"
            )
            f.write("MARNVLNAIDVLSRLETHLNGLIRLAVDKMDLSEVITSLPKMRRSFSENLNQLTKRVQEL\n")

        # Create a sample eggnog annotations file
        with open(self.eggnog_dir / f"{self.sample_id}.emapper.annotations", "w") as f:
            f.write(
                "# eggNOG-mapper - orthology assignment, gene annotation and functional interpretation\n"
            )
            f.write("# emapper version: 2.1.6\n")
            f.write(
                "v1_DLS_1234567890abcdef1234567890abcdef1234567890abcdef1234567890abcdef.p1\t-\t1e-10\t100.0\t-\t-\t-\tHypothetical protein\t-\tGO:0008150,GO:0003674,GO:0005575\tec:1.1.1.1\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n"
            )
            f.write(
                "v1_DLS_0987654321fedcba0987654321fedcba0987654321fedcba0987654321fedcba.p1\t-\t1e-5\t80.0\t-\t-\t-\tAnother protein\t-\tGO:0016491\tec:2.2.2.2\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n"
            )

        # Create a sample quant.json file
        # Note: for expression data, we use the gene ID (without .p1)
        quant_data = [
            {
                "Name": "v1_DLS_1234567890abcdef1234567890abcdef1234567890abcdef1234567890abcdef",
                "Length": 60,
                "EffectiveLength": 58.12,
                "TPM": 10.5,
                "NumReads": 100.0,
                "sample": self.sample_id,
            },
            {
                "Name": "v1_DLS_0987654321fedcba0987654321fedcba0987654321fedcba0987654321fedcba",
                "Length": 60,
                "EffectiveLength": 58.12,
                "TPM": 5.2,
                "NumReads": 50.0,
                "sample": self.sample_id,
            },
        ]
        with open(self.quants_dir / f"{self.sample_id}.quant.json", "w") as f:
            json.dump(quant_data, f)

    @patch("planter.database.builder.get_sra_info")
    def test_database_initialization(self, mock_get_sra_info):
        """Test that the database is properly initialized."""
        # Mock SRA info response
        mock_get_sra_info.return_value = {
            "organism": "Test Organism",
            "study_title": "Test Study",
            "study_abstract": "Abstract",
            "bioproject": "BIOPROJECT123",
            "biosample": "BIOSAMPLE123",
            "library_strategy": "RNA-Seq",
            "library_source": "TRANSCRIPTOMIC",
            "library_selection": "cDNA",
            "library_layout": "PAIRED",
            "instrument": "ILLUMINA",
            "run_spots": "1000000",
            "run_bases": "100000000",
            "run_published": "2023-01-01",
        }

        # Initialize database builder
        with SequenceDBBuilder(self.db_path, self.output_dir) as builder:
            # Build the database
            builder.build_database([self.sample_id])

            # Check that the required tables exist
            tables = builder.con.execute(
                """
                SELECT name FROM sqlite_master 
                WHERE type='table'
            """
            ).fetchall()

            table_names = [t[0] for t in tables]
            self.assertIn("sequences", table_names)
            self.assertIn("annotations", table_names)
            self.assertIn("go_terms", table_names)
            self.assertIn("ec_numbers", table_names)
            self.assertIn("expression", table_names)
            self.assertIn("sra_metadata", table_names)

            # Verify sequences were loaded
            sequences = builder.con.execute("SELECT * FROM sequences").fetchall()
            self.assertEqual(len(sequences), 2)

            # Verify the repseq_id field is populated correctly
            for seq in sequences:
                self.assertIsNotNone(seq[5], "repseq_id should not be None")

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
            expr = builder.con.execute(
                """
                SELECT tpm, num_reads, effective_length 
                FROM expression 
                WHERE gene_seqhash_id = ?
            """,
                [
                    "v1_DLS_1234567890abcdef1234567890abcdef1234567890abcdef1234567890abcdef"
                ],
            ).fetchone()

            self.assertAlmostEqual(expr[0], 10.5)  # TPM
            self.assertAlmostEqual(expr[1], 100.0)  # NumReads
            self.assertAlmostEqual(expr[2], 58.12)  # EffectiveLength

            # Verify gene_protein_map table was properly populated
            gene_protein_maps = builder.con.execute(
                "SELECT * FROM gene_protein_map"
            ).fetchall()
            self.assertEqual(
                len(gene_protein_maps), 2, "Should have two entries in gene_protein_map"
            )

            # Verify mapping for the first sequence
            gene_protein_map = builder.con.execute(
                """
                SELECT gene_seqhash_id, protein_seqhash_id 
                FROM gene_protein_map 
                WHERE gene_seqhash_id = ?
            """,
                [
                    "v1_DLS_1234567890abcdef1234567890abcdef1234567890abcdef1234567890abcdef"
                ],
            ).fetchone()

            self.assertEqual(
                gene_protein_map[0],
                "v1_DLS_1234567890abcdef1234567890abcdef1234567890abcdef1234567890abcdef",
            )
            self.assertEqual(
                gene_protein_map[1],
                "v1_DLS_1234567890abcdef1234567890abcdef1234567890abcdef1234567890abcdef.p1",
            )

    @patch("planter.database.builder.get_sra_info")
    def test_database_summary(self, mock_get_sra_info):
        """Test database summary functionality."""
        # Mock SRA info response
        mock_get_sra_info.return_value = {
            "organism": "Test Organism",
            "study_title": "Test Study",
        }

        # Initialize database builder and build database
        with SequenceDBBuilder(self.db_path, self.output_dir) as builder:
            builder.build_database([self.sample_id])

            # Get database summary
            summary = builder.get_database_summary()

            # Check summary statistics
            self.assertEqual(summary["total_sequences"].iloc[0], 2)
            self.assertEqual(summary["total_samples"].iloc[0], 1)
            self.assertEqual(summary["annotated_sequences"].iloc[0], 2)
            self.assertEqual(summary["sequences_with_go"].iloc[0], 2)
            self.assertEqual(summary["sequences_with_ec"].iloc[0], 2)
            self.assertEqual(summary["sequences_with_expression"].iloc[0], 2)

    def test_database_update_workflow(self):
        """Test the sequence clustering part of the database update workflow."""
        # Create a test database path
        db_path = os.path.join(self.temp_dir, "test_cluster.duckdb")

        # Initialize the database with tables and test data
        with duckdb.connect(db_path) as con:
            # Create the sequences table
            con.execute(
                """
                CREATE TABLE sequences (
                    seqhash_id VARCHAR PRIMARY KEY,
                    sequence VARCHAR,
                    sample_id VARCHAR,
                    assembly_date TIMESTAMP,
                    is_representative BOOLEAN,
                    repseq_id VARCHAR,
                    length INTEGER
                )
            """
            )

            # Create the clusters table
            con.execute(
                """
                CREATE TABLE clusters (
                    cluster_id VARCHAR PRIMARY KEY,
                    representative_seqhash_id VARCHAR,
                    size INTEGER
                )
            """
            )

            # Create the cluster_members table
            con.execute(
                """
                CREATE TABLE cluster_members (
                    seqhash_id VARCHAR PRIMARY KEY,
                    cluster_id VARCHAR
                )
            """
            )

            # Insert test sequences
            con.execute(
                """
                INSERT INTO sequences VALUES
                ('seq1', 'ACGT', 'sample1', CURRENT_TIMESTAMP, FALSE, 'seq1', 4),
                ('seq2', 'TGCA', 'sample1', CURRENT_TIMESTAMP, FALSE, 'seq2', 4),
                ('seq3', 'GCAT', 'sample2', CURRENT_TIMESTAMP, FALSE, 'seq3', 4)
            """
            )

        # Create a cluster TSV file (simulating MMSeqs2 clustering)
        cluster_tsv_path = os.path.join(self.temp_dir, "clusters.tsv")
        with open(cluster_tsv_path, "w") as f:
            # Make seq1 the representative for all sequences
            f.write("seq1\tseq1\n")
            f.write("seq1\tseq2\n")
            f.write("seq1\tseq3\n")

        # Update the database with clustering information
        from planter.database.utils.duckdb_utils import \
            update_duckdb_with_cluster_info

        update_duckdb_with_cluster_info(db_path, cluster_tsv_path)

        # Verify the results by connecting to the updated database
        with duckdb.connect(db_path) as con:
            # Check that all sequences have seq1 as their representative
            repseqs = con.execute(
                "SELECT seqhash_id, repseq_id FROM sequences"
            ).fetchall()
            self.assertEqual(len(repseqs), 3, "Should have 3 sequences")

            for seqhash_id, repseq_id in repseqs:
                self.assertEqual(
                    repseq_id,
                    "seq1",
                    f"Expected seq1 as repseq_id for {seqhash_id}, got {repseq_id}",
                )

            # Check cluster entry
            cluster = con.execute(
                "SELECT cluster_id, representative_seqhash_id, size FROM clusters"
            ).fetchone()
            self.assertIsNotNone(cluster, "Cluster should exist")
            self.assertEqual(cluster[1], "seq1", "Representative should be seq1")
            self.assertEqual(cluster[2], 3, "Cluster size should be 3")

            # Check cluster memberships
            members = con.execute(
                """
                SELECT cm.seqhash_id, c.representative_seqhash_id
                FROM cluster_members cm
                JOIN clusters c ON cm.cluster_id = c.cluster_id
            """
            ).fetchall()
            self.assertEqual(len(members), 3, "Should have 3 cluster members")

            # Check that all members point to seq1 as their representative
            for seqhash_id, representative_id in members:
                self.assertEqual(
                    representative_id,
                    "seq1",
                    f"Expected representative 'seq1' for {seqhash_id}",
                )

    @patch("planter.database.builder.get_sra_info")
    def test_database_expression_loading(self, mock_get_sra_info):
        """Test that expression data is properly loaded from quant.json files."""
        # Mock SRA info response
        mock_get_sra_info.return_value = {
            "organism": "Test Organism",
            "study_title": "Test Study",
        }

        # Create a database and process a sample
        with SequenceDBBuilder(self.db_path, self.output_dir) as builder:
            builder.init_database()
            builder.process_sample(self.sample_id)

            # Query the database for expression data
            expr_data = builder.con.execute(
                """
                SELECT e.gene_seqhash_id, e.tpm, e.num_reads, e.effective_length, 
                       g.protein_seqhash_id
                FROM expression e
                JOIN gene_protein_map g ON e.gene_seqhash_id = g.gene_seqhash_id
                WHERE e.sample_id = ?
            """,
                [self.sample_id],
            ).fetchall()

            # Verify the expression data was correctly loaded
            self.assertEqual(
                len(expr_data), 2, "Should have expression data for 2 genes"
            )

            # Check the first gene's expression data
            gene1 = "v1_DLS_1234567890abcdef1234567890abcdef1234567890abcdef1234567890abcdef"
            gene1_expr = builder.con.execute(
                """
                SELECT tpm, num_reads, effective_length FROM expression
                WHERE gene_seqhash_id = ? AND sample_id = ?
            """,
                [gene1, self.sample_id],
            ).fetchone()

            self.assertIsNotNone(gene1_expr, "Should have expression data for gene1")
            self.assertAlmostEqual(gene1_expr[0], 10.5, "TPM should match")
            self.assertAlmostEqual(gene1_expr[1], 100.0, "NumReads should match")
            self.assertAlmostEqual(gene1_expr[2], 58.12, "EffectiveLength should match")

            # Verify the gene-protein mapping
            gene_protein = builder.con.execute(
                """
                SELECT protein_seqhash_id FROM gene_protein_map
                WHERE gene_seqhash_id = ?
            """,
                [gene1],
            ).fetchone()

            self.assertIsNotNone(gene_protein, "Should have gene-protein mapping")
            self.assertEqual(
                gene_protein[0],
                f"{gene1}.p1",
                "Protein ID should match gene ID with .p1 suffix",
            )


if __name__ == "__main__":
    unittest.main()
