#!/usr/bin/env python3
"""
Tests for the ClusteringManager class.
"""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch, MagicMock

import duckdb
import pandas as pd

from planter.clustering.manager import ClusteringManager
from planter.clustering.cluster import ClusteringParams, ClusteringResult


class TestClusteringManagerSetup(unittest.TestCase):
    """Test ClusteringManager initialization and schema setup."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.db_path = Path(self.temp_dir) / "test.duckdb"

        # Create a basic database with sequences table
        con = duckdb.connect(str(self.db_path))
        con.execute("""
            CREATE TABLE sra_metadata (
                sample_id VARCHAR PRIMARY KEY,
                organism VARCHAR
            )
        """)
        con.execute("""
            CREATE TABLE sequences (
                seqhash_id VARCHAR PRIMARY KEY,
                sequence VARCHAR NOT NULL,
                sample_id VARCHAR NOT NULL,
                length INTEGER NOT NULL
            )
        """)
        con.execute("INSERT INTO sra_metadata VALUES ('sample1', 'Test organism')")
        con.execute("INSERT INTO sequences VALUES ('seq1', 'ACGT', 'sample1', 4)")
        con.execute("INSERT INTO sequences VALUES ('seq2', 'TGCA', 'sample1', 4)")
        con.execute("INSERT INTO sequences VALUES ('seq3', 'AAAA', 'sample1', 4)")
        con.close()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_initialization(self):
        """Test manager initialization."""
        manager = ClusteringManager(str(self.db_path))
        self.assertEqual(manager.db_path, str(self.db_path))

    def test_ensure_schema_creates_tables(self):
        """Test that schema is created correctly."""
        manager = ClusteringManager(str(self.db_path))

        # Trigger schema creation by creating a run
        run_id = manager.create_run()

        # Verify tables exist
        con = duckdb.connect(str(self.db_path))
        tables = {row[0] for row in con.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        ).fetchall()}
        con.close()

        self.assertIn("clustering_runs", tables)
        self.assertIn("versioned_clusters", tables)
        self.assertIn("versioned_cluster_members", tables)


class TestClusteringRuns(unittest.TestCase):
    """Test creating and managing clustering runs."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.db_path = Path(self.temp_dir) / "test.duckdb"

        # Create a basic database
        con = duckdb.connect(str(self.db_path))
        con.execute("""
            CREATE TABLE sra_metadata (
                sample_id VARCHAR PRIMARY KEY
            )
        """)
        con.execute("""
            CREATE TABLE sequences (
                seqhash_id VARCHAR PRIMARY KEY,
                sequence VARCHAR NOT NULL,
                sample_id VARCHAR NOT NULL,
                length INTEGER NOT NULL
            )
        """)
        con.execute("INSERT INTO sra_metadata VALUES ('sample1')")
        for i in range(10):
            con.execute(
                "INSERT INTO sequences VALUES (?, ?, 'sample1', 4)",
                [f"seq{i}", f"ACGT{i}"]
            )
        con.close()

        self.manager = ClusteringManager(str(self.db_path))

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_create_run(self):
        """Test creating a new clustering run."""
        params = ClusteringParams(min_seq_id=0.5)
        run_id = self.manager.create_run(params, notes="Test run")

        self.assertEqual(run_id, 1)

        # Verify run was created
        con = duckdb.connect(str(self.db_path))
        run = con.execute(
            "SELECT status, notes FROM clustering_runs WHERE run_id = ?",
            [run_id]
        ).fetchone()
        con.close()

        self.assertEqual(run[0], "pending")
        self.assertEqual(run[1], "Test run")

    def test_create_multiple_runs(self):
        """Test creating multiple runs increments run_id."""
        run1 = self.manager.create_run()
        run2 = self.manager.create_run()
        run3 = self.manager.create_run()

        self.assertEqual(run1, 1)
        self.assertEqual(run2, 2)
        self.assertEqual(run3, 3)

    def test_update_run_status(self):
        """Test updating run status."""
        run_id = self.manager.create_run()

        self.manager.update_run_status(run_id, "running")

        con = duckdb.connect(str(self.db_path))
        status = con.execute(
            "SELECT status FROM clustering_runs WHERE run_id = ?",
            [run_id]
        ).fetchone()[0]
        con.close()

        self.assertEqual(status, "running")

    def test_update_run_completed(self):
        """Test marking run as completed with stats."""
        run_id = self.manager.create_run()

        self.manager.update_run_status(
            run_id, "completed",
            sequence_count=100,
            cluster_count=50
        )

        con = duckdb.connect(str(self.db_path))
        run = con.execute(
            "SELECT status, sequence_count, cluster_count, completed_at FROM clustering_runs WHERE run_id = ?",
            [run_id]
        ).fetchone()
        con.close()

        self.assertEqual(run[0], "completed")
        self.assertEqual(run[1], 100)
        self.assertEqual(run[2], 50)
        self.assertIsNotNone(run[3])  # completed_at should be set


class TestStoreResults(unittest.TestCase):
    """Test storing clustering results."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.db_path = Path(self.temp_dir) / "test.duckdb"

        # Create a database with sequences
        con = duckdb.connect(str(self.db_path))
        con.execute("CREATE TABLE sra_metadata (sample_id VARCHAR PRIMARY KEY)")
        con.execute("""
            CREATE TABLE sequences (
                seqhash_id VARCHAR PRIMARY KEY,
                sequence VARCHAR NOT NULL,
                sample_id VARCHAR NOT NULL,
                length INTEGER NOT NULL
            )
        """)
        con.execute("INSERT INTO sra_metadata VALUES ('sample1')")
        for i in range(5):
            con.execute(
                "INSERT INTO sequences VALUES (?, ?, 'sample1', 4)",
                [f"seq{i}", f"ACGT"]
            )
        con.close()

        self.manager = ClusteringManager(str(self.db_path))

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_store_results(self):
        """Test storing clustering results."""
        run_id = self.manager.create_run()

        # Create a mock result
        assignments = pd.DataFrame({
            "seqhash_id": ["seq0", "seq1", "seq2", "seq3", "seq4"],
            "cluster_id": ["CLU_a", "CLU_a", "CLU_a", "CLU_b", "CLU_b"],
            "representative_seqhash_id": ["seq0", "seq0", "seq0", "seq3", "seq3"],
            "is_representative": [True, False, False, True, False]
        })

        result = ClusteringResult(
            assignments=assignments,
            sequence_count=5,
            cluster_count=2,
            params=ClusteringParams()
        )

        self.manager.store_results(run_id, result)

        # Verify clusters were stored
        con = duckdb.connect(str(self.db_path))

        clusters = con.execute(
            "SELECT cluster_id, representative_seqhash_id, size FROM versioned_clusters WHERE run_id = ?",
            [run_id]
        ).fetchall()

        self.assertEqual(len(clusters), 2)

        # Verify members were stored
        members = con.execute(
            "SELECT seqhash_id, cluster_id FROM versioned_cluster_members WHERE run_id = ?",
            [run_id]
        ).fetchall()

        self.assertEqual(len(members), 5)

        # Verify run was marked completed
        status = con.execute(
            "SELECT status FROM clustering_runs WHERE run_id = ?",
            [run_id]
        ).fetchone()[0]

        self.assertEqual(status, "completed")

        con.close()

    def test_store_results_fails_for_unsuccessful_result(self):
        """Test that storing fails for unsuccessful clustering result."""
        run_id = self.manager.create_run()

        result = ClusteringResult(
            assignments=pd.DataFrame(),
            sequence_count=0,
            cluster_count=0,
            params=ClusteringParams(),
            success=False,
            error_message="Test error"
        )

        with self.assertRaises(ValueError):
            self.manager.store_results(run_id, result)


class TestExportAndRetrieve(unittest.TestCase):
    """Test exporting sequences and retrieving results."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.db_path = Path(self.temp_dir) / "test.duckdb"

        # Create a database with sequences
        con = duckdb.connect(str(self.db_path))
        con.execute("CREATE TABLE sra_metadata (sample_id VARCHAR PRIMARY KEY)")
        con.execute("""
            CREATE TABLE sequences (
                seqhash_id VARCHAR PRIMARY KEY,
                sequence VARCHAR NOT NULL,
                sample_id VARCHAR NOT NULL,
                length INTEGER NOT NULL
            )
        """)
        con.execute("INSERT INTO sra_metadata VALUES ('sample1')")
        con.execute("INSERT INTO sequences VALUES ('seq1', 'ACGT', 'sample1', 4)")
        con.execute("INSERT INTO sequences VALUES ('seq2', 'TGCA', 'sample1', 4)")
        con.execute("INSERT INTO sequences VALUES ('seq3', 'AAAA', 'sample1', 4)")
        con.close()

        self.manager = ClusteringManager(str(self.db_path))

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_export_sequences_for_clustering(self):
        """Test exporting all sequences to FASTA."""
        output_path = Path(self.temp_dir) / "sequences.fasta"
        count = self.manager.export_sequences_for_clustering(output_path)

        self.assertEqual(count, 3)
        self.assertTrue(output_path.exists())

        content = output_path.read_text()
        self.assertIn(">seq1", content)
        self.assertIn("ACGT", content)

    def test_get_current_representatives_no_runs(self):
        """Test getting representatives when no clustering has been done."""
        reps = self.manager.get_current_representatives()
        self.assertTrue(reps.empty)

    def test_get_current_representatives_with_runs(self):
        """Test getting representatives after clustering."""
        # Create a run and store results
        run_id = self.manager.create_run()

        assignments = pd.DataFrame({
            "seqhash_id": ["seq1", "seq2", "seq3"],
            "cluster_id": ["CLU_a", "CLU_a", "CLU_b"],
            "representative_seqhash_id": ["seq1", "seq1", "seq3"],
            "is_representative": [True, False, True]
        })

        result = ClusteringResult(
            assignments=assignments,
            sequence_count=3,
            cluster_count=2,
            params=ClusteringParams()
        )

        self.manager.store_results(run_id, result)

        # Get representatives
        reps = self.manager.get_current_representatives()

        self.assertEqual(len(reps), 2)
        self.assertIn("seq1", reps["seqhash_id"].values)
        self.assertIn("seq3", reps["seqhash_id"].values)

    def test_export_representatives_fasta(self):
        """Test exporting representative sequences to FASTA."""
        # Create a run and store results
        run_id = self.manager.create_run()

        assignments = pd.DataFrame({
            "seqhash_id": ["seq1", "seq2", "seq3"],
            "cluster_id": ["CLU_a", "CLU_a", "CLU_b"],
            "representative_seqhash_id": ["seq1", "seq1", "seq3"],
            "is_representative": [True, False, True]
        })

        result = ClusteringResult(
            assignments=assignments,
            sequence_count=3,
            cluster_count=2,
            params=ClusteringParams()
        )

        self.manager.store_results(run_id, result)

        # Export representatives
        output_path = Path(self.temp_dir) / "reps.fasta"
        count = self.manager.export_representatives_fasta(output_path)

        self.assertEqual(count, 2)
        self.assertTrue(output_path.exists())

        content = output_path.read_text()
        self.assertIn(">seq1", content)
        self.assertIn(">seq3", content)
        self.assertNotIn(">seq2", content)


class TestClusteringHistory(unittest.TestCase):
    """Test clustering history and statistics."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.db_path = Path(self.temp_dir) / "test.duckdb"

        # Create a database with sequences
        con = duckdb.connect(str(self.db_path))
        con.execute("CREATE TABLE sra_metadata (sample_id VARCHAR PRIMARY KEY)")
        con.execute("""
            CREATE TABLE sequences (
                seqhash_id VARCHAR PRIMARY KEY,
                sequence VARCHAR NOT NULL,
                sample_id VARCHAR NOT NULL,
                length INTEGER NOT NULL
            )
        """)
        con.execute("INSERT INTO sra_metadata VALUES ('sample1')")
        for i in range(10):
            con.execute(
                "INSERT INTO sequences VALUES (?, ?, 'sample1', 4)",
                [f"seq{i}", f"ACGT"]
            )
        con.close()

        self.manager = ClusteringManager(str(self.db_path))

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_get_clustering_history(self):
        """Test getting clustering history."""
        # Create multiple runs
        run1 = self.manager.create_run(notes="First run")
        run2 = self.manager.create_run(notes="Second run")

        # Mark first as completed
        self.manager.update_run_status(run1, "completed", sequence_count=10, cluster_count=5)

        history = self.manager.get_clustering_history()

        self.assertEqual(len(history), 2)
        self.assertIn("run_id", history.columns)
        self.assertIn("status", history.columns)
        self.assertIn("notes", history.columns)

    def test_get_current_clustering(self):
        """Test getting current clustering assignments."""
        run_id = self.manager.create_run()

        # Create assignments for all sequences
        seq_ids = [f"seq{i}" for i in range(10)]
        assignments = pd.DataFrame({
            "seqhash_id": seq_ids,
            "cluster_id": ["CLU_a"] * 5 + ["CLU_b"] * 5,
            "representative_seqhash_id": ["seq0"] * 5 + ["seq5"] * 5,
            "is_representative": [True] + [False] * 4 + [True] + [False] * 4
        })

        result = ClusteringResult(
            assignments=assignments,
            sequence_count=10,
            cluster_count=2,
            params=ClusteringParams()
        )

        self.manager.store_results(run_id, result)

        # Get current clustering
        current = self.manager.get_current_clustering()

        self.assertEqual(len(current), 10)
        self.assertIn("seqhash_id", current.columns)
        self.assertIn("cluster_id", current.columns)
        self.assertIn("representative_seqhash_id", current.columns)

    def test_get_cluster_stats(self):
        """Test getting cluster statistics."""
        run_id = self.manager.create_run()

        # Create assignments with varied cluster sizes
        assignments = pd.DataFrame({
            "seqhash_id": [f"seq{i}" for i in range(10)],
            "cluster_id": ["CLU_a"] * 6 + ["CLU_b"] * 3 + ["CLU_c"] * 1,
            "representative_seqhash_id": ["seq0"] * 6 + ["seq6"] * 3 + ["seq9"] * 1,
            "is_representative": [True] + [False] * 5 + [True] + [False] * 2 + [True]
        })

        result = ClusteringResult(
            assignments=assignments,
            sequence_count=10,
            cluster_count=3,
            params=ClusteringParams()
        )

        self.manager.store_results(run_id, result)

        stats = self.manager.get_cluster_stats()

        self.assertEqual(stats["total_clusters"].values[0], 3)
        self.assertEqual(stats["total_sequences"].values[0], 10)
        self.assertEqual(stats["min_cluster_size"].values[0], 1)
        self.assertEqual(stats["max_cluster_size"].values[0], 6)


class TestFullClusteringWorkflow(unittest.TestCase):
    """Test the complete clustering workflow."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.db_path = Path(self.temp_dir) / "test.duckdb"

        # Create a database with sequences
        con = duckdb.connect(str(self.db_path))
        con.execute("CREATE TABLE sra_metadata (sample_id VARCHAR PRIMARY KEY)")
        con.execute("""
            CREATE TABLE sequences (
                seqhash_id VARCHAR PRIMARY KEY,
                sequence VARCHAR NOT NULL,
                sample_id VARCHAR NOT NULL,
                length INTEGER NOT NULL
            )
        """)
        con.execute("INSERT INTO sra_metadata VALUES ('sample1')")
        # Add some sequences
        sequences = [
            ("seq1", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAV"),
            ("seq2", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAV"),
            ("seq3", "MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNCN"),
        ]
        for seq_id, seq in sequences:
            con.execute(
                "INSERT INTO sequences VALUES (?, ?, 'sample1', ?)",
                [seq_id, seq, len(seq)]
            )
        con.close()

        self.manager = ClusteringManager(str(self.db_path))

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    @patch('planter.clustering.cluster.run_clustering')
    def test_run_full_clustering_mocked(self, mock_run_clustering):
        """Test full clustering workflow with mocked MMseqs2."""
        # Setup mock return
        assignments = pd.DataFrame({
            "seqhash_id": ["seq1", "seq2", "seq3"],
            "cluster_id": ["CLU_a", "CLU_a", "CLU_b"],
            "representative_seqhash_id": ["seq1", "seq1", "seq3"],
            "is_representative": [True, False, True]
        })

        mock_run_clustering.return_value = ClusteringResult(
            assignments=assignments,
            sequence_count=3,
            cluster_count=2,
            params=ClusteringParams()
        )

        # Run clustering
        run_id, result = self.manager.run_full_clustering(
            params=ClusteringParams(min_seq_id=0.9),
            notes="Test clustering"
        )

        # Verify
        self.assertTrue(result.success)
        self.assertEqual(result.cluster_count, 2)

        # Verify database state
        con = duckdb.connect(str(self.db_path))

        # Check run was created and completed
        run = con.execute(
            "SELECT status, notes FROM clustering_runs WHERE run_id = ?",
            [run_id]
        ).fetchone()
        self.assertEqual(run[0], "completed")
        self.assertEqual(run[1], "Test clustering")

        # Check clusters were stored
        clusters = con.execute(
            "SELECT COUNT(*) FROM versioned_clusters WHERE run_id = ?",
            [run_id]
        ).fetchone()[0]
        self.assertEqual(clusters, 2)

        con.close()


if __name__ == "__main__":
    unittest.main()
