#!/usr/bin/env python3
"""
Tests for the pure clustering functions.
"""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch, MagicMock

import pandas as pd

from planter.clustering.cluster import (
    ClusteringParams,
    ClusteringResult,
    run_clustering,
    export_sequences_to_fasta,
    _generate_cluster_id,
    _parse_cluster_tsv,
)


class TestClusteringParams(unittest.TestCase):
    """Test cases for ClusteringParams dataclass."""

    def test_default_params(self):
        """Test default parameter values."""
        params = ClusteringParams()
        self.assertEqual(params.min_seq_id, 0.3)
        self.assertEqual(params.coverage, 0.8)
        self.assertEqual(params.coverage_mode, 0)
        self.assertEqual(params.cluster_mode, 0)
        self.assertEqual(params.threads, 4)

    def test_custom_params(self):
        """Test custom parameter values."""
        params = ClusteringParams(
            min_seq_id=0.5,
            coverage=0.9,
            threads=8
        )
        self.assertEqual(params.min_seq_id, 0.5)
        self.assertEqual(params.coverage, 0.9)
        self.assertEqual(params.threads, 8)

    def test_to_mmseqs_args(self):
        """Test conversion to MMseqs2 command line arguments."""
        params = ClusteringParams(min_seq_id=0.5, coverage=0.9, threads=8)
        args = params.to_mmseqs_args()

        self.assertIn("--min-seq-id", args)
        self.assertIn("0.5", args)
        self.assertIn("-c", args)
        self.assertIn("0.9", args)
        self.assertIn("--threads", args)
        self.assertIn("8", args)

    def test_to_dict(self):
        """Test conversion to dictionary."""
        params = ClusteringParams(min_seq_id=0.5, coverage=0.9)
        d = params.to_dict()

        self.assertEqual(d["min_seq_id"], 0.5)
        self.assertEqual(d["coverage"], 0.9)
        self.assertIn("threads", d)


class TestClusteringResult(unittest.TestCase):
    """Test cases for ClusteringResult dataclass."""

    def test_successful_result(self):
        """Test a successful clustering result."""
        assignments = pd.DataFrame({
            "seqhash_id": ["seq1", "seq2", "seq3"],
            "cluster_id": ["CLU_abc", "CLU_abc", "CLU_def"],
            "representative_seqhash_id": ["seq1", "seq1", "seq3"],
            "is_representative": [True, False, True]
        })

        result = ClusteringResult(
            assignments=assignments,
            sequence_count=3,
            cluster_count=2,
            params=ClusteringParams(),
            success=True
        )

        self.assertTrue(result.success)
        self.assertEqual(result.sequence_count, 3)
        self.assertEqual(result.cluster_count, 2)
        self.assertIsNone(result.error_message)

    def test_failed_result(self):
        """Test a failed clustering result."""
        result = ClusteringResult(
            assignments=pd.DataFrame(),
            sequence_count=0,
            cluster_count=0,
            params=ClusteringParams(),
            success=False,
            error_message="MMseqs2 not found"
        )

        self.assertFalse(result.success)
        self.assertEqual(result.error_message, "MMseqs2 not found")

    def test_representatives_property(self):
        """Test the representatives property."""
        assignments = pd.DataFrame({
            "seqhash_id": ["seq1", "seq2", "seq3", "seq4"],
            "cluster_id": ["CLU_a", "CLU_a", "CLU_b", "CLU_b"],
            "representative_seqhash_id": ["seq1", "seq1", "seq3", "seq3"],
            "is_representative": [True, False, True, False]
        })

        result = ClusteringResult(
            assignments=assignments,
            sequence_count=4,
            cluster_count=2,
            params=ClusteringParams()
        )

        reps = result.representatives
        self.assertEqual(len(reps), 2)
        self.assertIn("seq1", reps["seqhash_id"].values)
        self.assertIn("seq3", reps["seqhash_id"].values)


class TestHelperFunctions(unittest.TestCase):
    """Test helper functions."""

    def test_generate_cluster_id(self):
        """Test cluster ID generation."""
        cluster_id = _generate_cluster_id("test_sequence_123")

        # Should start with CLU_ prefix
        self.assertTrue(cluster_id.startswith("CLU_"))

        # Should be deterministic
        self.assertEqual(
            _generate_cluster_id("test_sequence_123"),
            _generate_cluster_id("test_sequence_123")
        )

        # Different inputs should produce different IDs
        self.assertNotEqual(
            _generate_cluster_id("seq1"),
            _generate_cluster_id("seq2")
        )

    def test_parse_cluster_tsv(self):
        """Test parsing MMseqs2 cluster TSV output."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("rep1\tseq1\n")
            f.write("rep1\tseq2\n")
            f.write("rep1\trep1\n")
            f.write("rep2\tseq3\n")
            f.write("rep2\trep2\n")
            tsv_path = Path(f.name)

        try:
            df = _parse_cluster_tsv(tsv_path)

            # Check structure
            self.assertIn("seqhash_id", df.columns)
            self.assertIn("cluster_id", df.columns)
            self.assertIn("representative_seqhash_id", df.columns)
            self.assertIn("is_representative", df.columns)

            # Check counts
            self.assertEqual(len(df), 5)
            self.assertEqual(df["cluster_id"].nunique(), 2)

            # Check representative flags
            reps = df[df["is_representative"]]
            self.assertEqual(len(reps), 2)
            self.assertIn("rep1", reps["seqhash_id"].values)
            self.assertIn("rep2", reps["seqhash_id"].values)

        finally:
            tsv_path.unlink()


class TestExportSequences(unittest.TestCase):
    """Test sequence export functionality."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_export_sequences_to_fasta(self):
        """Test exporting DataFrame to FASTA format."""
        sequences = pd.DataFrame({
            "seqhash_id": ["seq1", "seq2", "seq3"],
            "sequence": ["ACGT", "TGCA", "AAAA"]
        })

        output_path = Path(self.temp_dir) / "output.fasta"
        count = export_sequences_to_fasta(sequences, output_path)

        self.assertEqual(count, 3)
        self.assertTrue(output_path.exists())

        # Verify content
        content = output_path.read_text()
        self.assertIn(">seq1", content)
        self.assertIn("ACGT", content)
        self.assertIn(">seq2", content)
        self.assertIn("TGCA", content)

    def test_export_creates_parent_dirs(self):
        """Test that export creates parent directories."""
        sequences = pd.DataFrame({
            "seqhash_id": ["seq1"],
            "sequence": ["ACGT"]
        })

        output_path = Path(self.temp_dir) / "subdir" / "nested" / "output.fasta"
        count = export_sequences_to_fasta(sequences, output_path)

        self.assertEqual(count, 1)
        self.assertTrue(output_path.exists())


class TestRunClustering(unittest.TestCase):
    """Test the main clustering function."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_missing_input_file(self):
        """Test handling of missing input file."""
        result = run_clustering(
            fasta_path=Path("/nonexistent/file.fasta"),
            params=ClusteringParams()
        )

        self.assertFalse(result.success)
        self.assertIn("not found", result.error_message)

    def test_empty_fasta(self):
        """Test handling of empty FASTA file."""
        fasta_path = Path(self.temp_dir) / "empty.fasta"
        fasta_path.touch()

        result = run_clustering(
            fasta_path=fasta_path,
            params=ClusteringParams(),
            work_dir=Path(self.temp_dir)
        )

        # Should fail or return empty result
        # (exact behavior depends on MMseqs2 handling of empty input)
        if result.success:
            self.assertEqual(result.sequence_count, 0)

    @patch('planter.clustering.cluster._run_mmseqs_command')
    def test_clustering_workflow_mocked(self, mock_run):
        """Test clustering workflow with mocked MMseqs2."""
        # Create a test FASTA
        fasta_path = Path(self.temp_dir) / "test.fasta"
        with open(fasta_path, 'w') as f:
            f.write(">seq1\nACGT\n")
            f.write(">seq2\nTGCA\n")
            f.write(">seq3\nAAAA\n")

        # Mock successful MMseqs2 runs
        mock_run.return_value = True

        # Create expected TSV output
        work_dir = Path(self.temp_dir) / "work"
        work_dir.mkdir()
        tsv_path = work_dir / "clusters.tsv"
        with open(tsv_path, 'w') as f:
            f.write("seq1\tseq1\n")
            f.write("seq1\tseq2\n")
            f.write("seq3\tseq3\n")

        result = run_clustering(
            fasta_path=fasta_path,
            params=ClusteringParams(),
            work_dir=work_dir,
            keep_temp=True
        )

        # Verify MMseqs2 commands were called
        self.assertTrue(mock_run.called)


class TestClusteringIntegration(unittest.TestCase):
    """Integration tests that require MMseqs2 to be installed."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

        # Check if MMseqs2 is available
        import subprocess
        try:
            subprocess.run(["mmseqs", "version"], capture_output=True, check=True)
            self.mmseqs_available = True
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.mmseqs_available = False

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_full_clustering_workflow(self):
        """Test full clustering workflow with real MMseqs2."""
        if not self.mmseqs_available:
            self.skipTest("MMseqs2 not available")

        # Create test sequences
        fasta_path = Path(self.temp_dir) / "test_sequences.fasta"
        with open(fasta_path, 'w') as f:
            # Create some sequences that should cluster together
            f.write(">seq1\n")
            f.write("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAV\n")
            f.write(">seq2\n")
            f.write("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAV\n")  # Same as seq1
            f.write(">seq3\n")
            f.write("MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNCN\n")  # Different
            f.write(">seq4\n")
            f.write("MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNCN\n")  # Same as seq3

        # Run clustering
        result = run_clustering(
            fasta_path=fasta_path,
            params=ClusteringParams(min_seq_id=0.9, coverage=0.8),
            work_dir=Path(self.temp_dir) / "mmseqs_work"
        )

        # Verify results
        self.assertTrue(result.success, f"Clustering failed: {result.error_message}")
        self.assertEqual(result.sequence_count, 4)
        self.assertEqual(result.cluster_count, 2)  # Should have 2 clusters

        # Verify assignments
        self.assertEqual(len(result.assignments), 4)
        self.assertEqual(result.assignments["is_representative"].sum(), 2)


if __name__ == "__main__":
    unittest.main()
