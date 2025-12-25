#!/usr/bin/env python3
"""
Unit tests for MMSeqs2 representative sequence selection.

Tests verify that clustering with --cluster-mode 2 and --cov-mode 1
correctly selects the longest sequence as the cluster representative.
"""
import pytest
import subprocess
import tempfile
import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class TestMMSeqsRepresentativeSelection:
    """Test MMSeqs2 representative selection with synthetic data."""

    def create_test_sequences(self, sequences):
        """
        Create a FASTA file with test sequences.

        Args:
            sequences: List of tuples (id, sequence)

        Returns:
            Path to temporary FASTA file
        """
        fasta_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)

        for seq_id, seq in sequences:
            fasta_file.write(f">{seq_id}\n{seq}\n")

        fasta_file.close()
        return fasta_file.name

    def run_mmseqs_cluster(self, fasta_file, cluster_mode=2, cov_mode=1, min_seq_id=0.9):
        """
        Run MMSeqs2 clustering on a FASTA file.

        Returns:
            Dictionary mapping cluster_rep -> [member_ids]
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create databases
            seqdb = os.path.join(tmpdir, "seqDB")
            clusterdb = os.path.join(tmpdir, "clusterDB")
            cluster_tsv = os.path.join(tmpdir, "cluster.tsv")

            # Create sequence database
            subprocess.run(
                ["mmseqs", "createdb", fasta_file, seqdb],
                check=True,
                capture_output=True
            )

            # Run clustering
            cmd = [
                "mmseqs", "cluster",
                seqdb, clusterdb, tmpdir,
                "--min-seq-id", str(min_seq_id),
                "-c", "0.8",
                "--cluster-mode", str(cluster_mode),
            ]

            if cov_mode is not None:
                cmd.extend(["--cov-mode", str(cov_mode)])

            subprocess.run(cmd, check=True, capture_output=True)

            # Create TSV output
            subprocess.run(
                ["mmseqs", "createtsv", seqdb, seqdb, clusterdb, cluster_tsv],
                check=True,
                capture_output=True
            )

            # Parse clusters
            clusters = {}
            with open(cluster_tsv, 'r') as f:
                for line in f:
                    rep, member = line.strip().split('\t')
                    if rep not in clusters:
                        clusters[rep] = []
                    clusters[rep].append(member)

            return clusters

    def test_simple_truncation_case(self):
        """
        Test that full-length sequence is chosen over truncated version.

        Setup: Two identical sequences except one is truncated
        Expected: Full-length sequence should be representative
        """
        print("\n" + "="*80)
        print("TEST: Simple truncation case")
        print("="*80)

        # Full-length sequence (100 aa)
        full_seq = "M" + "ACDEFGHIKLMNPQRSTVWY" * 5  # 100 residues

        # Truncated version (50 aa) - first half
        truncated_seq = full_seq[:50]

        sequences = [
            ("full_length", full_seq),
            ("truncated", truncated_seq),
        ]

        print(f"\nSequences:")
        print(f"  full_length: {len(full_seq)} aa")
        print(f"  truncated: {len(truncated_seq)} aa")

        # Test with cluster-mode 2, cov-mode 1
        fasta = self.create_test_sequences(sequences)
        try:
            clusters = self.run_mmseqs_cluster(fasta, cluster_mode=2, cov_mode=1)

            print(f"\nClusters (mode=2, cov=1):")
            for rep, members in clusters.items():
                print(f"  Representative: {rep}")
                print(f"  Members: {members}")

            # Should have 1 cluster with full_length as representative
            assert len(clusters) == 1, f"Expected 1 cluster, got {len(clusters)}"

            rep = list(clusters.keys())[0]
            assert rep == "full_length", f"Expected 'full_length' as rep, got '{rep}'"
            print("\n✓ PASS: Full-length sequence correctly selected as representative")

        finally:
            os.unlink(fasta)

    def test_multiple_lengths_in_cluster(self):
        """
        Test with multiple sequences of different lengths in same cluster.

        Setup: 5 sequences: 500, 400, 300, 200, 100 aa (all similar)
        Expected: Longest (500 aa) should be representative
        """
        print("\n" + "="*80)
        print("TEST: Multiple lengths in cluster")
        print("="*80)

        # Create base sequence
        base = "ACDEFGHIKLMNPQRSTVWY"

        sequences = [
            ("seq_500", "M" + base * 25),  # 501 aa
            ("seq_400", "M" + base * 20),  # 401 aa
            ("seq_300", "M" + base * 15),  # 301 aa
            ("seq_200", "M" + base * 10),  # 201 aa
            ("seq_100", "M" + base * 5),   # 101 aa
        ]

        print(f"\nSequences:")
        for seq_id, seq in sequences:
            print(f"  {seq_id}: {len(seq)} aa")

        fasta = self.create_test_sequences(sequences)
        try:
            clusters = self.run_mmseqs_cluster(fasta, cluster_mode=2, cov_mode=1, min_seq_id=0.8)

            print(f"\nClusters (mode=2, cov=1):")
            for rep, members in clusters.items():
                print(f"  Representative: {rep} ({len(members)} members)")
                print(f"  Members: {members}")

            # Should have 1 cluster with seq_500 as representative
            assert len(clusters) == 1, f"Expected 1 cluster, got {len(clusters)}"

            rep = list(clusters.keys())[0]
            assert rep == "seq_500", f"Expected 'seq_500' as rep, got '{rep}'"
            print("\n✓ PASS: Longest sequence (500 aa) correctly selected as representative")

        finally:
            os.unlink(fasta)

    def test_without_cov_mode(self):
        """
        Test clustering WITHOUT cov-mode 1 to see different behavior.

        This demonstrates why cov-mode 1 is necessary.
        """
        print("\n" + "="*80)
        print("TEST: Clustering WITHOUT cov-mode 1 (for comparison)")
        print("="*80)

        base = "ACDEFGHIKLMNPQRSTVWY"

        sequences = [
            ("seq_500", "M" + base * 25),
            ("seq_400", "M" + base * 20),
            ("seq_300", "M" + base * 15),
            ("seq_200", "M" + base * 10),
            ("seq_100", "M" + base * 5),
        ]

        print(f"\nSequences:")
        for seq_id, seq in sequences:
            print(f"  {seq_id}: {len(seq)} aa")

        fasta = self.create_test_sequences(sequences)
        try:
            # Test with cluster-mode 2 but NO cov-mode
            clusters = self.run_mmseqs_cluster(fasta, cluster_mode=2, cov_mode=None, min_seq_id=0.8)

            print(f"\nClusters (mode=2, NO cov-mode):")
            for rep, members in clusters.items():
                print(f"  Representative: {rep} ({len(members)} members)")

            rep = list(clusters.keys())[0]
            if rep != "seq_500":
                print(f"\n⚠ WARNING: Without cov-mode 1, representative is '{rep}' not 'seq_500'")
                print("This demonstrates why cov-mode 1 is needed!")
            else:
                print(f"\n✓ Representative is '{rep}' (may vary)")

        finally:
            os.unlink(fasta)

    def test_rhodiola_like_scenario(self):
        """
        Test scenario similar to the Rhodiola UGT problem.

        Setup: Cluster with 5 full-length (498 aa) and 1 truncated (241 aa)
        Expected: One of the full-length should be representative, NOT truncated
        """
        print("\n" + "="*80)
        print("TEST: Rhodiola-like scenario")
        print("="*80)

        # Simulate UGT sequences
        base = "MSLIEKPLTAIETREKPHAVCIPYPAQGHINPMMQLAKLLHHSGFH"
        full_length = base * 10 + "ENDING"  # ~480 aa

        sequences = [
            ("UGT_full_1", full_length),
            ("UGT_full_2", full_length.replace("MSLIEK", "MSLIKK")),  # Slight variation
            ("UGT_full_3", full_length.replace("ENDING", "FINISH")),
            ("UGT_full_4", full_length),
            ("UGT_full_5", full_length.replace("PLTAIE", "PLTAIV")),
            ("UGT_truncated", full_length[:241]),  # Truncated at 241 aa
        ]

        print(f"\nSequences:")
        for seq_id, seq in sequences:
            print(f"  {seq_id}: {len(seq)} aa")

        fasta = self.create_test_sequences(sequences)
        try:
            clusters = self.run_mmseqs_cluster(fasta, cluster_mode=2, cov_mode=1, min_seq_id=0.9)

            print(f"\nClusters (mode=2, cov=1):")
            for rep, members in clusters.items():
                print(f"  Representative: {rep} ({len(members)} members)")

                # Get length of representative
                rep_seq = [seq for sid, seq in sequences if sid == rep][0]
                print(f"  Representative length: {len(rep_seq)} aa")

            # Verify truncated is NOT representative
            for rep in clusters.keys():
                assert not rep.startswith("UGT_truncated"), \
                    f"FAIL: Truncated sequence became representative!"

            print("\n✓ PASS: Truncated sequence NOT selected as representative")

        finally:
            os.unlink(fasta)

    def test_cluster_mode_comparison(self):
        """
        Compare different cluster modes to understand their behavior.

        Tests cluster-mode 0, 2, and 3 to see differences.
        """
        print("\n" + "="*80)
        print("TEST: Cluster mode comparison")
        print("="*80)

        base = "ACDEFGHIKLMNPQRSTVWY"
        sequences = [
            ("seq_500", "M" + base * 25),
            ("seq_300", "M" + base * 15),
            ("seq_100", "M" + base * 5),
        ]

        print(f"\nSequences:")
        for seq_id, seq in sequences:
            print(f"  {seq_id}: {len(seq)} aa")

        fasta = self.create_test_sequences(sequences)

        try:
            for mode in [0, 2, 3]:
                print(f"\n--- Cluster mode {mode} ---")
                clusters = self.run_mmseqs_cluster(fasta, cluster_mode=mode, cov_mode=1, min_seq_id=0.8)

                for rep, members in clusters.items():
                    print(f"  Representative: {rep}")
                    print(f"  Members: {members}")

        finally:
            os.unlink(fasta)

    def test_coverage_threshold_impact(self):
        """
        Test how coverage threshold impacts clustering with different length sequences.

        This tests the hypothesis that high coverage threshold with cov-mode 1
        prevents clustering of different-length sequences.
        """
        print("\n" + "="*80)
        print("TEST: Coverage threshold impact on clustering")
        print("="*80)

        base = "ACDEFGHIKLMNPQRSTVWY"
        sequences = [
            ("seq_500", "M" + base * 25),  # 501 aa
            ("seq_400", "M" + base * 20),  # 401 aa
            ("seq_300", "M" + base * 15),  # 301 aa
            ("seq_200", "M" + base * 10),  # 201 aa
            ("seq_100", "M" + base * 5),   # 101 aa
        ]

        print(f"\nSequences:")
        for seq_id, seq in sequences:
            print(f"  {seq_id}: {len(seq)} aa")

        fasta = self.create_test_sequences(sequences)

        try:
            # Test 1: cluster-mode 2, cov-mode 1, DEFAULT coverage (0.8)
            print(f"\n--- Test 1: cluster-mode 2, cov-mode 1, -c 0.8 (DEFAULT) ---")

            with tempfile.TemporaryDirectory() as tmpdir:
                seqdb = os.path.join(tmpdir, "seqDB")
                clusterdb = os.path.join(tmpdir, "clusterDB")
                cluster_tsv = os.path.join(tmpdir, "cluster.tsv")

                subprocess.run(["mmseqs", "createdb", fasta, seqdb], check=True, capture_output=True)

                cmd = [
                    "mmseqs", "cluster", seqdb, clusterdb, tmpdir,
                    "--cluster-mode", "2",
                    "--cov-mode", "1",
                    "-c", "0.8",  # Default
                    "--min-seq-id", "0.8"
                ]
                subprocess.run(cmd, check=True, capture_output=True)
                subprocess.run(["mmseqs", "createtsv", seqdb, seqdb, clusterdb, cluster_tsv], check=True, capture_output=True)

                clusters_high_cov = {}
                with open(cluster_tsv, 'r') as f:
                    for line in f:
                        rep, member = line.strip().split('\t')
                        if rep not in clusters_high_cov:
                            clusters_high_cov[rep] = []
                        clusters_high_cov[rep].append(member)

            print(f"Clusters with -c 0.8:")
            for rep, members in clusters_high_cov.items():
                print(f"  Representative: {rep} ({len(members)} members)")
                print(f"    Members: {members}")

            # Test 2: cluster-mode 2, cov-mode 1, LOW coverage (0.0)
            print(f"\n--- Test 2: cluster-mode 2, cov-mode 1, -c 0.0 (PERMISSIVE) ---")

            with tempfile.TemporaryDirectory() as tmpdir:
                seqdb = os.path.join(tmpdir, "seqDB")
                clusterdb = os.path.join(tmpdir, "clusterDB")
                cluster_tsv = os.path.join(tmpdir, "cluster.tsv")

                subprocess.run(["mmseqs", "createdb", fasta, seqdb], check=True, capture_output=True)

                cmd = [
                    "mmseqs", "cluster", seqdb, clusterdb, tmpdir,
                    "--cluster-mode", "2",
                    "--cov-mode", "1",
                    "-c", "0.0",  # Permissive
                    "--min-seq-id", "0.8"
                ]
                subprocess.run(cmd, check=True, capture_output=True)
                subprocess.run(["mmseqs", "createtsv", seqdb, seqdb, clusterdb, cluster_tsv], check=True, capture_output=True)

                clusters_low_cov = {}
                with open(cluster_tsv, 'r') as f:
                    for line in f:
                        rep, member = line.strip().split('\t')
                        if rep not in clusters_low_cov:
                            clusters_low_cov[rep] = []
                        clusters_low_cov[rep].append(member)

            print(f"Clusters with -c 0.0:")
            for rep, members in clusters_low_cov.items():
                print(f"  Representative: {rep} ({len(members)} members)")
                print(f"    Members: {members}")

            # Analysis
            print(f"\n{'='*80}")
            print("ANALYSIS")
            print(f"{'='*80}")

            print(f"\nWith -c 0.8 (high coverage threshold):")
            print(f"  Total clusters: {len(clusters_high_cov)}")

            print(f"\nWith -c 0.0 (no coverage requirement):")
            print(f"  Total clusters: {len(clusters_low_cov)}")

            if len(clusters_low_cov) < len(clusters_high_cov):
                print(f"\n✓ CONFIRMED: Lower coverage (-c 0.0) creates fewer, larger clusters")
                print(f"  This allows sequences of different lengths to cluster together")

            # Check if longest is representative in both cases
            for label, clusters in [("High coverage", clusters_high_cov), ("Low coverage", clusters_low_cov)]:
                if len(clusters) == 1:
                    rep = list(clusters.keys())[0]
                    if rep == "seq_500":
                        print(f"\n✓ {label}: Longest sequence (seq_500) is representative")
                    else:
                        print(f"\n✗ {label}: Representative is {rep}, not seq_500!")

        finally:
            os.unlink(fasta)


def test_suite():
    """Run all tests."""
    tester = TestMMSeqsRepresentativeSelection()

    print("\n" + "="*80)
    print("MMSeqs2 Representative Selection Test Suite")
    print("="*80)

    tester.test_simple_truncation_case()
    tester.test_multiple_lengths_in_cluster()
    tester.test_without_cov_mode()
    tester.test_rhodiola_like_scenario()
    tester.test_cluster_mode_comparison()
    tester.test_coverage_threshold_impact()

    print("\n" + "="*80)
    print("All tests completed!")
    print("="*80)


    def test_coverage_threshold_impact(self):
        """
        Test how coverage threshold impacts clustering with different length sequences.

        This tests the hypothesis that high coverage threshold with cov-mode 1
        prevents clustering of different-length sequences.
        """
        print("\n" + "="*80)
        print("TEST: Coverage threshold impact on clustering")
        print("="*80)

        base = "ACDEFGHIKLMNPQRSTVWY"
        sequences = [
            ("seq_500", "M" + base * 25),  # 501 aa
            ("seq_400", "M" + base * 20),  # 401 aa
            ("seq_300", "M" + base * 15),  # 301 aa
            ("seq_200", "M" + base * 10),  # 201 aa
            ("seq_100", "M" + base * 5),   # 101 aa
        ]

        print(f"\nSequences:")
        for seq_id, seq in sequences:
            print(f"  {seq_id}: {len(seq)} aa")

        fasta = self.create_test_sequences(sequences)

        try:
            # Test 1: cluster-mode 2, cov-mode 1, DEFAULT coverage (0.8)
            print(f"\n--- Test 1: cluster-mode 2, cov-mode 1, -c 0.8 (DEFAULT) ---")

            with tempfile.TemporaryDirectory() as tmpdir:
                seqdb = os.path.join(tmpdir, "seqDB")
                clusterdb = os.path.join(tmpdir, "clusterDB")
                cluster_tsv = os.path.join(tmpdir, "cluster.tsv")

                subprocess.run(["mmseqs", "createdb", fasta, seqdb], check=True, capture_output=True)

                cmd = [
                    "mmseqs", "cluster", seqdb, clusterdb, tmpdir,
                    "--cluster-mode", "2",
                    "--cov-mode", "1",
                    "-c", "0.8",  # Default
                    "--min-seq-id", "0.8"
                ]
                subprocess.run(cmd, check=True, capture_output=True)
                subprocess.run(["mmseqs", "createtsv", seqdb, seqdb, clusterdb, cluster_tsv], check=True, capture_output=True)

                clusters_high_cov = {}
                with open(cluster_tsv, 'r') as f:
                    for line in f:
                        rep, member = line.strip().split('\t')
                        if rep not in clusters_high_cov:
                            clusters_high_cov[rep] = []
                        clusters_high_cov[rep].append(member)

            print(f"Clusters with -c 0.8:")
            for rep, members in clusters_high_cov.items():
                print(f"  Representative: {rep} ({len(members)} members)")
                print(f"    Members: {members}")

            # Test 2: cluster-mode 2, cov-mode 1, LOW coverage (0.0)
            print(f"\n--- Test 2: cluster-mode 2, cov-mode 1, -c 0.0 (PERMISSIVE) ---")

            with tempfile.TemporaryDirectory() as tmpdir:
                seqdb = os.path.join(tmpdir, "seqDB")
                clusterdb = os.path.join(tmpdir, "clusterDB")
                cluster_tsv = os.path.join(tmpdir, "cluster.tsv")

                subprocess.run(["mmseqs", "createdb", fasta, seqdb], check=True, capture_output=True)

                cmd = [
                    "mmseqs", "cluster", seqdb, clusterdb, tmpdir,
                    "--cluster-mode", "2",
                    "--cov-mode", "1",
                    "-c", "0.0",  # Permissive
                    "--min-seq-id", "0.8"
                ]
                subprocess.run(cmd, check=True, capture_output=True)
                subprocess.run(["mmseqs", "createtsv", seqdb, seqdb, clusterdb, cluster_tsv], check=True, capture_output=True)

                clusters_low_cov = {}
                with open(cluster_tsv, 'r') as f:
                    for line in f:
                        rep, member = line.strip().split('\t')
                        if rep not in clusters_low_cov:
                            clusters_low_cov[rep] = []
                        clusters_low_cov[rep].append(member)

            print(f"Clusters with -c 0.0:")
            for rep, members in clusters_low_cov.items():
                print(f"  Representative: {rep} ({len(members)} members)")
                print(f"    Members: {members}")

            # Analysis
            print(f"\n{'='*80}")
            print("ANALYSIS")
            print(f"{'='*80}")

            print(f"\nWith -c 0.8 (high coverage threshold):")
            print(f"  Total clusters: {len(clusters_high_cov)}")

            print(f"\nWith -c 0.0 (no coverage requirement):")
            print(f"  Total clusters: {len(clusters_low_cov)}")

            if len(clusters_low_cov) < len(clusters_high_cov):
                print(f"\n✓ CONFIRMED: Lower coverage (-c 0.0) creates fewer, larger clusters")
                print(f"  This allows sequences of different lengths to cluster together")

            # Check if longest is representative in both cases
            for label, clusters in [("High coverage", clusters_high_cov), ("Low coverage", clusters_low_cov)]:
                if len(clusters) == 1:
                    rep = list(clusters.keys())[0]
                    if rep == "seq_500":
                        print(f"\n✓ {label}: Longest sequence (seq_500) is representative")
                    else:
                        print(f"\n✗ {label}: Representative is {rep}, not seq_500!")

        finally:
            os.unlink(fasta)


if __name__ == "__main__":
    test_suite()
