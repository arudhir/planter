#!/usr/bin/env python3
"""
Tests for the DuckDB utility functions.
"""
import os
import tempfile
import unittest
from pathlib import Path
import duckdb
import pandas as pd

from planter.database.utils.duckdb_utils import merge_duckdbs, update_duckdb_with_cluster_info


class TestDuckDBUtils(unittest.TestCase):
    """Test cases for the DuckDB utility functions."""

    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        
        # Create paths for the test databases
        self.db1_path = Path(self.temp_dir) / "db1.duckdb"
        self.db2_path = Path(self.temp_dir) / "db2.duckdb"
        self.master_db_path = Path(self.temp_dir) / "master.duckdb"
        
        # Create a schema SQL file
        self.schema_sql_path = Path(self.temp_dir) / "schema.sql"
        with open(self.schema_sql_path, "w") as f:
            f.write("""
            CREATE TABLE IF NOT EXISTS sra_metadata (
                sample_id VARCHAR PRIMARY KEY,
                organism VARCHAR,
                study_title VARCHAR
            );
            
            CREATE TABLE IF NOT EXISTS sequences (
                seqhash_id VARCHAR PRIMARY KEY,
                sample_id VARCHAR,
                sequence VARCHAR,
                length INTEGER,
                description VARCHAR,
                repseq_id VARCHAR
            );
            
            CREATE TABLE IF NOT EXISTS annotations (
                seqhash_id VARCHAR,
                annotation_source VARCHAR,
                annotation_value VARCHAR,
                PRIMARY KEY (seqhash_id, annotation_source, annotation_value)
            );
            
            CREATE TABLE IF NOT EXISTS go_terms (
                seqhash_id VARCHAR,
                go_term VARCHAR,
                PRIMARY KEY (seqhash_id, go_term)
            );
            
            CREATE TABLE IF NOT EXISTS ec_numbers (
                seqhash_id VARCHAR,
                ec_number VARCHAR,
                PRIMARY KEY (seqhash_id, ec_number)
            );
            
            CREATE TABLE IF NOT EXISTS kegg_info (
                seqhash_id VARCHAR,
                kegg_id VARCHAR,
                PRIMARY KEY (seqhash_id, kegg_id)
            );
            
            CREATE TABLE IF NOT EXISTS clusters (
                cluster_id VARCHAR PRIMARY KEY,
                representative_seqhash_id VARCHAR,
                size INTEGER
            );
            
            CREATE TABLE IF NOT EXISTS cluster_members (
                seqhash_id VARCHAR PRIMARY KEY,
                cluster_id VARCHAR
            );
            
            CREATE TABLE IF NOT EXISTS expression (
                gene_seqhash_id VARCHAR,
                sample_id VARCHAR,
                tpm FLOAT,
                num_reads FLOAT,
                effective_length FLOAT,
                PRIMARY KEY (gene_seqhash_id, sample_id)
            );
            
            CREATE TABLE IF NOT EXISTS gene_protein_map (
                gene_seqhash_id VARCHAR PRIMARY KEY,
                protein_seqhash_id VARCHAR
            );
            """)
        
        # Create the first test database
        con1 = duckdb.connect(str(self.db1_path))
        con1.execute(open(self.schema_sql_path).read())
        con1.execute("INSERT INTO sra_metadata VALUES ('SRR1', 'Organism1', 'Study1');")
        con1.execute("INSERT INTO sequences VALUES ('seq1', 'SRR1', 'ACGT', 4, 'Seq1', 'seq1');")
        con1.execute("INSERT INTO annotations VALUES ('seq1', 'source1', 'value1');")
        con1.execute("INSERT INTO go_terms VALUES ('seq1', 'GO:0001');")
        con1.execute("INSERT INTO ec_numbers VALUES ('seq1', 'EC:1.1.1.1');")
        con1.execute("INSERT INTO kegg_info VALUES ('seq1', 'K00001');")
        con1.execute("INSERT INTO clusters VALUES ('cluster1', 'seq1', 1);")
        con1.execute("INSERT INTO cluster_members VALUES ('seq1', 'cluster1');")
        con1.execute("INSERT INTO expression VALUES ('seq1', 'SRR1', 10.5, 100.0, 58.12);")
        con1.execute("INSERT INTO gene_protein_map VALUES ('seq1', 'seq1.p1');")
        con1.close()
        
        # Create the second test database
        con2 = duckdb.connect(str(self.db2_path))
        con2.execute(open(self.schema_sql_path).read())
        con2.execute("INSERT INTO sra_metadata VALUES ('SRR2', 'Organism2', 'Study2');")
        con2.execute("INSERT INTO sequences VALUES ('seq2', 'SRR2', 'TGCA', 4, 'Seq2', 'seq2');")
        con2.execute("INSERT INTO annotations VALUES ('seq2', 'source2', 'value2');")
        con2.execute("INSERT INTO go_terms VALUES ('seq2', 'GO:0002');")
        con2.execute("INSERT INTO ec_numbers VALUES ('seq2', 'EC:2.2.2.2');")
        con2.execute("INSERT INTO kegg_info VALUES ('seq2', 'K00002');")
        con2.execute("INSERT INTO clusters VALUES ('cluster2', 'seq2', 1);")
        con2.execute("INSERT INTO cluster_members VALUES ('seq2', 'cluster2');")
        con2.execute("INSERT INTO expression VALUES ('seq2', 'SRR2', 5.2, 50.0, 58.12);")
        con2.execute("INSERT INTO gene_protein_map VALUES ('seq2', 'seq2.p1');")
        con2.close()
        
        # Create a TSV file for cluster info
        self.cluster_tsv_path = Path(self.temp_dir) / "clusters.tsv"
        with open(self.cluster_tsv_path, "w") as f:
            f.write("rep1\tseq1\n")
            f.write("rep1\tseq2\n")

    def tearDown(self):
        """Clean up after tests."""
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_merge_duckdbs(self):
        """Test merging DuckDB databases."""
        # Merge the databases
        merged_db_path = merge_duckdbs(
            duckdb_paths=[self.db1_path, self.db2_path],
            master_db_path=self.master_db_path,
            schema_sql_path=self.schema_sql_path
        )
        
        # Verify the merged database
        con = duckdb.connect(merged_db_path)
        
        # Check sra_metadata
        result = con.execute("SELECT count(*) FROM sra_metadata").fetchone()
        self.assertEqual(result[0], 2)
        
        # Check sequences
        result = con.execute("SELECT count(*) FROM sequences").fetchone()
        self.assertEqual(result[0], 2)
        
        # Check annotations
        result = con.execute("SELECT count(*) FROM annotations").fetchone()
        self.assertEqual(result[0], 2)
        
        # Check go_terms
        result = con.execute("SELECT count(*) FROM go_terms").fetchone()
        self.assertEqual(result[0], 2)
        
        # Check ec_numbers
        result = con.execute("SELECT count(*) FROM ec_numbers").fetchone()
        self.assertEqual(result[0], 2)
        
        # Check kegg_info
        result = con.execute("SELECT count(*) FROM kegg_info").fetchone()
        self.assertEqual(result[0], 2)
        
        # Check clusters
        result = con.execute("SELECT count(*) FROM clusters").fetchone()
        self.assertEqual(result[0], 2)
        
        # Check cluster_members
        result = con.execute("SELECT count(*) FROM cluster_members").fetchone()
        self.assertEqual(result[0], 2)
        
        # Check expression
        result = con.execute("SELECT count(*) FROM expression").fetchone()
        self.assertEqual(result[0], 2)
        
        # Check gene_protein_map
        result = con.execute("SELECT count(*) FROM gene_protein_map").fetchone()
        self.assertEqual(result[0], 2)
        
        con.close()

    def test_update_duckdb_with_cluster_info(self):
        """Test updating DuckDB with cluster info."""
        # First, create a database to update
        shutil.copy(str(self.db1_path), str(self.master_db_path))
        
        # Update the database with cluster info
        update_duckdb_with_cluster_info(self.master_db_path, self.cluster_tsv_path)
        
        # Verify the updated database
        con = duckdb.connect(str(self.master_db_path))
        
        # Check if the repseq_id in sequences table is updated
        result = con.execute("SELECT repseq_id FROM sequences WHERE seqhash_id = 'seq1'").fetchone()
        self.assertEqual(result[0], "rep1")
        
        # Check if the clusters table is updated
        result = con.execute("SELECT count(*) FROM clusters").fetchone()
        self.assertEqual(result[0], 2)  # Original cluster1 and the new rep1 cluster
        
        # Check if the cluster_members table is updated
        result = con.execute("SELECT cluster_id FROM cluster_members WHERE seqhash_id = 'seq1'").fetchone()
        self.assertEqual(result[0], "rep1")
        
        con.close()

    def test_with_real_fixture_data(self):
        """Test with real fixture data from the project."""
        import shutil
        
        # Define fixture data paths
        fixtures_dir = Path("/home/ubuntu/planter/tests/fixtures")
        src_db1 = fixtures_dir / "SRR12068547/SRR12068547.duckdb"
        src_db2 = fixtures_dir / "SRR12068548/SRR12068548.duckdb"
        
        # Copy fixture databases to temp directory
        fixture_db1 = Path(self.temp_dir) / "SRR12068547.duckdb"
        fixture_db2 = Path(self.temp_dir) / "SRR12068548.duckdb"
        fixture_master = Path(self.temp_dir) / "fixture_master.duckdb"
        
        shutil.copy(src_db1, fixture_db1)
        shutil.copy(src_db2, fixture_db2)
        
        # Get the real schema SQL path
        schema_sql_path = Path("/home/ubuntu/planter/planter/database/schema/migrations/001_initial_schema.sql")
        
        # Make sure we have the fixture files
        self.assertTrue(fixture_db1.exists(), "Fixture SRR12068547.duckdb not found")
        self.assertTrue(fixture_db2.exists(), "Fixture SRR12068548.duckdb not found")
        self.assertTrue(schema_sql_path.exists(), "Schema SQL file not found")
        
        # 1. Test merging real databases
        merged_db_path = merge_duckdbs(
            duckdb_paths=[fixture_db1, fixture_db2],
            master_db_path=fixture_master,
            schema_sql_path=schema_sql_path
        )
        
        # Verify merge results - check table counts
        con = duckdb.connect(merged_db_path)
        
        # Get counts from source databases
        con1 = duckdb.connect(str(fixture_db1))
        sample1_id = con1.execute("SELECT sample_id FROM sra_metadata LIMIT 1").fetchone()[0]
        seq_count1 = con1.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
        con1.close()
        
        con2 = duckdb.connect(str(fixture_db2))
        sample2_id = con2.execute("SELECT sample_id FROM sra_metadata LIMIT 1").fetchone()[0]
        seq_count2 = con2.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
        con2.close()
        
        # Verify merged counts
        samples = con.execute("SELECT sample_id FROM sra_metadata").fetchall()
        self.assertEqual(len(samples), 2, "Should have 2 samples after merge")
        sample_ids = [s[0] for s in samples]
        self.assertIn(sample1_id, sample_ids, f"Sample {sample1_id} should be in merged database")
        self.assertIn(sample2_id, sample_ids, f"Sample {sample2_id} should be in merged database")
        
        merged_seq_count = con.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
        self.assertEqual(merged_seq_count, seq_count1 + seq_count2, 
                        f"Expected {seq_count1 + seq_count2} sequences, got {merged_seq_count}")
        
        con.close()
        
        # 2. Test cluster update with real data
        cluster_tsv = Path(self.temp_dir) / "cluster_test.tsv"
        
        # Create a simplified cluster mapping for testing
        # Extract a few sequence IDs from the merged database
        con = duckdb.connect(merged_db_path)
        seq_ids = con.execute("SELECT seqhash_id FROM sequences LIMIT 10").fetchall()
        seq_ids = [seq[0] for seq in seq_ids]
        con.close()
        
        # Create a cluster mapping (first seq is rep for all others)
        with open(cluster_tsv, "w") as f:
            rep_id = seq_ids[0] 
            for seq_id in seq_ids:
                f.write(f"{rep_id}\t{seq_id}\n")
        
        # Update merged database with cluster info
        update_duckdb_with_cluster_info(merged_db_path, cluster_tsv)
        
        # Verify cluster update results
        con = duckdb.connect(merged_db_path)
        
        # Check that all test sequences now have the rep_id as their repseq_id
        for seq_id in seq_ids:
            repseq = con.execute(
                "SELECT repseq_id FROM sequences WHERE seqhash_id = ?", 
                [seq_id]
            ).fetchone()
            self.assertEqual(repseq[0], rep_id, 
                            f"Expected repseq_id={rep_id} for seq {seq_id}, got {repseq[0]}")
        
        # Verify the cluster was created with correct size
        cluster = con.execute(
            "SELECT size FROM clusters WHERE cluster_id = ?",
            [rep_id]
        ).fetchone()
        self.assertIsNotNone(cluster, f"Cluster {rep_id} should exist")
        self.assertEqual(cluster[0], len(seq_ids), 
                        f"Cluster should have {len(seq_ids)} members, has {cluster[0]}")
        
        con.close()
    
    def test_mmseqs_integration(self):
        """Test integration with MMSeqs2 clustering results."""
        import shutil
        import subprocess
        import tempfile
        
        # Skip if mmseqs is not installed
        try:
            subprocess.run(["which", "mmseqs"], check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.skipTest("MMSeqs2 not installed, skipping test")
        
        # Setup test environment with real peptide files
        fixtures_dir = Path("/home/ubuntu/planter/tests/fixtures")
        pep1_path = fixtures_dir / "SRR12068547/transdecoder/SRR12068547.pep"
        pep2_path = fixtures_dir / "SRR12068548/transdecoder/SRR12068548.pep"
        
        self.assertTrue(pep1_path.exists(), "Peptide file 1 not found")
        self.assertTrue(pep2_path.exists(), "Peptide file 2 not found")
        
        # Create a mmseqs work directory
        mmseqs_dir = Path(self.temp_dir) / "mmseqs_test"
        mmseqs_dir.mkdir()
        
        # Create a concatenated peptide file
        concat_pep = mmseqs_dir / "combined.pep"
        with open(concat_pep, "wb") as outfile:
            for pep_file in [pep1_path, pep2_path]:
                with open(pep_file, "rb") as infile:
                    shutil.copyfileobj(infile, outfile)
        
        # Run mmseqs2 clustering
        try:
            # Create sequence database
            subprocess.run([
                "mmseqs", "createdb", str(concat_pep), 
                str(mmseqs_dir / "seqDB")
            ], check=True, capture_output=True)
            
            # Cluster sequences
            subprocess.run([
                "mmseqs", "cluster", str(mmseqs_dir / "seqDB"), 
                str(mmseqs_dir / "clusterDB"), str(mmseqs_dir),
                "--min-seq-id", "0.9", "--threads", "1"
            ], check=True, capture_output=True)
            
            # Create tsv from clustering results
            subprocess.run([
                "mmseqs", "createtsv", str(mmseqs_dir / "seqDB"), 
                str(mmseqs_dir / "seqDB"), str(mmseqs_dir / "clusterDB"),
                str(mmseqs_dir / "clusters.tsv")
            ], check=True, capture_output=True)
            
            # Check that tsv file exists and has content
            tsv_path = mmseqs_dir / "clusters.tsv"
            self.assertTrue(tsv_path.exists(), "Cluster TSV file not created")
            self.assertGreater(tsv_path.stat().st_size, 0, "Cluster TSV file is empty")
            
            # Create a test database and populate it with some of the sequences
            db_path = mmseqs_dir / "test.duckdb"
            con = duckdb.connect(str(db_path))
            con.execute(open(self.schema_sql_path).read())
            
            # Get a few sample sequence IDs from the pep files
            seq_ids = []
            with open(pep1_path, "r") as f:
                for line in f:
                    if line.startswith(">"):
                        seq_id = line.strip().lstrip(">").split()[0]
                        seq_ids.append(seq_id)
                        if len(seq_ids) >= 5:
                            break
            
            # Insert test sequences
            for i, seq_id in enumerate(seq_ids):
                con.execute(
                    "INSERT INTO sequences VALUES (?, ?, ?, ?, ?, ?)",
                    [seq_id, f"SRR{i}", "ACGT", 4, f"Seq{i}", seq_id]
                )
            
            con.close()
            
            # Update database with clustering results
            update_duckdb_with_cluster_info(db_path, tsv_path)
            
            # Verify database was updated with clustering info
            con = duckdb.connect(str(db_path))
            
            # Check that clusters were created
            clusters = con.execute("SELECT COUNT(*) FROM clusters").fetchone()[0]
            self.assertGreater(clusters, 0, "No clusters were created")
            
            # Check that cluster members were added
            members = con.execute("SELECT COUNT(*) FROM cluster_members").fetchone()[0]
            self.assertGreater(members, 0, "No cluster members were added")
            
            # Check repseq_id updates - at least some of our sequences should have
            # their repseq_id different from their seqhash_id after clustering
            con.close()
            
        except subprocess.CalledProcessError as e:
            self.fail(f"MMSeqs2 command failed: {e.stderr}")


if __name__ == "__main__":
    unittest.main()