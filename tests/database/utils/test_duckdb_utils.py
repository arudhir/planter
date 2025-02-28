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


if __name__ == "__main__":
    unittest.main()