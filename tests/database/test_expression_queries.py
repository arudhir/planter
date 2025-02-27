#!/usr/bin/env python3
"""
Tests for the expression-related queries in the database.
"""
import os
import tempfile
import unittest
from pathlib import Path
import json
import duckdb
import pandas as pd
from unittest.mock import patch

from planter.database.query_manager import QueryManager
from planter.database.builder import SequenceDBBuilder


class TestExpressionQueries(unittest.TestCase):
    """Test cases for expression-related database queries."""
    
    def setUp(self):
        """Set up a test database with expression data."""
        # Create a temporary database
        self.temp_dir = tempfile.mkdtemp()
        self.db_path = Path(self.temp_dir) / "test_expr.duckdb"
        
        # Create a connection to the database
        self.con = duckdb.connect(str(self.db_path))
        
        # Set up the database schema
        self.sample_id = 'SRR12068547'
        self._create_test_database()
        
        # Create the query manager
        query_dir = Path(__file__).parents[2] / 'planter' / 'database' / 'queries' / 'sql'
        self.query_manager = QueryManager(self.con, str(query_dir))
    
    def tearDown(self):
        """Clean up after tests."""
        self.con.close()
        os.remove(self.db_path)
        os.rmdir(self.temp_dir)
    
    def _create_test_database(self):
        """Create a test database with schema and sample data."""
        # Create the schema
        schema_sql = """
        -- SRA metadata table
        CREATE TABLE sra_metadata (
            sample_id VARCHAR PRIMARY KEY,
            organism VARCHAR,
            study_title VARCHAR,
            study_abstract VARCHAR,
            bioproject VARCHAR,
            biosample VARCHAR,
            library_strategy VARCHAR,
            library_source VARCHAR,
            library_selection VARCHAR,
            library_layout VARCHAR,
            instrument VARCHAR,
            run_spots INTEGER,
            run_bases INTEGER,
            run_published TIMESTAMP
        );
        
        -- Sequences table
        CREATE TABLE sequences (
            seqhash_id VARCHAR PRIMARY KEY,
            sequence VARCHAR NOT NULL,
            sample_id VARCHAR NOT NULL,
            assembly_date TIMESTAMP,
            is_representative BOOLEAN DEFAULT FALSE,
            length INTEGER,
            FOREIGN KEY (sample_id) REFERENCES sra_metadata(sample_id)
        );
        
        -- Annotations table
        CREATE TABLE annotations (
            seqhash_id VARCHAR,
            seed_ortholog VARCHAR,
            evalue DOUBLE,
            score DOUBLE,
            eggnog_ogs VARCHAR,
            max_annot_lvl VARCHAR,
            cog_category VARCHAR,
            description VARCHAR,
            preferred_name VARCHAR,
            sample_id VARCHAR,
            PRIMARY KEY (seqhash_id, sample_id),
            FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id),
            FOREIGN KEY (sample_id) REFERENCES sra_metadata(sample_id)
        );
        
        -- GO terms table
        CREATE TABLE go_terms (
            seqhash_id VARCHAR,
            go_term VARCHAR,
            PRIMARY KEY (seqhash_id, go_term),
            FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id)
        );
        
        -- EC numbers table
        CREATE TABLE ec_numbers (
            seqhash_id VARCHAR,
            ec_number VARCHAR,
            PRIMARY KEY (seqhash_id, ec_number),
            FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id)
        );
        
        -- KEGG info table
        CREATE TABLE kegg_info (
            seqhash_id VARCHAR,
            kegg_id VARCHAR,
            kegg_pathway VARCHAR,
            PRIMARY KEY (seqhash_id, kegg_id, kegg_pathway),
            FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id)
        );
        
        -- Clusters table
        CREATE TABLE clusters (
            cluster_id VARCHAR PRIMARY KEY,
            representative_seqhash_id VARCHAR,
            size INTEGER,
            FOREIGN KEY (representative_seqhash_id) REFERENCES sequences(seqhash_id)
        );
        
        -- Cluster members table
        CREATE TABLE cluster_members (
            seqhash_id VARCHAR,
            cluster_id VARCHAR,
            PRIMARY KEY (seqhash_id, cluster_id),
            FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id),
            FOREIGN KEY (cluster_id) REFERENCES clusters(cluster_id)
        );
        
        -- Expression table
        CREATE TABLE expression (
            seqhash_id VARCHAR NOT NULL,
            sample_id VARCHAR NOT NULL,
            tpm DOUBLE NOT NULL,
            num_reads DOUBLE NOT NULL,
            effective_length DOUBLE NOT NULL,
            PRIMARY KEY (seqhash_id, sample_id),
            FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id),
            FOREIGN KEY (sample_id) REFERENCES sra_metadata(sample_id)
        );
        """
        
        self.con.execute(schema_sql)
        
        # Insert sample data
        self.con.execute("""
        INSERT INTO sra_metadata (
            sample_id, organism, study_title, study_abstract
        ) VALUES (
            'SRR12068547', 'Mesoplasma florum', 'Test Study', 'Abstract'
        );
        """)
        
        self.con.execute("""
        INSERT INTO sequences VALUES
        ('v1_DLS_1326316412ebf3de1b3287ad9d63156b914d59fac8091a0ace01b5460d43e49c', 'ACTG', 'SRR12068547', CURRENT_TIMESTAMP, FALSE, 5415),
        ('v1_DLS_bdfe3e5075584cd087ebd251e8ccdb41d07f712fef6076743611cc829c909ace', 'GGCC', 'SRR12068547', CURRENT_TIMESTAMP, FALSE, 2485),
        ('v1_DLS_6f6388f460dca325a7c0c54982c1f0a2e701b872afcf81d06fb166db02ef64cf', 'AATT', 'SRR12068547', CURRENT_TIMESTAMP, TRUE, 2114);
        """)
        
        self.con.execute("""
        INSERT INTO annotations (
            seqhash_id, seed_ortholog, evalue, score, description, preferred_name, sample_id
        ) VALUES 
        ('v1_DLS_1326316412ebf3de1b3287ad9d63156b914d59fac8091a0ace01b5460d43e49c', 'test_ortholog', 0.001, 100, 'Phosphoglycerate dehydrogenase', 'serA', 'SRR12068547'),
        ('v1_DLS_bdfe3e5075584cd087ebd251e8ccdb41d07f712fef6076743611cc829c909ace', 'test_ortholog2', 0.002, 90, 'Preflagellin peptidase', 'PilD', 'SRR12068547');
        """)
        
        # Insert expression data
        self.con.execute("""
        INSERT INTO expression (
            seqhash_id, sample_id, tpm, num_reads, effective_length
        ) VALUES 
        ('v1_DLS_1326316412ebf3de1b3287ad9d63156b914d59fac8091a0ace01b5460d43e49c', 'SRR12068547', 3.762602, 1670.012, 5125.11),
        ('v1_DLS_bdfe3e5075584cd087ebd251e8ccdb41d07f712fef6076743611cc829c909ace', 'SRR12068547', 614.765489, 116867.548, 2195.11),
        ('v1_DLS_6f6388f460dca325a7c0c54982c1f0a2e701b872afcf81d06fb166db02ef64cf', 'SRR12068547', 654.468328, 103387.45, 1824.11);
        """)
    
    def test_expression_summary(self):
        """Test the expression summary query."""
        result = self.query_manager.sequences.get_expression_summary()
        
        # Check for expected columns
        self.assertIn('sample_id', result.columns)
        self.assertIn('mean_tpm', result.columns)
        self.assertIn('max_tpm', result.columns)
        
        # Check the values
        self.assertEqual(result['sample_id'].iloc[0], 'all')
        self.assertEqual(result['expression_records'].iloc[0], 3)
        self.assertGreater(result['max_tpm'].iloc[0], 600)
        
        # Test with sample filter
        result_sample = self.query_manager.sequences.get_expression_summary(self.sample_id)
        self.assertEqual(result_sample['sample_id'].iloc[0], self.sample_id)
    
    def test_expression_distribution(self):
        """Test the expression distribution query."""
        result = self.query_manager.sequences.get_expression_distribution()
        
        # We should have categories in the result
        self.assertIn('expression_level', result.columns)
        self.assertIn('count', result.columns)
        self.assertIn('percentage', result.columns)
        
        # With our test data, we should have High and Very High categories
        levels = result['expression_level'].tolist()
        self.assertIn('Very High (>100 TPM)', levels)
        
        # Test with sample filter
        result_sample = self.query_manager.sequences.get_expression_distribution(self.sample_id)
        self.assertEqual(len(result_sample), len(result))
    
    def test_top_expressed_sequences(self):
        """Test the top expressed sequences query."""
        result = self.query_manager.sequences.get_top_expressed_sequences(limit=2)
        
        # Check the columns
        self.assertIn('seqhash_id', result.columns)
        self.assertIn('tpm', result.columns)
        
        # We should have at most 2 rows due to the limit
        self.assertLessEqual(len(result), 2)
        
        # The rows should be ordered by TPM descending
        tpm_values = result['tpm'].tolist()
        self.assertEqual(tpm_values, sorted(tpm_values, reverse=True))
        
        # The top sequence should be the one with highest TPM
        top_seqhash = result['seqhash_id'].iloc[0]
        self.assertEqual(top_seqhash, 'v1_DLS_6f6388f460dca325a7c0c54982c1f0a2e701b872afcf81d06fb166db02ef64cf')
        
        # Test with sample filter
        result_sample = self.query_manager.sequences.get_top_expressed_sequences(limit=10, sample_id=self.sample_id)
        self.assertEqual(len(result_sample), 3)  # We have 3 total sequences
    
    def test_sequence_expression(self):
        """Test the get expression for sequence query."""
        # Get expression for a specific sequence
        seqhash = 'v1_DLS_bdfe3e5075584cd087ebd251e8ccdb41d07f712fef6076743611cc829c909ace'
        result = self.query_manager.sequences.get_expression_for_sequence(seqhash)
        
        # Should have the expression data
        self.assertIn('tpm', result.columns)
        self.assertIn('sample_id', result.columns)
        
        # Should have the specific TPM value we inserted
        self.assertAlmostEqual(result['tpm'].iloc[0], 614.765489)
        
        # Should have the metadata joined
        self.assertEqual(result['organism'].iloc[0], 'Mesoplasma florum')


if __name__ == '__main__':
    unittest.main()