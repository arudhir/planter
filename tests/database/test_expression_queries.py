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
        # Create tables one by one to avoid foreign key issues
        
        # SRA metadata table
        self.con.execute("""
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
        """)
        
        # Sequences table (protein sequences)
        self.con.execute("""
        CREATE TABLE sequences (
            seqhash_id VARCHAR PRIMARY KEY,
            sequence VARCHAR NOT NULL,
            sample_id VARCHAR NOT NULL,
            assembly_date TIMESTAMP,
            is_representative BOOLEAN DEFAULT FALSE,
            length INTEGER,
            FOREIGN KEY (sample_id) REFERENCES sra_metadata(sample_id)
        );
        """)
        
        # Gene-Protein mapping table
        self.con.execute("""
        CREATE TABLE gene_protein_map (
            gene_seqhash_id VARCHAR PRIMARY KEY,
            protein_seqhash_id VARCHAR NOT NULL,
            FOREIGN KEY (protein_seqhash_id) REFERENCES sequences(seqhash_id)
        );
        """)
        
        # Annotations table
        self.con.execute("""
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
        """)
        
        # GO terms table
        self.con.execute("""
        CREATE TABLE go_terms (
            seqhash_id VARCHAR,
            go_term VARCHAR,
            PRIMARY KEY (seqhash_id, go_term),
            FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id)
        );
        """)
        
        # EC numbers table
        self.con.execute("""
        CREATE TABLE ec_numbers (
            seqhash_id VARCHAR,
            ec_number VARCHAR,
            PRIMARY KEY (seqhash_id, ec_number),
            FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id)
        );
        """)
        
        # KEGG info table
        self.con.execute("""
        CREATE TABLE kegg_info (
            seqhash_id VARCHAR,
            kegg_id VARCHAR,
            kegg_pathway VARCHAR,
            PRIMARY KEY (seqhash_id, kegg_id, kegg_pathway),
            FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id)
        );
        """)
        
        # Clusters table
        self.con.execute("""
        CREATE TABLE clusters (
            cluster_id VARCHAR PRIMARY KEY,
            representative_seqhash_id VARCHAR,
            size INTEGER,
            FOREIGN KEY (representative_seqhash_id) REFERENCES sequences(seqhash_id)
        );
        """)
        
        # Cluster members table
        self.con.execute("""
        CREATE TABLE cluster_members (
            seqhash_id VARCHAR,
            cluster_id VARCHAR,
            PRIMARY KEY (seqhash_id, cluster_id),
            FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id),
            FOREIGN KEY (cluster_id) REFERENCES clusters(cluster_id)
        );
        """)
        
        # Expression table
        self.con.execute("""
        CREATE TABLE expression (
            gene_seqhash_id VARCHAR NOT NULL,
            sample_id VARCHAR NOT NULL,
            tpm DOUBLE NOT NULL,
            num_reads DOUBLE NOT NULL,
            effective_length DOUBLE NOT NULL,
            PRIMARY KEY (gene_seqhash_id, sample_id),
            FOREIGN KEY (gene_seqhash_id) REFERENCES gene_protein_map(gene_seqhash_id),
            FOREIGN KEY (sample_id) REFERENCES sra_metadata(sample_id)
        );
        """)
        
        # Insert sample data
        self.con.execute("""
        INSERT INTO sra_metadata (
            sample_id, organism, study_title, study_abstract
        ) VALUES (
            'SRR12068547', 'Mesoplasma florum', 'Test Study', 'Abstract'
        );
        """)
        
        # Insert protein sequences
        self.con.execute("""
        INSERT INTO sequences VALUES
        ('v1_DLS_1326316412ebf3de1b3287ad9d63156b914d59fac8091a0ace01b5460d43e49c.p1', 'ACTG', 'SRR12068547', CURRENT_TIMESTAMP, FALSE, 5415),
        ('v1_DLS_bdfe3e5075584cd087ebd251e8ccdb41d07f712fef6076743611cc829c909ace.p1', 'GGCC', 'SRR12068547', CURRENT_TIMESTAMP, FALSE, 2485),
        ('v1_DLS_6f6388f460dca325a7c0c54982c1f0a2e701b872afcf81d06fb166db02ef64cf.p1', 'AATT', 'SRR12068547', CURRENT_TIMESTAMP, TRUE, 2114);
        """)
        
        # Insert gene-protein mappings
        self.con.execute("""
        INSERT INTO gene_protein_map (
            gene_seqhash_id, protein_seqhash_id
        ) VALUES 
        ('v1_DLS_1326316412ebf3de1b3287ad9d63156b914d59fac8091a0ace01b5460d43e49c', 'v1_DLS_1326316412ebf3de1b3287ad9d63156b914d59fac8091a0ace01b5460d43e49c.p1'),
        ('v1_DLS_bdfe3e5075584cd087ebd251e8ccdb41d07f712fef6076743611cc829c909ace', 'v1_DLS_bdfe3e5075584cd087ebd251e8ccdb41d07f712fef6076743611cc829c909ace.p1'),
        ('v1_DLS_6f6388f460dca325a7c0c54982c1f0a2e701b872afcf81d06fb166db02ef64cf', 'v1_DLS_6f6388f460dca325a7c0c54982c1f0a2e701b872afcf81d06fb166db02ef64cf.p1');
        """)
        
        # Insert annotations (for protein sequences)
        self.con.execute("""
        INSERT INTO annotations (
            seqhash_id, seed_ortholog, evalue, score, description, preferred_name, sample_id
        ) VALUES 
        ('v1_DLS_1326316412ebf3de1b3287ad9d63156b914d59fac8091a0ace01b5460d43e49c.p1', 'test_ortholog', 0.001, 100, 'Phosphoglycerate dehydrogenase', 'serA', 'SRR12068547'),
        ('v1_DLS_bdfe3e5075584cd087ebd251e8ccdb41d07f712fef6076743611cc829c909ace.p1', 'test_ortholog2', 0.002, 90, 'Preflagellin peptidase', 'PilD', 'SRR12068547');
        """)
        
        # Insert expression data (for gene sequences)
        self.con.execute("""
        INSERT INTO expression (
            gene_seqhash_id, sample_id, tpm, num_reads, effective_length
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
        
        # Print the result for verification
        print(f"\nExpression Summary:\n{result}")
        
        # Test with sample filter
        result_sample = self.query_manager.sequences.get_expression_summary(self.sample_id)
        self.assertEqual(result_sample['sample_id'].iloc[0], self.sample_id)
        print(f"\nExpression Summary for sample {self.sample_id}:\n{result_sample}")
    
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
        
        # Print the distribution
        print(f"\nExpression Distribution:\n{result}")
        
        # Test with sample filter
        result_sample = self.query_manager.sequences.get_expression_distribution(self.sample_id)
        self.assertEqual(len(result_sample), len(result))
        print(f"\nExpression Distribution for sample {self.sample_id}:\n{result_sample}")
    
    def test_top_expressed_sequences(self):
        """Test the top expressed sequences query."""
        result = self.query_manager.sequences.get_top_expressed_sequences(limit=2)
        
        # Check the columns
        self.assertIn('gene_seqhash_id', result.columns)
        self.assertIn('tpm', result.columns)
        
        # We should have at most 2 rows due to the limit
        self.assertLessEqual(len(result), 2)
        
        # The rows should be ordered by TPM descending
        tpm_values = result['tpm'].tolist()
        self.assertEqual(tpm_values, sorted(tpm_values, reverse=True))
        
        # The top sequence should be the one with highest TPM
        top_seqhash = result['gene_seqhash_id'].iloc[0]
        self.assertEqual(top_seqhash, 'v1_DLS_6f6388f460dca325a7c0c54982c1f0a2e701b872afcf81d06fb166db02ef64cf')
        
        # Print the results
        print(f"\nTop Expressed Sequences (limit=2):\n{result}")
        
        # Test with sample filter
        result_sample = self.query_manager.sequences.get_top_expressed_sequences(limit=10, sample_id=self.sample_id)
        self.assertEqual(len(result_sample), 3)  # We have 3 total sequences
        print(f"\nTop Expressed Sequences for sample {self.sample_id}:\n{result_sample}")
    
    def test_sequence_expression(self):
        """Test the get expression for sequence query."""
        # Get expression for a specific gene sequence
        gene_seqhash = 'v1_DLS_bdfe3e5075584cd087ebd251e8ccdb41d07f712fef6076743611cc829c909ace'
        result = self.query_manager.sequences.get_expression_for_sequence(gene_seqhash)
        
        # Should have the expression data
        self.assertIn('tpm', result.columns)
        self.assertIn('sample_id', result.columns)
        
        # Should have the specific TPM value we inserted
        self.assertAlmostEqual(result['tpm'].iloc[0], 614.765489)
        
        # Should have the metadata joined
        self.assertEqual(result['organism'].iloc[0], 'Mesoplasma florum')
        
        # Print gene sequence expression results
        print(f"\nExpression for gene sequence {gene_seqhash}:\n{result}")
        
        # Now test with protein seqhash - should find expression through gene_protein_map
        protein_seqhash = 'v1_DLS_bdfe3e5075584cd087ebd251e8ccdb41d07f712fef6076743611cc829c909ace.p1'
        result = self.query_manager.sequences.get_expression_for_sequence(protein_seqhash)
        
        # Should still find the expression data via the protein->gene mapping
        self.assertAlmostEqual(result['tpm'].iloc[0], 614.765489)
        
        # Print protein sequence expression results
        print(f"\nExpression for protein sequence {protein_seqhash}:\n{result}")
    
    def test_annotation_to_expression(self):
        """Test getting expression data for annotated sequences."""
        # First, get an annotated sequence
        annotated_sequences = self.con.execute("""
            SELECT a.seqhash_id, a.preferred_name 
            FROM annotations a 
            LIMIT 1
        """).fetchall()
        
        if len(annotated_sequences) == 0:
            self.skipTest("No annotated sequences available")
        
        protein_seqhash = annotated_sequences[0][0]
        protein_name = annotated_sequences[0][1]
        
        # Get gene_seqhash from protein_seqhash
        gene_seqhash = self.con.execute("""
            SELECT gene_seqhash_id 
            FROM gene_protein_map 
            WHERE protein_seqhash_id = ?
        """, [protein_seqhash]).fetchone()[0]
        
        # Now get expression data for this gene
        expression_data = self.con.execute("""
            SELECT e.gene_seqhash_id, e.tpm, e.num_reads
            FROM expression e
            WHERE e.gene_seqhash_id = ?
        """, [gene_seqhash]).fetchall()
        
        self.assertGreater(len(expression_data), 0, "Should find expression data for annotated gene")
        
        # Print the results
        print(f"\nJoin across tables: Annotation to Expression")
        print(f"Protein: {protein_seqhash} ({protein_name})")
        print(f"Gene: {gene_seqhash}")
        print(f"Expression TPM: {expression_data[0][1]}")
        
        # Test a full join query that combines annotations, gene-protein mapping, and expression
        joined_data = self.con.execute("""
            SELECT 
                a.seqhash_id as protein_id,
                a.preferred_name,
                a.description,
                gpm.gene_seqhash_id as gene_id,
                e.tpm,
                e.num_reads
            FROM annotations a
            JOIN gene_protein_map gpm ON a.seqhash_id = gpm.protein_seqhash_id
            JOIN expression e ON gpm.gene_seqhash_id = e.gene_seqhash_id
            ORDER BY e.tpm DESC
        """).fetchall()
        
        self.assertGreater(len(joined_data), 0, "Should find joined data across all tables")
        
        # Print the joined results
        print("\nAnnotation-Expression joined data (all tables):")
        for row in joined_data:
            print(f"Protein: {row[0]} | Gene: {row[3]} | Name: {row[1]} | TPM: {row[4]}")
            
    def test_expression_level_categories(self):
        """Test expression level categories across different genes."""
        # Query that categorizes expression levels
        result = self.con.execute("""
            SELECT
                CASE 
                    WHEN tpm < 1 THEN 'Low (<1 TPM)'
                    WHEN tpm < 10 THEN 'Medium (1-10 TPM)'
                    WHEN tpm < 100 THEN 'High (10-100 TPM)'
                    ELSE 'Very High (>100 TPM)'
                END as expression_level,
                COUNT(*) as gene_count,
                ROUND(AVG(tpm), 2) as avg_tpm
            FROM expression
            GROUP BY expression_level
            ORDER BY avg_tpm DESC
        """).fetchall()
        
        # Print the distribution
        print("\nExpression level distribution with statistics:")
        for row in result:
            print(f"{row[0]}: {row[1]} genes (avg TPM: {row[2]})")
            
        self.assertGreater(len(result), 1, "Should have multiple expression level categories")


if __name__ == '__main__':
    unittest.main()