#!/usr/bin/env python3
"""
End-to-end tests for the Snakemake workflow using Mesoplasma test data.
These tests are marked as slow and will simulate running parts of the workflow.
"""
import os
import tempfile
import shutil
import json
import unittest
from pathlib import Path
import yaml
import pytest
import duckdb
from unittest.mock import patch, MagicMock, mock_open

# Import needed code from planter
from planter.database.builder import SequenceDBBuilder
from planter.database.utils.duckdb_utils import merge_duckdbs, update_duckdb_with_cluster_info
from planter.database.schema.schema_version import get_db_schema_version, ensure_compatibility


@pytest.mark.slow
class TestSnakemakeWorkflow(unittest.TestCase):
    """Test the Snakemake workflow end-to-end with test data."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test environment once for all tests."""
        # Create a temporary directory structure for the workflow
        cls.root_temp_dir = tempfile.mkdtemp()
        
        # Define our test sample ID
        cls.sample_id = "MESOPLASMA"
        
        # Create output directory structure
        cls.output_dir = Path(cls.root_temp_dir) / "output"
        cls.output_dir.mkdir(exist_ok=True)
        cls.sample_dir = cls.output_dir / cls.sample_id
        cls.sample_dir.mkdir(exist_ok=True)
        
        # Create subdirectories
        (cls.sample_dir / "transdecoder").mkdir(exist_ok=True)
        (cls.sample_dir / "eggnog").mkdir(exist_ok=True)
        (cls.sample_dir / "quants").mkdir(exist_ok=True)
        (cls.sample_dir / "logs").mkdir(exist_ok=True)
        
        # Create a minimal test configuration file
        cls.config_path = Path(cls.root_temp_dir) / "config.yaml"
        config = {
            "outdir": str(cls.output_dir),
            "samples": [cls.sample_id],
            "s3_bucket": "test-bucket"
        }
        
        with open(cls.config_path, 'w') as f:
            yaml.dump(config, f)
        
        # Create test data files
        cls._create_test_data()
    
    @classmethod
    def tearDownClass(cls):
        """Clean up after all tests."""
        shutil.rmtree(cls.root_temp_dir)
    
    @classmethod
    def _create_test_data(cls):
        """Create test data files for the workflow."""
        # 1. Create a test protein file (TransDecoder output)
        transdecoder_file = cls.sample_dir / "transdecoder" / f"{cls.sample_id}.pep"
        with open(transdecoder_file, 'w') as f:
            f.write(">MESOPLASMA_seq1.p1 Type=Internal ORF=1 len=100\n")
            f.write("MEPKSLYVGDLHGAFYSIWKLVKSENDNRYLLLVDLNEKIAEMIFNLHGKLDVLSQLPQKK\n")
            f.write(">MESOPLASMA_seq2.p1 Type=Complete ORF=2 len=100\n")
            f.write("MARNVLNAIDVLSRLETHLNGLIRLAVDKMDLSEVITSLPKMRRSFSENLNQLTKRVQEL\n")
            f.write(">MESOPLASMA_seq3.p1 Type=5Prime ORF=3 len=100\n")
            f.write("MVDLPFGWKASKALLRERIKAAQAEWEASARAGPSSPKGVSAPLSPLDITPIAGSVCPVT\n")
        
        # 2. Create a test EggNOG annotations file
        eggnog_file = cls.sample_dir / "eggnog" / f"{cls.sample_id}.emapper.annotations"
        with open(eggnog_file, 'w') as f:
            f.write("# eggNOG-mapper - orthology assignment, gene annotation and functional interpretation\n")
            f.write("# emapper version: 2.1.6\n")
            f.write("# Annotation Report generated on: 2025-03-05 19:49:57\n")
            f.write("#query\tseeds_ortholog\tevalue\tscore\ttaxonomic_range\tOGs\tCOG\tdescription\tpreferred_name\tGOs\tEC\tKEGG\tTC\tCAZy\tBiGG\ttax_scope\teggNOG OGs\tbestOG\tCOG cat\tEggnog score\n")
            f.write("MESOPLASMA_seq1.p1\tNOV1L7TZ_05168\t1e-10\t100.0\t-\t-\t-\tHypothetical protein\t-\tGO:0008150,GO:0003674,GO:0005575\tec:1.1.1.1\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n")
            f.write("MESOPLASMA_seq2.p1\tNOV17X89_08773\t1e-5\t80.0\t-\t-\t-\tAnother protein\t-\tGO:0016491\tec:2.2.2.2\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n")
            f.write("MESOPLASMA_seq3.p1\tNOV15KO2_03347\t1e-8\t90.0\t-\t-\t-\tThird protein\t-\tGO:0005524\tec:3.3.3.3\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n")
        
        # 3. Create a test Quant.json file
        quant_file = cls.sample_dir / "quants" / f"{cls.sample_id}.quant.json"
        quant_data = [
            {
                "Name": "MESOPLASMA_seq1",
                "Length": 60,
                "EffectiveLength": 58.12,
                "TPM": 10.5,
                "NumReads": 100.0,
                "sample": cls.sample_id
            },
            {
                "Name": "MESOPLASMA_seq2",
                "Length": 60,
                "EffectiveLength": 58.12,
                "TPM": 5.2,
                "NumReads": 50.0,
                "sample": cls.sample_id
            },
            {
                "Name": "MESOPLASMA_seq3",
                "Length": 60,
                "EffectiveLength": 58.12,
                "TPM": 20.8,
                "NumReads": 200.0,
                "sample": cls.sample_id
            }
        ]
        with open(quant_file, 'w') as f:
            json.dump(quant_data, f)
    
    def test_01_create_duckdb_step(self):
        """Test the DuckDB creation step of the workflow."""
        # Create path for the sample DuckDB
        duckdb_path = self.sample_dir / f"{self.sample_id}.duckdb"
        
        # Mock SRA info for the builder
        with patch('planter.database.builder.get_sra_info') as mock_get_sra_info:
            mock_get_sra_info.return_value = {
                'organism': 'Mesoplasma florum',
                'study_title': 'Test Mesoplasma Study',
                'study_abstract': 'Test study for workflow testing',
                'bioproject': 'TESTPROJ',
                'biosample': 'TESTSAMPLE',
                'library_strategy': 'RNA-Seq',
                'library_source': 'TRANSCRIPTOMIC',
                'library_selection': 'cDNA',
                'library_layout': 'PAIRED',
                'instrument': 'ILLUMINA',
                'run_spots': '1000000',
                'run_bases': '100000000',
                'run_published': '2025-03-01'
            }
            
            # Run the DuckDB creation (similar to create_duckdb rule)
            with SequenceDBBuilder(str(duckdb_path), output_dir=str(self.output_dir)) as builder:
                results = builder.build_database([self.sample_id])
                
                # Get database summary
                summary = builder.get_database_summary()
            
            # Verify DuckDB was created
            self.assertTrue(duckdb_path.exists(), "DuckDB file should exist")
            
            # Connect to database and verify content
            con = duckdb.connect(str(duckdb_path))
            
            # Check tables
            tables = con.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
            table_names = [t[0] for t in tables]
            
            expected_tables = [
                'schema_version',
                'sra_metadata',
                'sequences',
                'annotations',
                'go_terms',
                'ec_numbers',
                'expression'
            ]
            
            for table in expected_tables:
                self.assertIn(table, table_names, f"Table {table} should exist")
            
            # Verify sequences were loaded
            seq_count = con.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
            self.assertEqual(seq_count, 3, "Should have 3 sequences")
            
            # Verify annotations were loaded
            annot_count = con.execute("SELECT COUNT(*) FROM annotations").fetchone()[0]
            self.assertEqual(annot_count, 3, "Should have 3 annotations")
            
            # Verify GO terms were loaded
            go_count = con.execute("SELECT COUNT(*) FROM go_terms").fetchone()[0]
            self.assertEqual(go_count, 5, "Should have 5 GO terms")
            
            # Verify EC numbers were loaded
            ec_count = con.execute("SELECT COUNT(*) FROM ec_numbers").fetchone()[0]
            self.assertEqual(ec_count, 3, "Should have 3 EC numbers")
            
            # Verify expression data was loaded
            expr_count = con.execute("SELECT COUNT(*) FROM expression").fetchone()[0]
            self.assertEqual(expr_count, 3, "Should have 3 expression records")
            
            # Verify gene-protein mappings were created
            gene_protein_count = con.execute("SELECT COUNT(*) FROM gene_protein_map").fetchone()[0]
            self.assertEqual(gene_protein_count, 3, "Should have 3 gene-protein mappings")
            
            # Close connection
            con.close()
    
    def test_02_database_merge_step(self):
        """Test the database merge step of the workflow."""
        # Create a master database
        master_db_path = self.output_dir / "test_master.duckdb"
        
        # Create a sample database
        sample_db_path = self.sample_dir / f"{self.sample_id}.duckdb"
        
        # Ensure the sample database exists (reuse from previous test or create if needed)
        if not sample_db_path.exists():
            self.test_01_create_duckdb_step()
        
        # Create a minimal master database
        with SequenceDBBuilder(str(master_db_path), output_dir=str(self.output_dir)) as builder:
            # Initialize empty database
            builder.init_database()
            
            # Add a few sequences with different IDs
            with builder.transaction():
                builder.con.execute("""
                    INSERT INTO sequences VALUES
                    ('MASTER_seq1.p1', 'MSTGKV', 'MASTER', CURRENT_TIMESTAMP, TRUE, 'MASTER_seq1.p1', 6),
                    ('MASTER_seq2.p1', 'MVLKSD', 'MASTER', CURRENT_TIMESTAMP, FALSE, 'MASTER_seq1.p1', 6);
                """)
                
                builder.con.execute("""
                    INSERT INTO gene_protein_map VALUES
                    ('MASTER_seq1', 'MASTER_seq1.p1'),
                    ('MASTER_seq2', 'MASTER_seq2.p1');
                """)
                
                builder.con.execute("""
                    INSERT INTO sra_metadata (
                        sample_id, organism, study_title
                    ) VALUES (
                        'MASTER', 'Master Test Organism', 'Master Test Study'
                    );
                """)
        
        # Get initial sequence count in master
        con = duckdb.connect(str(master_db_path))
        initial_seq_count = con.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
        con.close()
        
        # Get initial sequence count in sample
        con = duckdb.connect(str(sample_db_path))
        sample_seq_count = con.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
        con.close()
        
        # Create a merged database path
        merged_db_path = self.output_dir / "merged.duckdb"
        shutil.copy(str(master_db_path), str(merged_db_path))
        
        # Merge the databases
        merge_duckdbs(
            duckdb_paths=[sample_db_path],
            master_db_path=merged_db_path,
            upgrade_schema=True
        )
        
        # Verify the merge
        con = duckdb.connect(str(merged_db_path))
        
        # Check sequence count
        merged_seq_count = con.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
        self.assertEqual(merged_seq_count, initial_seq_count + sample_seq_count, 
                        f"Merged database should have {initial_seq_count + sample_seq_count} sequences")
        
        # Check sample count
        sample_count = con.execute("SELECT COUNT(DISTINCT sample_id) FROM sra_metadata").fetchone()[0]
        self.assertEqual(sample_count, 2, "Merged database should have 2 samples")
        
        # Verify expression data was merged
        expr_count = con.execute("SELECT COUNT(*) FROM expression").fetchone()[0]
        self.assertGreaterEqual(expr_count, 3, "Should have at least 3 expression records")
        
        # Check schema version
        schema_version = get_db_schema_version(merged_db_path)
        self.assertTrue(schema_version >= 2, f"Schema version should be at least 2, got {schema_version}")
        
        con.close()
    
    def test_03_cluster_update_step(self):
        """Test the cluster update step of the workflow."""
        # Create a master database
        master_db_path = self.output_dir / "cluster_test.duckdb"
        
        # Create a basic database with some sequences
        with SequenceDBBuilder(str(master_db_path), output_dir=str(self.output_dir)) as builder:
            # Initialize empty database
            builder.init_database()
            
            # Add some sequences
            with builder.transaction():
                builder.con.execute("""
                    INSERT INTO sequences VALUES
                    ('MESOPLASMA_seq1.p1', 'MEPKSL', 'MESOPLASMA', CURRENT_TIMESTAMP, TRUE, 'MESOPLASMA_seq1.p1', 6),
                    ('MESOPLASMA_seq2.p1', 'MARNVL', 'MESOPLASMA', CURRENT_TIMESTAMP, TRUE, 'MESOPLASMA_seq2.p1', 6),
                    ('MESOPLASMA_seq3.p1', 'MVDLPF', 'MESOPLASMA', CURRENT_TIMESTAMP, TRUE, 'MESOPLASMA_seq3.p1', 6),
                    ('MASTER_seq1.p1', 'MEPKSL', 'MASTER', CURRENT_TIMESTAMP, TRUE, 'MASTER_seq1.p1', 6);
                """)
        
        # Create a clusters.tsv file (similar to MMSeqs2 output)
        clusters_file = self.output_dir / "newClusterDB.tsv"
        with open(clusters_file, 'w') as f:
            # Format: representative_id<tab>member_id
            # Group MESOPLASMA_seq1 and MASTER_seq1 in one cluster
            f.write("MESOPLASMA_seq1.p1\tMESOPLASMA_seq1.p1\n")
            f.write("MESOPLASMA_seq1.p1\tMASTER_seq1.p1\n")
            # MESOPLASMA_seq2 in its own cluster
            f.write("MESOPLASMA_seq2.p1\tMESOPLASMA_seq2.p1\n")
            # MESOPLASMA_seq3 in its own cluster
            f.write("MESOPLASMA_seq3.p1\tMESOPLASMA_seq3.p1\n")
        
        # Update the database with cluster info
        update_duckdb_with_cluster_info(
            master_db_path,
            clusters_file,
            upgrade_schema=True
        )
        
        # Verify the clusters were created
        con = duckdb.connect(str(master_db_path))
        
        # Check cluster count
        cluster_count = con.execute("SELECT COUNT(*) FROM clusters").fetchone()[0]
        self.assertEqual(cluster_count, 3, "Should have 3 clusters")
        
        # Check cluster members
        member_count = con.execute("SELECT COUNT(*) FROM cluster_members").fetchone()[0]
        self.assertEqual(member_count, 4, "Should have 4 cluster members")
        
        # Check that MASTER_seq1 now has MESOPLASMA_seq1 as its repseq_id
        master_seq = con.execute("""
            SELECT repseq_id, is_representative 
            FROM sequences 
            WHERE seqhash_id = 'MASTER_seq1.p1'
        """).fetchone()
        
        self.assertEqual(master_seq[0], "MESOPLASMA_seq1.p1", 
                        "MASTER_seq1.p1 should have MESOPLASMA_seq1.p1 as its repseq_id")
        self.assertEqual(master_seq[1], False, 
                        "MASTER_seq1.p1 should not be marked as representative")
        
        # Verify that MESOPLASMA_seq1 is marked as representative
        mesoplasma_seq = con.execute("""
            SELECT repseq_id, is_representative 
            FROM sequences 
            WHERE seqhash_id = 'MESOPLASMA_seq1.p1'
        """).fetchone()
        
        self.assertEqual(mesoplasma_seq[0], "MESOPLASMA_seq1.p1", 
                        "MESOPLASMA_seq1.p1 should have itself as repseq_id")
        self.assertEqual(mesoplasma_seq[1], True, 
                        "MESOPLASMA_seq1.p1 should be marked as representative")
        
        con.close()
    
    def test_04_schema_compatibility(self):
        """Test schema compatibility handling in the workflow."""
        # Create a v1 schema database
        v1_db_path = self.output_dir / "v1_schema.duckdb"
        
        # Create v1 schema manually
        con = duckdb.connect(str(v1_db_path))
        con.execute("""
            CREATE TABLE sequences (
                seqhash_id VARCHAR PRIMARY KEY,
                sequence VARCHAR,
                sample_id VARCHAR,
                assembly_date TIMESTAMP,
                repseq_id VARCHAR,
                length INTEGER
            )
        """)
        
        con.execute("""
            INSERT INTO sequences VALUES
            ('v1_seq1.p1', 'MSTGKV', 'v1_sample', CURRENT_TIMESTAMP, 'v1_seq1.p1', 6);
        """)
        con.close()
        
        # Verify it's detected as v1
        schema_version = get_db_schema_version(v1_db_path)
        self.assertEqual(schema_version, 1, f"Should detect schema version 1, got {schema_version}")
        
        # Use ensure_compatibility to upgrade to v2
        new_version, was_upgraded = ensure_compatibility(v1_db_path, required_version=2)
        self.assertEqual(new_version, 2, "Should upgrade to version 2")
        self.assertTrue(was_upgraded, "Should report that upgrade was performed")
        
        # Verify the upgrade
        con = duckdb.connect(str(v1_db_path))
        
        # Check schema_version table exists
        has_version_table = con.execute("""
            SELECT COUNT(*) FROM sqlite_master 
            WHERE type='table' AND name='schema_version'
        """).fetchone()[0]
        
        self.assertTrue(has_version_table, "schema_version table should exist")
        
        # Check is_representative column was added
        has_rep_column = con.execute("""
            SELECT COUNT(*) FROM pragma_table_info('sequences') 
            WHERE name = 'is_representative'
        """).fetchone()[0]
        
        self.assertTrue(has_rep_column, "is_representative column should exist")
        
        # Check schema version in table
        db_version = con.execute("SELECT MAX(version) FROM schema_version").fetchone()[0]
        self.assertEqual(db_version, 2, "Recorded schema version should be 2")
        
        con.close()
    
    @patch('planter.database.builder.get_sra_info')
    @patch('builtins.open', new_callable=mock_open)
    @patch('subprocess.run')
    def test_05_full_workflow_simulation(self, mock_subprocess, mock_file, mock_get_sra_info):
        """Simulate full workflow execution with mocks."""
        # This test simulates the full workflow by mocking external calls
        
        # Mock SRA info response
        mock_get_sra_info.return_value = {
            'organism': 'Mesoplasma florum',
            'study_title': 'Full Workflow Test',
            'study_abstract': 'Simulating the full workflow',
            'bioproject': 'TESTPROJ',
            'biosample': 'TESTSAMPLE',
            'library_strategy': 'RNA-Seq',
            'library_source': 'TRANSCRIPTOMIC',
            'library_selection': 'cDNA',
            'library_layout': 'PAIRED',
            'instrument': 'ILLUMINA',
            'run_spots': '1000000',
            'run_bases': '100000000',
            'run_published': '2025-03-01'
        }
        
        # Set up paths
        output_dir = self.output_dir
        sample_id = self.sample_id
        workflow_dir = Path(self.root_temp_dir) / "workflow"
        workflow_dir.mkdir(exist_ok=True)
        
        # Create a Snakefile for testing
        snakefile_path = workflow_dir / "Snakefile"
        
        # Create output paths that would be generated by workflow
        sample_db_path = output_dir / sample_id / f"{sample_id}.duckdb"
        master_db_path = output_dir / "updated_master.duckdb"
        
        # Ensure the sample database exists from previous tests
        if not sample_db_path.exists():
            self.test_01_create_duckdb_step()
        
        # Create a mock master database
        if not master_db_path.exists():
            with SequenceDBBuilder(str(master_db_path), output_dir=str(output_dir)) as builder:
                builder.init_database()
        
        # Simulate the update_database rule steps
        
        # 1. Create log directory
        log_dir = output_dir / "logs"
        log_dir.mkdir(exist_ok=True)
        
        # 2. Ensure schema compatibility
        schema_version, was_upgraded = ensure_compatibility(
            master_db_path,
            required_version=None  # Use latest
        )
        
        # 3. Create temp directory
        temp_dir = output_dir / "tmp"
        temp_dir.mkdir(exist_ok=True)
        
        # 4. Create cluster file
        cluster_file = temp_dir / "newClusterDB.tsv"
        with open(cluster_file, 'w') as f:
            f.write("MESOPLASMA_seq1.p1\tMESOPLASMA_seq1.p1\n")
            f.write("MESOPLASMA_seq2.p1\tMESOPLASMA_seq2.p1\n")
            f.write("MESOPLASMA_seq3.p1\tMESOPLASMA_seq3.p1\n")
        
        # 5. Merge databases
        merge_duckdbs(
            duckdb_paths=[sample_db_path],
            master_db_path=master_db_path,
            upgrade_schema=True
        )
        
        # 6. Update with cluster info
        update_duckdb_with_cluster_info(
            master_db_path,
            cluster_file,
            upgrade_schema=True
        )
        
        # 7. Verify final database
        con = duckdb.connect(str(master_db_path))
        
        # Get schema version
        final_schema_version = get_db_schema_version(master_db_path)
        self.assertTrue(final_schema_version >= 2, 
                        f"Final schema version should be at least 2, got {final_schema_version}")
        
        # Check that data was properly merged
        sequence_count = con.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
        self.assertGreaterEqual(sequence_count, 3, 
                                f"Should have at least 3 sequences, got {sequence_count}")
        
        # Check expression data
        expr_count = con.execute("SELECT COUNT(*) FROM expression").fetchone()[0]
        self.assertGreaterEqual(expr_count, 3, 
                                f"Should have at least 3 expression records, got {expr_count}")
        
        # Check cluster data
        cluster_count = con.execute("SELECT COUNT(*) FROM clusters").fetchone()[0]
        self.assertGreaterEqual(cluster_count, 3, 
                                f"Should have at least 3 clusters, got {cluster_count}")
        
        con.close()
        
        # Mark test as successful - simulating a complete workflow run
        self.assertTrue(True, "Full workflow simulation completed successfully")


if __name__ == "__main__":
    # To run individual tests from command line:
    # python -m unittest tests.workflow.test_snakemake_workflow.TestSnakemakeWorkflow.test_01_create_duckdb_step
    unittest.main()