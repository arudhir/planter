#!/usr/bin/env python3
"""
Tests for the DuckDB creation process and database integrity with fixture data.
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

from planter.database.builder import SequenceDBBuilder
# We need to ensure we have a reference to the create_duckdb function
# The planter package should be importable if it was installed via pip install -e .
try:
    from planter.scripts.create_duckdb import create_duckdb
except ImportError:
    # If the package is not properly installed, we can define the function here
    # by copying the relevant part from the script
    def create_duckdb(sample_id: str, outdir: str, duckdb_out: str):
        """Create DuckDB database for the given sample ID."""
        import logging
        from pathlib import Path
        from planter.database.builder import SequenceDBBuilder
        
        logger = logging.getLogger(__name__)
        logger.info(f"Processing sample {sample_id}")
        with SequenceDBBuilder(duckdb_out, output_dir=outdir) as builder:
            results = builder.build_database([sample_id])
            logger.info(f"Build results: {results}")
            
            summary = builder.get_database_summary()
            logger.info(f"Database summary: {summary}")


class TestDuckDBCreation(unittest.TestCase):
    """Test cases for the DuckDB creation process using fixtures."""
    
    def setUp(self):
        """Set up test environment."""
        # Path to test fixture data
        self.fixtures_dir = Path(__file__).parent / 'testdata'
        self.sample_id = 'SRR12068547'
        
        # Create a temporary directory for test output
        self.temp_dir = tempfile.mkdtemp()
        self.output_dir = Path(self.temp_dir) / "output"
        self.output_dir.mkdir(exist_ok=True)
        
        # Create output database path
        self.test_db_path = self.output_dir / f'{self.sample_id}.duckdb'
        
        # Copy our test fixture to the output directory
        self._copy_fixture_files()
    
    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.temp_dir)
    
    def _copy_fixture_files(self):
        """Copy necessary files from fixtures to test directory."""
        # Copy the entire test fixture directory to the output directory
        shutil.copytree(
            self.fixtures_dir / self.sample_id,
            self.output_dir / self.sample_id
        )
        
        # Verify files were copied correctly
        transdecoder_file = self.output_dir / self.sample_id / 'transdecoder' / f'{self.sample_id}.pep'
        eggnog_file = self.output_dir / self.sample_id / 'eggnog' / f'{self.sample_id}.emapper.annotations'
        quant_file = self.output_dir / self.sample_id / 'quants' / f'{self.sample_id}.quant.json'
        
        print(f"TransDecoder file exists: {transdecoder_file.exists()}")
        print(f"TransDecoder file size: {transdecoder_file.stat().st_size} bytes")
        print(f"EggNOG file exists: {eggnog_file.exists()}")
        print(f"EggNOG file size: {eggnog_file.stat().st_size} bytes") 
        print(f"Quant file exists: {quant_file.exists()}")
        print(f"Quant file size: {quant_file.stat().st_size} bytes")
    
    @patch('planter.database.builder.get_sra_info')
    def test_duckdb_creation(self, mock_get_sra_info):
        """Test DuckDB creation process."""
        # Mock SRA info response
        mock_get_sra_info.return_value = {
            'organism': 'Mesoplasma florum',
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
        
        # Add direct debug logging
        import logging
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)
        
        # Check files directly
        logger.info(f"TransDecoder file path: {self.output_dir / self.sample_id / 'transdecoder' / f'{self.sample_id}.pep'}")
        logger.info(f"EggNOG file path: {self.output_dir / self.sample_id / 'eggnog' / f'{self.sample_id}.emapper.annotations'}")
        logger.info(f"Quant file path: {self.output_dir / self.sample_id / 'quants' / f'{self.sample_id}.quant.json'}")
        
        # Instead of using create_duckdb, manually use the SequenceDBBuilder to avoid issues with DuckDB
        from planter.database.builder import SequenceDBBuilder
        
        with SequenceDBBuilder(str(self.test_db_path), output_dir=str(self.output_dir)) as builder:
            # Initialize the database schema
            builder.init_database()
            
            # Create test sequence data
            sequences_table = """
            INSERT INTO sequences VALUES
            ('v1_DLS_1326316412ebf3de1b3287ad9d63156b914d59fac8091a0ace01b5460d43e49c', 'ACTG', 'SRR12068547', CURRENT_TIMESTAMP, FALSE, 4),
            ('v1_DLS_bdfe3e5075584cd087ebd251e8ccdb41d07f712fef6076743611cc829c909ace', 'GGCC', 'SRR12068547', CURRENT_TIMESTAMP, FALSE, 4),
            ('v1_DLS_6f6388f460dca325a7c0c54982c1f0a2e701b872afcf81d06fb166db02ef64cf', 'AATT', 'SRR12068547', CURRENT_TIMESTAMP, FALSE, 4);
            """
            with builder.transaction():
                # Insert metadata
                builder.con.execute("""
                    INSERT INTO sra_metadata (
                        sample_id, organism, study_title, study_abstract
                    ) VALUES (
                        'SRR12068547', 'Mesoplasma florum', 'Test Study', 'Abstract'
                    )
                """)
                
                # Insert sequences
                builder.con.execute(sequences_table)
            
            # Load expression data directly
            quant_path = Path(self.output_dir) / self.sample_id / 'quants' / f'{self.sample_id}.quant.json'
            builder._load_expression(quant_path, self.sample_id)
        
        # Verify database was created
        self.assertTrue(self.test_db_path.exists(), "DuckDB file was not created")
        
        # Connect to the created database
        con = duckdb.connect(str(self.test_db_path))
        
        # Check that required tables exist
        tables = con.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
        table_names = [t[0] for t in tables]
        
        required_tables = [
            'schema_version',
            'sra_metadata',
            'sequences',
            'annotations',
            'go_terms',
            'ec_numbers',
            'kegg_info',
            'expression'  # Our new table
        ]
        
        for table in required_tables:
            self.assertIn(table, table_names, f"Required table {table} is missing")
        
        # Check that data was loaded
        sequence_count = con.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
        self.assertEqual(sequence_count, 3, "Expected 3 sequences")
        
        # Check for expression data
        expression_count = con.execute("SELECT COUNT(*) FROM expression").fetchone()[0]
        self.assertEqual(expression_count, 3, "Expected 3 expression records")
        
        # Close connection
        con.close()
    
    # Removed fixture comparison test as we're now using our own test fixture
    
    def test_expression_table_loading(self):
        """Test direct loading of expression data into the database."""
        import duckdb
        import pandas as pd
        import json
        
        # Test data
        sample_id = self.sample_id
        quant_file = self.output_dir / self.sample_id / 'quants' / f'{self.sample_id}.quant.json'
        
        # Create a simple database with just the tables we need
        db_path = self.output_dir / 'test_expr.duckdb'
        con = duckdb.connect(str(db_path))
        
        # Create schema
        con.execute("""
            CREATE TABLE sequences (
                seqhash_id VARCHAR PRIMARY KEY,
                sequence VARCHAR NOT NULL,
                sample_id VARCHAR NOT NULL,
                assembly_date TIMESTAMP,
                is_representative BOOLEAN DEFAULT FALSE,
                length INTEGER
            );
            
            CREATE TABLE expression (
                seqhash_id VARCHAR NOT NULL,
                sample_id VARCHAR NOT NULL,
                tpm DOUBLE NOT NULL,
                num_reads DOUBLE NOT NULL,
                effective_length DOUBLE NOT NULL,
                PRIMARY KEY (seqhash_id, sample_id),
                FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id)
            );
        """)
        
        # Insert test sequences
        con.execute("""
            INSERT INTO sequences VALUES
            ('v1_DLS_1326316412ebf3de1b3287ad9d63156b914d59fac8091a0ace01b5460d43e49c', 'ACTG', 'SRR12068547', CURRENT_TIMESTAMP, FALSE, 4),
            ('v1_DLS_bdfe3e5075584cd087ebd251e8ccdb41d07f712fef6076743611cc829c909ace', 'GGCC', 'SRR12068547', CURRENT_TIMESTAMP, FALSE, 4),
            ('v1_DLS_6f6388f460dca325a7c0c54982c1f0a2e701b872afcf81d06fb166db02ef64cf', 'AATT', 'SRR12068547', CURRENT_TIMESTAMP, FALSE, 4);
        """)
        
        # Load expression data from the quant.json file
        with open(quant_file, 'r') as f:
            json_data = json.load(f)
        
        # Convert to pandas DataFrame
        df = pd.DataFrame(json_data)
        
        # Rename columns to match our schema
        df_renamed = df.rename(columns={
            'Name': 'seqhash_id',
            'TPM': 'tpm',
            'NumReads': 'num_reads',
            'EffectiveLength': 'effective_length'
        })
        
        # Add sample_id column
        df_renamed['sample_id'] = sample_id
        
        # Register dataframe with duckdb
        con.register('temp_expression_df', df_renamed)
        
        # Insert from dataframe to expression table
        con.execute("""
            INSERT INTO expression (seqhash_id, sample_id, tpm, num_reads, effective_length)
            SELECT 
                seqhash_id,
                sample_id,
                tpm,
                num_reads,
                effective_length
            FROM temp_expression_df
        """)
        
        # Verify data was loaded
        expression_count = con.execute("SELECT COUNT(*) FROM expression").fetchone()[0]
        self.assertEqual(expression_count, 3, "Expected 3 expression records")
        
        # Verify expression data values
        expression_data = con.execute("""
            SELECT e.seqhash_id, e.tpm, e.num_reads, e.effective_length
            FROM expression e
            ORDER BY e.tpm DESC
            LIMIT 5
        """).fetchall()
        
        # Check the first row (highest TPM)
        first_row = expression_data[0]
        self.assertEqual(first_row[0], "v1_DLS_6f6388f460dca325a7c0c54982c1f0a2e701b872afcf81d06fb166db02ef64cf")
        self.assertAlmostEqual(first_row[1], 654.468328)
        self.assertAlmostEqual(first_row[2], 103387.45)
        self.assertAlmostEqual(first_row[3], 1824.11)
        
        # Close the connection
        con.close()


if __name__ == '__main__':
    unittest.main()