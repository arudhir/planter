import logging
import time
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import duckdb
import pandas as pd
from Bio import SeqIO

from .utils.sra import get_sra_info  # Add this import


@dataclass
class SamplePaths:
    """Paths for sample data files"""

    sequences: Path  # transdecoder .pep file
    annotations: Path  # eggnog annotations file
    expression: Path  # salmon quantification json file


@dataclass
class ProcessingResult:
    """Result of processing a single sample"""

    sample_id: str
    status: str
    sequences_loaded: int = 0
    annotations_loaded: int = 0
    duplicates: int = 0
    error: Optional[str] = None


class SequenceDBBuilder:
    """Builds and populates sequence database from SRA IDs."""

    def __init__(self, db_path: str, output_dir: Path):
        self.db_path = db_path
        self.output_dir = Path(output_dir)
        self.con = None
        self.logger = logging.getLogger(__name__)
        self.schema_dir = Path(__file__).parent / "schema" / "migrations"

    def __enter__(self):
        self.con = duckdb.connect(self.db_path)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.con:
            self.con.close()

    # Keep track of transaction nesting
    _transaction_level = 0

    @contextmanager
    def transaction(self):
        """Context manager for database transactions that supports nesting"""
        is_outermost = self._transaction_level == 0
        self._transaction_level += 1

        try:
            if is_outermost:
                self.con.execute("BEGIN")
            yield
            if is_outermost:
                self.con.execute("COMMIT")
        except Exception as e:
            if is_outermost:
                self.con.execute("ROLLBACK")
            self.logger.error(f"Transaction failed: {str(e)}")
            raise
        finally:
            self._transaction_level -= 1

    def execute_transaction(self, func, *args, **kwargs):
        """Execute a function within a transaction"""
        with self.transaction():
            return func(*args, **kwargs)

    def _init_version_tracking(self):
        """Initialize schema version tracking"""
        version_sql = """
        CREATE TABLE IF NOT EXISTS schema_version (
            version INTEGER PRIMARY KEY,
            migration_name VARCHAR NOT NULL,
            applied_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
        """
        self.con.execute(version_sql)

    def _get_applied_migrations(self) -> List[str]:
        """Get list of already applied migrations"""
        return [
            row[0]
            for row in self.con.execute(
                """
            SELECT migration_name 
            FROM schema_version 
            ORDER BY version
        """
            ).fetchall()
        ]

    def clean_database(self):
        """Drop all existing tables in the correct order to respect foreign key constraints."""
        self.logger.info("Cleaning up any existing tables...")
        
        # We need to drop tables in the correct order based on foreign key dependencies
        # First, drop tables that reference other tables
        cleanup_sql1 = """
        DROP TABLE IF EXISTS cluster_members;
        DROP TABLE IF EXISTS annotations;
        DROP TABLE IF EXISTS go_terms;
        DROP TABLE IF EXISTS ec_numbers;
        DROP TABLE IF EXISTS expression;
        DROP TABLE IF EXISTS kegg_info;
        """
        
        # Then drop tables that are referenced by other tables
        cleanup_sql2 = """
        DROP TABLE IF EXISTS clusters;
        DROP TABLE IF EXISTS gene_protein_map;
        DROP TABLE IF EXISTS sequences;
        DROP TABLE IF EXISTS sra_metadata;
        """
        
        # Finally drop any temporary tables
        cleanup_sql3 = """
        DROP TABLE IF EXISTS schema_version;
        DROP TABLE IF EXISTS temp_annotations;
        DROP TABLE IF EXISTS temp_clusters;
        """
        
        with self.transaction():
            try:
                # Execute each batch of drops separately
                self.con.execute(cleanup_sql1)
                self.con.execute(cleanup_sql2)
                self.con.execute(cleanup_sql3)
            except Exception as e:
                self.logger.error(f"Error during database cleanup: {str(e)}")
                # If we can't drop tables, try to disable foreign key constraints and retry
                try:
                    self.con.execute("PRAGMA foreign_keys=OFF;")
                    self.con.execute(cleanup_sql1)
                    self.con.execute(cleanup_sql2)
                    self.con.execute(cleanup_sql3)
                    self.con.execute("PRAGMA foreign_keys=ON;")
                except Exception as inner_e:
                    self.logger.error(f"Failed to clean database even with foreign_keys off: {str(inner_e)}")
                    # Last resort - try to drop each table individually
                    for table in [
                        "cluster_members", "annotations", "go_terms", "ec_numbers", 
                        "expression", "kegg_info", "clusters", "gene_protein_map", 
                        "sequences", "sra_metadata", "schema_version", "temp_annotations", 
                        "temp_clusters"
                    ]:
                        try:
                            self.con.execute(f"DROP TABLE IF EXISTS {table}")
                        except Exception:
                            pass
    
    def init_database(self):
        """Initialize database with clean schema from migration files."""
        # Attempt to clean up any incomplete migrations first
        try:
            self._cleanup_incomplete_migrations()
        except Exception as e:
            self.logger.error(f"Error cleaning up incomplete migrations: {str(e)}")
        
        # Clean database 
        try:
            self.clean_database()
        except Exception as e:
            self.logger.error(f"Error during database cleanup: {str(e)}")
            # Try more aggressive cleanup
            try:
                self.logger.info("Trying more aggressive cleanup...")
                self.con.close()
                import os
                if os.path.exists(self.db_path):
                    os.remove(self.db_path)
                self.con = duckdb.connect(self.db_path)
                self.logger.info("Recreated database file from scratch")
            except Exception as e2:
                self.logger.error(f"Failed to recreate database: {str(e2)}")

        # First, create schema version table
        try:
            self.con.execute("""
                CREATE TABLE IF NOT EXISTS schema_version (
                    version INTEGER PRIMARY KEY,
                    migration_name VARCHAR NOT NULL,
                    applied_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            """)
        except Exception as e:
            self.logger.error(f"Failed to create schema_version table: {str(e)}")
            return False
            
        # Apply migrations directly with DuckDB-compatible SQL
        self._apply_duckdb_migrations_directly()
        # Add this to your init_database method, right after creating the tables:
        self._fix_gene_protein_schema()

    def _fix_gene_protein_schema(self):
        """
        Fix the gene_protein_map table schema if it has gene_seqhash_id as a PRIMARY KEY.
        This is a one-time fix to update an existing database.
        """
        try:
            # Check if the table exists and what its schema is
            table_info = self.con.execute("""
                PRAGMA table_info(gene_protein_map)
            """).fetchall()
            
            # If the table exists, examine its schema
            if table_info:
                has_primary_key = False
                primary_key_column = None
                
                for col_info in table_info:
                    if col_info[5] > 0:  # pk column in sqlite_master
                        has_primary_key = True
                        primary_key_column = col_info[1]  # column name
                
                if has_primary_key and primary_key_column == "gene_seqhash_id":
                    self.logger.warning("Detected gene_seqhash_id as PRIMARY KEY. Fixing schema...")
                    
                    # Create a new table with the correct schema
                    self.con.execute("""
                        CREATE TABLE gene_protein_map_new (
                            gene_seqhash_id VARCHAR,
                            protein_seqhash_id VARCHAR,
                            PRIMARY KEY (gene_seqhash_id, protein_seqhash_id),
                            FOREIGN KEY (protein_seqhash_id) REFERENCES sequences(seqhash_id)
                        )
                    """)
                    
                    # Copy data from old table to new table
                    self.con.execute("""
                        INSERT INTO gene_protein_map_new
                        SELECT gene_seqhash_id, protein_seqhash_id
                        FROM gene_protein_map
                    """)
                    
                    # Drop the old table
                    self.con.execute("DROP TABLE gene_protein_map")
                    
                    # Rename the new table
                    self.con.execute("ALTER TABLE gene_protein_map_new RENAME TO gene_protein_map")
                    
                    self.logger.info("Successfully fixed gene_protein_map schema")
        except Exception as e:
            self.logger.error(f"Error fixing gene_protein_map schema: {str(e)}")

    def _cleanup_incomplete_migrations(self):
        """Clean up any incomplete migrations by dropping backup tables."""
        try:
            tables = self.con.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name LIKE '%_backup'"
            ).fetchall()
            
            for table in tables:
                table_name = table[0]
                self.logger.warning(f"Found backup table {table_name}, attempting to drop")
                try:
                    self.con.execute(f"DROP TABLE IF EXISTS {table_name}")
                    self.logger.info(f"Successfully dropped backup table {table_name}")
                except Exception as e:
                    self.logger.error(f"Failed to drop backup table {table_name}: {str(e)}")
        except Exception as e:
            self.logger.error(f"Error checking for backup tables: {str(e)}")

    def _apply_duckdb_migrations_directly(self):
        """Apply migrations with DuckDB-specific syntax."""
        # Create core tables manually rather than relying on migration files
        
        # 1. Create sra_metadata table
        try:
            self.con.execute("""
                CREATE TABLE IF NOT EXISTS sra_metadata (
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
                    run_bases BIGINT,
                    run_published TIMESTAMP
                )
            """)
            self.logger.info("Created sra_metadata table")
        except Exception as e:
            self.logger.error(f"Failed to create sra_metadata table: {str(e)}")
        
        # 2. Create sequences table
        try:
            self.con.execute("""
                CREATE TABLE IF NOT EXISTS sequences (
                    seqhash_id VARCHAR PRIMARY KEY,
                    sequence VARCHAR NOT NULL,
                    sample_id VARCHAR NOT NULL,
                    assembly_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    is_representative BOOLEAN DEFAULT FALSE,
                    repseq_id VARCHAR NOT NULL,
                    length INTEGER NOT NULL,
                    FOREIGN KEY (sample_id) REFERENCES sra_metadata(sample_id)
                )
            """)
            self.logger.info("Created sequences table")
        except Exception as e:
            self.logger.error(f"Failed to create sequences table: {str(e)}")
        
        # 3. Create gene_protein_map table - without partial unique index
        try:
            self.con.execute("""
                CREATE TABLE IF NOT EXISTS gene_protein_map (
                    gene_seqhash_id VARCHAR NOT NULL,
                    protein_seqhash_id VARCHAR NOT NULL,
                    PRIMARY KEY (gene_seqhash_id, protein_seqhash_id),
                    FOREIGN KEY (protein_seqhash_id) REFERENCES sequences(seqhash_id)
                )
            """)
            
            # Create a regular index on gene_seqhash_id (no uniqueness constraints)
            self.con.execute("""
                CREATE INDEX IF NOT EXISTS idx_gene_protein_gene 
                ON gene_protein_map(gene_seqhash_id)
            """)
            
            self.logger.info("Created gene_protein_map table")
        except Exception as e:
            self.logger.error(f"Failed to create gene_protein_map table: {str(e)}")
        
        # 4. Create expression table - without foreign key to gene_protein_map
        try:
            self.con.execute("""
                CREATE TABLE IF NOT EXISTS expression (
                    gene_seqhash_id VARCHAR NOT NULL,
                    sample_id VARCHAR NOT NULL,
                    tpm DOUBLE,
                    num_reads DOUBLE,
                    effective_length DOUBLE,
                    PRIMARY KEY (gene_seqhash_id, sample_id),
                    FOREIGN KEY (sample_id) REFERENCES sra_metadata(sample_id)
                )
            """)
            self.logger.info("Created expression table")
        except Exception as e:
            self.logger.error(f"Failed to create expression table: {str(e)}")
            
        # 5. Create annotations table
        try:
            self.con.execute("""
                CREATE TABLE IF NOT EXISTS annotations (
                    seqhash_id VARCHAR PRIMARY KEY,
                    seed_ortholog VARCHAR,
                    evalue DOUBLE,
                    score DOUBLE,
                    eggnog_ogs VARCHAR,
                    max_annot_lvl VARCHAR,
                    cog_category VARCHAR,
                    description VARCHAR,
                    preferred_name VARCHAR,
                    sample_id VARCHAR,
                    FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id),
                    FOREIGN KEY (sample_id) REFERENCES sra_metadata(sample_id)
                )
            """)
            self.logger.info("Created annotations table")
        except Exception as e:
            self.logger.error(f"Failed to create annotations table: {str(e)}")
        
        # 6. Create go_terms table
        try:
            self.con.execute("""
                CREATE TABLE IF NOT EXISTS go_terms (
                    seqhash_id VARCHAR,
                    go_term VARCHAR,
                    PRIMARY KEY (seqhash_id, go_term),
                    FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id)
                )
            """)
            self.logger.info("Created go_terms table")
        except Exception as e:
            self.logger.error(f"Failed to create go_terms table: {str(e)}")
        
        # 7. Create ec_numbers table
        try:
            self.con.execute("""
                CREATE TABLE IF NOT EXISTS ec_numbers (
                    seqhash_id VARCHAR,
                    ec_number VARCHAR,
                    PRIMARY KEY (seqhash_id, ec_number),
                    FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id)
                )
            """)
            self.logger.info("Created ec_numbers table")
        except Exception as e:
            self.logger.error(f"Failed to create ec_numbers table: {str(e)}")
        
        # 8. Create clusters table
        try:
            self.con.execute("""
                CREATE TABLE IF NOT EXISTS clusters (
                    cluster_id VARCHAR PRIMARY KEY,
                    representative_seqhash_id VARCHAR NOT NULL,
                    size INTEGER NOT NULL,
                    FOREIGN KEY (representative_seqhash_id) REFERENCES sequences(seqhash_id)
                )
            """)
            self.logger.info("Created clusters table")
        except Exception as e:
            self.logger.error(f"Failed to create clusters table: {str(e)}")
        
        # 9. Create cluster_members table
        try:
            self.con.execute("""
                CREATE TABLE IF NOT EXISTS cluster_members (
                    seqhash_id VARCHAR,
                    cluster_id VARCHAR,
                    PRIMARY KEY (seqhash_id, cluster_id),
                    FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id),
                    FOREIGN KEY (cluster_id) REFERENCES clusters(cluster_id)
                )
            """)
            self.logger.info("Created cluster_members table")
        except Exception as e:
            self.logger.error(f"Failed to create cluster_members table: {str(e)}")
        
        # 10. Create kegg_info table
        try:
            self.con.execute("""
                CREATE TABLE IF NOT EXISTS kegg_info (
                    seqhash_id VARCHAR PRIMARY KEY,
                    kegg_ko VARCHAR,
                    kegg_pathway VARCHAR,
                    kegg_module VARCHAR,
                    kegg_reaction VARCHAR,
                    kegg_rclass VARCHAR,
                    brite VARCHAR,
                    kegg_tc VARCHAR,
                    cazy VARCHAR,
                    bigg_reaction VARCHAR,
                    pfams VARCHAR,
                    sample_id VARCHAR,
                    FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id),
                    FOREIGN KEY (sample_id) REFERENCES sra_metadata(sample_id)
                )
            """)
            self.logger.info("Created kegg_info table")
        except Exception as e:
            self.logger.error(f"Failed to create kegg_info table: {str(e)}")


    def _check_incomplete_migration(self):
        """Check if there's an incomplete migration by looking for backup tables."""
        try:
            tables = self.con.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name LIKE '%_backup'"
            ).fetchall()
            
            if tables:
                return tables[0][0]  # Return the first backup table name
            return None
        except Exception:
            return None

    def _cleanup_incomplete_migrations(self):
        """Clean up any incomplete migrations by dropping backup tables."""
        try:
            tables = self.con.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name LIKE '%_backup'"
            ).fetchall()
            
            for table in tables:
                table_name = table[0]
                self.logger.warning(f"Found backup table {table_name}, attempting to drop")
                try:
                    self.con.execute(f"DROP TABLE IF EXISTS {table_name}")
                    self.logger.info(f"Successfully dropped backup table {table_name}")
                except Exception as e:
                    self.logger.error(f"Failed to drop backup table {table_name}: {str(e)}")
        except Exception as e:
            self.logger.error(f"Error checking for backup tables: {str(e)}")

            
    def _get_sample_paths(self, sample_id: str) -> SamplePaths:
        """Get paths to sample data files."""
        sample_dir = self.output_dir / sample_id
        return SamplePaths(
            sequences=sample_dir / f"transdecoder/{sample_id}.pep",
            annotations=sample_dir / f"eggnog/{sample_id}.emapper.annotations",
            expression=sample_dir / f"quants/{sample_id}.quant.json",
        )

    def _process_migration_file(self, migration_file):
        """Process a single migration file with error handling."""
        self.logger.info(f"Processing migration: {migration_file.name}")
        
        try:
            with open(migration_file) as f:
                migration_sql = f.read()
            
            # Split the migration into individual statements
            statements = self._split_sql_statements(migration_sql)
            
            # Execute each statement separately
            for i, statement in enumerate(statements):
                if not statement.strip():
                    continue
                    
                try:
                    # Add safety to CREATE TABLE statements
                    if statement.strip().upper().startswith("CREATE TABLE") and "IF NOT EXISTS" not in statement.upper():
                        table_name = self._extract_table_name(statement)
                        
                        # Check if table already exists
                        table_exists = self.con.execute(
                            f"SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='{table_name}'"
                        ).fetchone()[0] > 0
                        
                        if table_exists:
                            self.logger.warning(f"Table {table_name} already exists, skipping creation")
                            continue
                        else:
                            # Add IF NOT EXISTS
                            statement = statement.replace("CREATE TABLE", "CREATE TABLE IF NOT EXISTS")
                            statement = statement.replace("create table", "CREATE TABLE IF NOT EXISTS")
                    
                    # Execute the statement
                    self.con.execute(statement)
                    self.logger.info(f"Executed statement {i+1}/{len(statements)}")
                    
                except Exception as e:
                    self.logger.error(f"Error executing statement {i+1}: {str(e)}")
                    self.logger.error(f"Statement: {statement}")
                    
                    # Check if this is a "table already exists" error
                    if "already exists" in str(e):
                        self.logger.warning("Continuing despite 'already exists' error")
                        continue
                    else:
                        raise
            
            # Record that we applied this migration
            self.con.execute(
                """
                INSERT INTO schema_version (version, migration_name)
                VALUES (
                    (SELECT COALESCE(MAX(version), 0) + 1 FROM schema_version),
                    ?
                )
            """,
                [migration_file.name],
            )
            self.logger.info(f"Successfully recorded migration: {migration_file.name}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to process migration {migration_file.name}: {str(e)}")
            return False

    def _split_sql_statements(self, sql):
        """Split a SQL script into individual statements."""
        # Simple splitting by semicolon - may need to be more sophisticated
        # for complex SQL with semicolons in strings, etc.
        statements = []
        current = []
        
        for line in sql.split('\n'):
            # Skip comments
            if line.strip().startswith('--'):
                continue
                
            current.append(line)
            
            if ';' in line:
                statements.append('\n'.join(current))
                current = []
        
        # Add any remaining content without a semicolon
        if current:
            statements.append('\n'.join(current))
        
        return statements

    def _extract_table_name(self, create_statement):
        """Extract the table name from a CREATE TABLE statement."""
        # This is a simple approach - might need to be more robust
        statement = create_statement.strip().upper()
        
        if "CREATE TABLE IF NOT EXISTS" in statement:
            after_create = statement.split("CREATE TABLE IF NOT EXISTS")[1]
        else:
            after_create = statement.split("CREATE TABLE")[1]
        
        # Extract the table name (the first word after CREATE TABLE)
        table_name = after_create.strip().split()[0].strip('"').strip('`').strip("'")
        
        # Remove any schema prefix
        if '.' in table_name:
            table_name = table_name.split('.')[-1]
            
        return table_name

    def _fetch_sra_metadata(self, sample_id: str) -> Dict:
        """Fetch metadata for a single SRA ID."""
        self.logger.info(f"Fetching metadata for {sample_id}")
        try:
            info = get_sra_info(sample_id)
            if not isinstance(info, dict):
                raise ValueError(f"Invalid response for {sample_id}")

            return {
                "sample_id": sample_id,
                "organism": info.get("organism"),
                "study_title": info.get("study_title"),
                "study_abstract": info.get("study_abstract"),
                "bioproject": info.get("bioproject"),
                "biosample": info.get("biosample"),
                "library_strategy": info.get("library", {}).get("strategy"),
                "library_source": info.get("library", {}).get("source"),
                "library_selection": info.get("library", {}).get("selection"),
                "library_layout": info.get("library", {}).get("layout"),
                "instrument": info.get("instrument"),
                "run_spots": info.get("run", {}).get("spots"),
                "run_bases": info.get("run", {}).get("bases"),
                "run_published": info.get("run", {}).get("published"),
            }
        except Exception as e:
            self.logger.error(f"Failed to fetch metadata for {sample_id}: {e}")
            raise

    def _batch_insert_sequences(self, sequences: List[dict]):
        """Helper method to insert sequences in batches."""
        if not sequences:
            return
        
        try:
            # Create a pandas DataFrame
            df = pd.DataFrame(sequences)
            
            # Register the DataFrame with DuckDB
            self.con.register("temp_sequences_df", df)
            
            # Insert using SQL
            self.con.execute("""
                INSERT OR IGNORE INTO sequences 
                SELECT * FROM temp_sequences_df
            """)
            
            self.logger.info(f"Inserted {len(sequences)} sequences")
        except Exception as e:
            self.logger.error(f"Error in batch insert sequences: {str(e)}")
            # Re-raise to handle in the calling method
            raise

    def _extract_gene_id(self, protein_id):
        """
        Extract gene ID from protein ID by removing .p1, .p2, etc. suffixes.
        Uses string operations instead of regex to avoid escape sequence issues.
        """
        # Split by '.p' and take the first part
        if '.p' in protein_id:
            return protein_id.split('.p')[0]
        return protein_id

    def _load_sequences(self, sequence_path: Path, sample_id: str, duplicates: list) -> int:
        """Load sequences from FASTA file."""
        self.logger.info(f"Loading sequences from {sequence_path}")

        existing_seqs = set(
            self.con.execute("SELECT seqhash_id FROM sequences").df()["seqhash_id"]
        )
        sequences = []
        gene_protein_maps = []
        sequences_loaded = 0
        skipped_count = 0

        # Use with statement to properly close the file
        with open(sequence_path, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                # Example ID: v1_DLS_xxxxx.p1
                protein_seqhash_id = record.id.split()[0]

                # Extract gene seqhash from protein seqhash by removing .p1, .p2, etc.
                # Using string operations instead of regex
                gene_seqhash_id = self._extract_gene_id(protein_seqhash_id)

                if protein_seqhash_id in existing_seqs:
                    duplicates.append(
                        {
                            "seqhash_id": protein_seqhash_id,
                            "sample_id": sample_id,
                            "length": len(record.seq),
                        }
                    )
                    skipped_count += 1
                    continue

                sequences.append(
                    {
                        "seqhash_id": protein_seqhash_id,
                        "sequence": str(record.seq),
                        "sample_id": sample_id,
                        "assembly_date": datetime.now(),
                        "is_representative": False,
                        "repseq_id": protein_seqhash_id,  # initialize the repseq_id to be the same as seqhash_id
                        "length": len(record.seq),
                    }
                )

                # Add the gene-protein mapping
                gene_protein_maps.append(
                    {
                        "gene_seqhash_id": gene_seqhash_id,
                        "protein_seqhash_id": protein_seqhash_id,
                    }
                )

                sequences_loaded += 1

                if len(sequences) % 1000 == 0:
                    self._batch_insert_sequences(sequences)
                    self._batch_insert_gene_protein_maps(gene_protein_maps)
                    sequences = []
                    gene_protein_maps = []

        if sequences:
            self._batch_insert_sequences(sequences)

        if gene_protein_maps:
            self._batch_insert_gene_protein_maps(gene_protein_maps)

        self.logger.info(f"Loaded {sequences_loaded} new sequences for {sample_id}")
        if skipped_count > 0:
            self.logger.info(f"Skipped {skipped_count} duplicate sequences")

        return sequences_loaded

    def _batch_insert_gene_protein_maps(self, gene_protein_maps: List[dict]):
        """Helper method to insert gene-protein mappings in batches."""
        if not gene_protein_maps:
            return

        try:
            # Create a pandas DataFrame
            df = pd.DataFrame(gene_protein_maps)
            
            # Register the DataFrame with DuckDB
            self.con.register("temp_gene_protein_map_df", df)
            
            # First, we'll create a table of existing mappings
            self.con.execute("""
                CREATE TEMP TABLE existing_mappings AS
                SELECT gene_seqhash_id, protein_seqhash_id
                FROM gene_protein_map
            """)
            
            # Then, create a table with only the new mappings
            # In DuckDB, we can't do a tuple comparison like "(a,b) NOT IN (SELECT c,d...)"
            # So we'll use a LEFT JOIN with a NULL filter instead
            self.con.execute("""
                CREATE TEMP TABLE new_mappings AS
                SELECT t.gene_seqhash_id, t.protein_seqhash_id
                FROM temp_gene_protein_map_df t
                LEFT JOIN existing_mappings e ON 
                    t.gene_seqhash_id = e.gene_seqhash_id AND
                    t.protein_seqhash_id = e.protein_seqhash_id
                WHERE e.gene_seqhash_id IS NULL
            """)
            
            # Get count of new mappings
            new_count = self.con.execute("SELECT COUNT(*) FROM new_mappings").fetchone()[0]
            
            if new_count > 0:
                # Insert the new mappings
                self.con.execute("""
                    INSERT INTO gene_protein_map
                    SELECT * FROM new_mappings
                """)
                self.logger.info(f"Inserted {new_count} new gene-protein mappings")
            
            # Clean up
            self.con.execute("DROP TABLE IF EXISTS existing_mappings")
            self.con.execute("DROP TABLE IF EXISTS new_mappings")
            
        except Exception as e:
            self.logger.error(f"Error in batch insert gene-protein maps: {str(e)}")
            # Clean up temporary tables if they exist
            try:
                self.con.execute("DROP TABLE IF EXISTS existing_mappings")
                self.con.execute("DROP TABLE IF EXISTS new_mappings")
            except:
                pass
            # Re-raise to handle in the calling method
            raise
        
    def _load_annotations(self, annotation_path: Path, sample_id: str) -> int:
        """Load and process annotations from eggNOG file."""
        self.logger.info(f"Loading annotations from {annotation_path}")

        column_names = [
            "query",
            "seed_ortholog",
            "evalue",
            "score",
            "eggNOG_OGs",
            "max_annot_lvl",
            "COG_category",
            "Description",
            "Preferred_name",
            "GOs",
            "EC",
            "KEGG_ko",
            "KEGG_Pathway",
            "KEGG_Module",
            "KEGG_Reaction",
            "KEGG_rclass",
            "BRITE",
            "KEGG_TC",
            "CAZy",
            "BiGG_Reaction",
            "PFAMs",
        ]

        try:
            # Create temporary table
            columns_def = ", ".join([f'"{name}" VARCHAR' for name in column_names])
            self.con.execute(f"CREATE TABLE temp_annotations ({columns_def})")

            # Load annotation data
            self.con.execute(
                f"""
                INSERT INTO temp_annotations
                SELECT * FROM read_csv_auto(
                    '{annotation_path}',
                    sep='\t',
                    header=False,
                    names={column_names},
                    comment='#'
                )
            """
            )

            # Process annotations and get count
            return self._process_annotations(sample_id)

        finally:
            self.con.execute("DROP TABLE IF EXISTS temp_annotations")

    def _process_annotations(self, sample_id: str) -> int:
        """Process annotations from temporary table into final tables."""
        # Make sure the kegg_info table exists
        self._ensure_kegg_info_table_exists()
        
        # Main annotations
        self.con.execute(
            f"""
            INSERT INTO annotations
            SELECT 
                query as seqhash_id,
                seed_ortholog,
                TRY_CAST(evalue AS DOUBLE) as evalue,
                TRY_CAST(score AS DOUBLE) as score,
                "eggNOG_OGs" as eggnog_ogs,
                max_annot_lvl,
                "COG_category" as cog_category,
                "Description" as description,
                "Preferred_name" as preferred_name,
                '{sample_id}' as sample_id
            FROM temp_annotations
            WHERE query IN (SELECT seqhash_id FROM sequences)
            AND query NOT IN (SELECT seqhash_id FROM annotations)
        """
        )

        annotations_loaded = self.con.execute(
            "SELECT COUNT(*) FROM annotations WHERE sample_id = ?", [sample_id]
        ).fetchone()[0]

        # GO terms
        self.con.execute(
            """
            INSERT INTO go_terms
            SELECT DISTINCT
                query as seqhash_id,
                UNNEST(STRING_SPLIT(NULLIF("GOs", '-'), ',')) as go_term
            FROM temp_annotations
            WHERE query IN (SELECT seqhash_id FROM sequences)
            AND query NOT IN (SELECT seqhash_id FROM go_terms)
            AND "GOs" IS NOT NULL AND "GOs" != '-'
        """
        )

        # EC numbers
        self.con.execute(
            """
            INSERT INTO ec_numbers
            SELECT DISTINCT
                query as seqhash_id,
                UNNEST(STRING_SPLIT(NULLIF("EC", '-'), ',')) as ec_number
            FROM temp_annotations
            WHERE query IN (SELECT seqhash_id FROM sequences)
            AND query NOT IN (SELECT seqhash_id FROM ec_numbers)
            AND "EC" IS NOT NULL AND "EC" != '-'
        """
        )

        # Add KEGG information in a separate transaction
        try:
            with self.transaction():
                self.con.execute(
                    f"""
                    INSERT INTO kegg_info
                    SELECT DISTINCT
                        query as seqhash_id,
                        "KEGG_ko" as kegg_ko,
                        "KEGG_Pathway" as kegg_pathway,
                        "KEGG_Module" as kegg_module,
                        "KEGG_Reaction" as kegg_reaction,
                        "KEGG_rclass" as kegg_rclass,
                        "BRITE" as brite,
                        "KEGG_TC" as kegg_tc,
                        "CAZy" as cazy,
                        "BiGG_Reaction" as bigg_reaction,
                        "PFAMs" as pfams,
                        '{sample_id}' as sample_id
                    FROM temp_annotations
                    WHERE query IN (SELECT seqhash_id FROM sequences)
                    AND query NOT IN (SELECT seqhash_id FROM kegg_info)
                    AND ("KEGG_ko" != '-' OR "KEGG_Pathway" != '-' OR 
                        "KEGG_Module" != '-' OR "KEGG_Reaction" != '-' OR
                        "KEGG_rclass" != '-' OR "BRITE" != '-' OR 
                        "KEGG_TC" != '-' OR "CAZy" != '-' OR 
                        "BiGG_Reaction" != '-' OR "PFAMs" != '-')
                """
                )
                self.logger.info(f"Added KEGG information for {sample_id}")
        except Exception as e:
            self.logger.error(f"Failed to add KEGG information: {str(e)}")
            # Continue despite KEGG insertion failure
        
        return annotations_loaded

    def _ensure_kegg_info_table_exists(self):
        """Create the kegg_info table if it doesn't exist."""
        try:
            with self.transaction():
                self.con.execute("""
                    CREATE TABLE IF NOT EXISTS kegg_info (
                        seqhash_id VARCHAR PRIMARY KEY,
                        kegg_ko VARCHAR,
                        kegg_pathway VARCHAR,
                        kegg_module VARCHAR,
                        kegg_reaction VARCHAR,
                        kegg_rclass VARCHAR,
                        brite VARCHAR,
                        kegg_tc VARCHAR,
                        cazy VARCHAR,
                        bigg_reaction VARCHAR,
                        pfams VARCHAR,
                        sample_id VARCHAR,
                        FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id),
                        FOREIGN KEY (sample_id) REFERENCES sra_metadata(sample_id)
                    )
                """)
                self.logger.info("Created kegg_info table")
        except Exception as e:
            self.logger.error(f"Failed to create kegg_info table: {str(e)}")



    def _check_sample_exists(self, sample_id: str) -> bool:
        """Check if a sample ID already exists in the database."""
        result = self.con.execute(
            "SELECT COUNT(*) FROM sra_metadata WHERE sample_id = ?", [sample_id]
        ).fetchone()[0]
        return result > 0

    # def process_sample(self, sample_id: str) -> ProcessingResult:
    #     """Process a single sample: fetch metadata and load sequence data."""
    #     try:
    #         paths = self._get_sample_paths(sample_id)

    #         if not paths.sequences.exists():
    #             raise FileNotFoundError(f"Sequence file not found: {paths.sequences}")

    #         metadata = self._fetch_sra_metadata(sample_id)

    #         def _process():
    #             # Insert metadata
    #             if metadata:
    #                 self.con.execute("DELETE FROM sra_metadata WHERE sample_id = ?", [sample_id])
    #                 self.con.execute("""
    #                     INSERT INTO sra_metadata (
    #                         sample_id, organism, study_title, study_abstract,
    #                         bioproject, biosample, library_strategy, library_source,
    #                         library_selection, library_layout, instrument,
    #                         run_spots, run_bases, run_published
    #                     ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    #                 """, [
    #                     sample_id,
    #                     metadata.get('organism'),
    #                     metadata.get('study_title'),
    #                     metadata.get('study_abstract'),
    #                     metadata.get('bioproject'),
    #                     metadata.get('biosample'),
    #                     metadata.get('library', {}).get('strategy'),
    #                     metadata.get('library', {}).get('source'),
    #                     metadata.get('library', {}).get('selection'),
    #                     metadata.get('library', {}).get('layout'),
    #                     metadata.get('instrument'),
    #                     metadata.get('run', {}).get('spots'),
    #                     metadata.get('run', {}).get('bases'),
    #                     metadata.get('run', {}).get('published')
    #                 ])

    #             # Load sequences
    #             duplicates = []
    #             sequences_loaded = self._load_sequences(paths.sequences, sample_id, duplicates)

    #             # Load annotations if available
    #             annotations_loaded = 0
    #             if paths.annotations.exists():
    #                 annotations_loaded = self._load_annotations(paths.annotations, sample_id)
    #             else:
    #                 self.logger.warning(f"No annotation file found: {paths.annotations}")

    #             return sequences_loaded, annotations_loaded, duplicates

    #         # Execute everything in a single transaction
    #         sequences_loaded, annotations_loaded, duplicates = self.execute_transaction(_process)

    #         return ProcessingResult(
    #             sample_id=sample_id,
    #             status='success',
    #             sequences_loaded=sequences_loaded,
    #             annotations_loaded=annotations_loaded,
    #             duplicates=len(duplicates)
    #         )

    #     except Exception as e:
    #         self.logger.error(f"Error processing {sample_id}: {str(e)}")
    #         return ProcessingResult(
    #             sample_id=sample_id,
    #             status='error',
    #             error=str(e)
    #         )
    def _load_expression(self, expression_path: Path, sample_id: str) -> int:
        """Load expression data from Salmon quantification JSON file."""
        self.logger.info(f"Loading expression data from {expression_path}")

        if not expression_path.exists():
            self.logger.warning(f"Expression file not found: {expression_path}")
            return 0

        try:
            # First, try to validate the JSON file
            try:
                import json
                import pandas as pd

                with open(expression_path, "r") as f:
                    data = json.load(f)
                    self.logger.info(f"Verified JSON file format: {len(data)} records")

                # Debug log the first few entries to understand structure
                if len(data) > 0:
                    self.logger.info(f"First entry structure: {list(data[0].keys())}")
            except Exception as e:
                self.logger.error(f"Error validating JSON file: {e}")
                return 0

            # Load JSON file and insert into expression table
            import json
            import pandas as pd

            # Read JSON from file
            with open(expression_path, "r") as f:
                json_data = json.load(f)

            # Convert to pandas DataFrame
            df = pd.DataFrame(json_data)

            # Rename columns to match our schema
            rename_cols = {}
            if "Name" in df.columns:
                rename_cols["Name"] = "gene_seqhash_id"
            if "TPM" in df.columns:
                rename_cols["TPM"] = "tpm"
            if "NumReads" in df.columns:
                rename_cols["NumReads"] = "num_reads"
            if "EffectiveLength" in df.columns:
                rename_cols["EffectiveLength"] = "effective_length"
            
            df_renamed = df.rename(columns=rename_cols)

            # Add sample_id column
            df_renamed["sample_id"] = sample_id

            # Step 1: Create a mapping table for gene IDs to protein IDs
            with self.transaction():
                self.logger.info("Creating gene-to-protein mapping table")
                self.con.execute("""
                    CREATE TEMP TABLE genes_from_proteins AS
                    SELECT DISTINCT
                        CASE 
                            WHEN CONTAINS(seqhash_id, '.p') 
                            THEN SPLIT_PART(seqhash_id, '.p', 1) 
                            ELSE seqhash_id 
                        END as gene_seqhash_id,
                        seqhash_id as protein_seqhash_id
                    FROM sequences
                    WHERE sample_id = ?
                """, [sample_id])

            # Step 2: Filter expression data to only include genes with sequence entries
            self.logger.info("Filtering expression data to match genes in sequences")
            
            # Get list of gene IDs from sequences
            gene_ids_from_sequences = self.con.execute("""
                SELECT DISTINCT gene_seqhash_id
                FROM genes_from_proteins
            """).fetchall()
            gene_ids_set = {row[0] for row in gene_ids_from_sequences}
            
            # Filter expression data to only include genes in the sequences table
            df_filtered = df_renamed[df_renamed['gene_seqhash_id'].isin(gene_ids_set)]
            self.logger.info(f"Filtered expression data: {len(df_filtered)} of {len(df_renamed)} entries match genes")
            
            if len(df_filtered) == 0:
                self.logger.warning("No matching genes found in expression data")
                self.con.execute("DROP TABLE IF EXISTS genes_from_proteins")
                return 0

            # Step 3: Get all gene-protein mappings (keeps ALL proteins for each gene)
            gene_protein_mappings = self.con.execute("""
                SELECT gene_seqhash_id, protein_seqhash_id
                FROM genes_from_proteins
                ORDER BY gene_seqhash_id, protein_seqhash_id
            """).fetchall()
            
            # Create DataFrame of all gene-protein mappings
            mappings_list = [
                {"gene_seqhash_id": gene_id, "protein_seqhash_id": protein_id} 
                for gene_id, protein_id in gene_protein_mappings
            ]
            
            # Step 4: Insert gene-protein mappings
            if mappings_list:
                self.logger.info(f"Found {len(mappings_list)} total gene-protein mappings")
                
                # Check which mappings already exist
                existing_mappings = self.con.execute("""
                    SELECT gene_seqhash_id, protein_seqhash_id
                    FROM gene_protein_map
                """).fetchall()
                
                # Create a set of existing mapping tuples for faster lookup
                existing_mapping_set = {(gene, protein) for gene, protein in existing_mappings}
                
                # Filter to only new mappings
                new_mappings = [
                    {"gene_seqhash_id": m["gene_seqhash_id"], "protein_seqhash_id": m["protein_seqhash_id"]}
                    for m in mappings_list
                    if (m["gene_seqhash_id"], m["protein_seqhash_id"]) not in existing_mapping_set
                ]
                
                if new_mappings:
                    # Insert new mappings in a single transaction
                    with self.transaction():
                        df_mappings = pd.DataFrame(new_mappings)
                        self.con.register("new_gene_protein_mappings", df_mappings)
                        self.con.execute("""
                            INSERT INTO gene_protein_map
                            SELECT * FROM new_gene_protein_mappings
                        """)
                        self.logger.info(f"Inserted {len(new_mappings)} new gene-protein mappings")
            
            # Step 5: Insert expression data
            self.logger.info(f"Preparing to insert expression data")
            # Check which expression records already exist
            existing_expressions = self.con.execute("""
                SELECT gene_seqhash_id
                FROM expression
                WHERE sample_id = ?
            """, [sample_id]).fetchall()
            existing_expr_genes = {row[0] for row in existing_expressions}
            
            # Filter to only include new expression records
            df_new_expr = df_filtered[~df_filtered['gene_seqhash_id'].isin(existing_expr_genes)]
            
            if len(df_new_expr) > 0:
                self.logger.info(f"Inserting {len(df_new_expr)} new expression records")
                with self.transaction():
                    self.con.register("new_expression_data", df_new_expr)
                    self.con.execute("""
                        INSERT INTO expression (gene_seqhash_id, sample_id, tpm, num_reads, effective_length)
                        SELECT gene_seqhash_id, sample_id, tpm, num_reads, effective_length
                        FROM new_expression_data
                    """)
            else:
                self.logger.info("No new expression records to insert")
            
            # Clean up
            self.con.execute("DROP TABLE IF EXISTS genes_from_proteins")

            # Count how many expression records were loaded
            expression_loaded = self.con.execute(
                "SELECT COUNT(*) FROM expression WHERE sample_id = ?", [sample_id]
            ).fetchone()[0]

            self.logger.info(
                f"Loaded {expression_loaded} expression records for {sample_id}"
            )
            return expression_loaded

        except Exception as e:
            self.logger.error(f"Error loading expression data: {e}")
            # Clean up any temporary tables
            try:
                self.con.execute("DROP TABLE IF EXISTS genes_from_proteins")
            except:
                pass
            # Return 0 instead of raising to allow processing to continue
            return 0
    def process_sample(self, sample_id: str) -> ProcessingResult:
        """Process a single sample: fetch metadata and load sequence data."""
        try:
            self.logger.info(f"Starting to process sample {sample_id}")
            
            # Check if sample already exists - handle if table doesn't exist yet
            try:
                exists = self._check_sample_exists(sample_id)
                if exists:
                    self.logger.info(f"Sample {sample_id} already exists in database")
                    return ProcessingResult(
                        sample_id=sample_id,
                        status="error",
                        error=f"Sample {sample_id} already exists in the database",
                    )
            except Exception as e:
                self.logger.warning(f"Error checking if sample exists: {str(e)}")
                # Continue anyway - if table doesn't exist, sample can't already exist

            paths = self._get_sample_paths(sample_id)
            self.logger.info(f"Sample paths: {paths}")

            if not paths.sequences.exists():
                self.logger.error(f"Sequence file not found: {paths.sequences}")
                raise FileNotFoundError(f"Sequence file not found: {paths.sequences}")

            # Fetch metadata
            try:
                self.logger.info(f"Fetching metadata for {sample_id}")
                metadata = self._fetch_sra_metadata(sample_id)
            except Exception as e:
                self.logger.error(f"Failed to fetch metadata: {str(e)}")
                # Create minimal metadata to allow processing to continue
                metadata = {"sample_id": sample_id}

            # Process in multiple steps with separate transactions for robustness
            
            # Step 1: Insert metadata
            try:
                with self.transaction():
                    self.logger.info(f"Inserting metadata for {sample_id}")
                    self.con.execute("DELETE FROM sra_metadata WHERE sample_id = ?", [sample_id])
                    self.con.execute(
                        """
                        INSERT INTO sra_metadata (
                            sample_id, organism, study_title, study_abstract,
                            bioproject, biosample, library_strategy, library_source,
                            library_selection, library_layout, instrument,
                            run_spots, run_bases, run_published
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                        [
                            sample_id,
                            metadata.get("organism"),
                            metadata.get("study_title"),
                            metadata.get("study_abstract"),
                            metadata.get("bioproject"),
                            metadata.get("biosample"),
                            metadata.get("library_strategy"),
                            metadata.get("library_source"),
                            metadata.get("library_selection"),
                            metadata.get("library_layout"),
                            metadata.get("instrument"),
                            metadata.get("run_spots"),
                            metadata.get("run_bases"),
                            metadata.get("run_published"),
                        ],
                    )
                    self.logger.info(f"Successfully inserted metadata for {sample_id}")
            except Exception as e:
                self.logger.error(f"Failed to insert metadata: {str(e)}")
                return ProcessingResult(
                    sample_id=sample_id,
                    status="error",
                    error=f"Failed to insert metadata: {str(e)}",
                )

            # Step 2: Load sequences
            duplicates = []
            sequences_loaded = 0
            try:
                with self.transaction():
                    self.logger.info(f"Loading sequences for {sample_id}")
                    sequences_loaded = self._load_sequences(
                        paths.sequences, sample_id, duplicates
                    )
                    self.logger.info(f"Loaded {sequences_loaded} sequences for {sample_id}")
            except Exception as e:
                self.logger.error(f"Failed to load sequences: {str(e)}")
                return ProcessingResult(
                    sample_id=sample_id,
                    status="error",
                    error=f"Failed to load sequences: {str(e)}",
                )

            # Check if we actually loaded sequences successfully
            seq_count = 0
            try:
                seq_count = self.con.execute(
                    "SELECT COUNT(*) FROM sequences WHERE sample_id = ?", [sample_id]
                ).fetchone()[0]
                self.logger.info(f"Verified sequence loading: {seq_count} sequences in database")
            except Exception as e:
                self.logger.error(f"Failed to verify sequence count: {str(e)}")

            # Step 3: Load annotations if available
            annotations_loaded = 0
            if paths.annotations.exists():
                try:
                    with self.transaction():
                        self.logger.info(f"Loading annotations for {sample_id}")
                        annotations_loaded = self._load_annotations(
                            paths.annotations, sample_id
                        )
                        self.logger.info(f"Loaded {annotations_loaded} annotations for {sample_id}")
                except Exception as e:
                    self.logger.error(f"Failed to load annotations: {str(e)}")
                    # Continue anyway - annotations are optional
            else:
                self.logger.warning(f"No annotation file found: {paths.annotations}")

            # Step 4: Load expression data if available
            expression_loaded = 0
            if paths.expression.exists() and seq_count > 0:
                try:
                    with self.transaction():
                        self.logger.info(f"Loading expression data for {sample_id}")
                        expression_loaded = self._load_expression(
                            paths.expression, sample_id
                        )
                        self.logger.info(f"Loaded {expression_loaded} expression records for {sample_id}")
                except Exception as e:
                    self.logger.error(f"Failed to load expression data: {str(e)}")
                    # Continue anyway - expression data is optional
            else:
                if not paths.expression.exists():
                    self.logger.warning(f"No expression file found: {paths.expression}")
                elif seq_count == 0:
                    self.logger.warning("Cannot load expression data: no sequences loaded")

            self.logger.info(f"Sample processing complete: {sample_id}")
            return ProcessingResult(
                sample_id=sample_id,
                status="success",
                sequences_loaded=sequences_loaded,
                annotations_loaded=annotations_loaded,
                duplicates=len(duplicates),
            )

        except Exception as e:
            self.logger.error(f"Error processing {sample_id}: {str(e)}")
            return ProcessingResult(
                sample_id=sample_id, 
                status="error", 
                error=str(e)
            )
    def build_database(self, sample_ids: List[str]) -> pd.DataFrame:
        """Build complete database from list of SRA IDs."""
        # Check if database already has data
        if self._check_database_exists():
            raise ValueError(
                "Database already contains data. Use update_database() to add new samples."
            )

        # Initialize database
        self.init_database()

        # Process each sample
        results = []
        for sample_id in sample_ids:
            result = self.process_sample(sample_id)
            results.append(result)

            # Log progress
            if result.status == "success":
                self.logger.info(
                    f"Processed {sample_id}: "
                    f"{result.sequences_loaded} sequences, "
                    f"{result.annotations_loaded} annotations, "
                    f"{result.duplicates} duplicates"
                )
            else:
                self.logger.error(f"Failed to process {sample_id}: {result.error}")

        return pd.DataFrame([vars(r) for r in results])

    def _check_database_exists(self) -> bool:
        """Check if database has existing data."""
        tables = self.con.execute(
            """
            SELECT name FROM sqlite_master 
            WHERE type='table' AND name='sra_metadata'
        """
        ).fetchall()
        if not tables:
            return False
        count = self.con.execute("SELECT COUNT(*) FROM sra_metadata").fetchone()[0]
        return count > 0

    def update_database(self, sample_ids: List[str]) -> pd.DataFrame:
        """Update existing database with new samples."""
        if not self._check_database_exists():
            raise ValueError(
                "Database not initialized. Use build_database() for new databases."
            )

        results = []
        for sample_id in sample_ids:
            # Check for duplicates before processing
            if self._check_sample_exists(sample_id):
                results.append(
                    ProcessingResult(
                        sample_id=sample_id,
                        status="error",
                        error=f"Sample {sample_id} already exists in the database",
                    )
                )
                self.logger.warning(f"Skipping {sample_id}: already exists in database")
                continue

            result = self.process_sample(sample_id)
            results.append(result)

            if result.status == "success":
                self.logger.info(
                    f"Added {sample_id}: "
                    f"{result.sequences_loaded} sequences, "
                    f"{result.annotations_loaded} annotations"
                )
            else:
                self.logger.error(f"Failed to add {sample_id}: {result.error}")

        return pd.DataFrame([vars(r) for r in results])

    def load_clusters_from_tsv(self, tsv_path: str):
        """Load cluster information from MMseqs2 cluster update TSV."""
        self.logger.info(f"Loading cluster data from {tsv_path}")

        with self.transaction():
            try:
                # Create temporary table for TSV data
                self.con.execute(
                    """
                    CREATE TEMP TABLE temp_clusters AS 
                    SELECT 
                        representative as representative_seqhash_id,
                        member as seqhash_id
                    FROM read_csv_auto(?, sep='\t', header=False, 
                                    names=['representative', 'member'])
                    WHERE representative IN (SELECT seqhash_id FROM sequences)
                    AND member IN (SELECT seqhash_id FROM sequences)
                    """,
                    [tsv_path],
                )

                # Insert clusters
                self.con.execute(
                    """
                    WITH cluster_info AS (
                        SELECT 
                            representative_seqhash_id,
                            ROW_NUMBER() OVER (ORDER BY representative_seqhash_id) as cluster_num,
                            COUNT(*) as size
                        FROM temp_clusters
                        GROUP BY representative_seqhash_id
                    )
                    INSERT INTO clusters (cluster_id, representative_seqhash_id, size)
                    SELECT 
                        'CLUSTER_' || cluster_num as cluster_id,
                        representative_seqhash_id,
                        size
                    FROM cluster_info
                    """
                )

                # Insert cluster members
                self.con.execute(
                    """
                    INSERT INTO cluster_members (seqhash_id, cluster_id)
                    SELECT 
                        tc.seqhash_id,
                        c.cluster_id
                    FROM temp_clusters tc
                    JOIN clusters c ON tc.representative_seqhash_id = c.representative_seqhash_id
                    """
                )

                # Update representative status
                self.con.execute(
                    """
                    UPDATE sequences
                    SET is_representative = TRUE
                    WHERE seqhash_id IN (SELECT representative_seqhash_id FROM clusters)
                    """
                )

            finally:
                self.con.execute("DROP TABLE IF EXISTS temp_clusters")

    def get_database_summary(self) -> pd.DataFrame:
        """Get summary statistics about the database."""
        try:
            # Check if all tables exist first
            tables = self.con.execute("""
                SELECT name FROM sqlite_master 
                WHERE type='table' 
                AND name IN ('sequences', 'sra_metadata', 'annotations', 'go_terms', 
                            'ec_numbers', 'gene_protein_map', 'expression', 
                            'cluster_members', 'clusters')
            """).fetchall()
            
            table_names = [t[0] for t in tables]
            
            if not all(t in table_names for t in ['sequences', 'sra_metadata']):
                self.logger.warning("Core tables don't exist. Cannot generate summary.")
                return pd.DataFrame([{
                    'status': 'error',
                    'message': 'Database schema not fully initialized'
                }])
            
            # Construct the query based on which tables exist
            query_parts = [
                "SELECT",
                "COUNT(DISTINCT s.seqhash_id) as total_sequences,",
                "COUNT(DISTINCT s.sample_id) as total_samples,"
            ]
            
            # Add counts for tables that exist
            if 'sequences' in table_names:
                query_parts.append("COUNT(DISTINCT CASE WHEN s.is_representative THEN s.seqhash_id END) as representative_sequences,")
            
            if 'annotations' in table_names:
                query_parts.append("COUNT(DISTINCT CASE WHEN a.seqhash_id IS NOT NULL THEN s.seqhash_id END) as annotated_sequences,")
            
            if 'go_terms' in table_names:
                query_parts.append("COUNT(DISTINCT CASE WHEN g.seqhash_id IS NOT NULL THEN s.seqhash_id END) as sequences_with_go,")
            
            if 'ec_numbers' in table_names:
                query_parts.append("COUNT(DISTINCT CASE WHEN e.seqhash_id IS NOT NULL THEN s.seqhash_id END) as sequences_with_ec,")
            
            if 'expression' in table_names and 'gene_protein_map' in table_names:
                query_parts.append("COUNT(DISTINCT CASE WHEN expr.gene_seqhash_id IS NOT NULL THEN gpm.protein_seqhash_id END) as sequences_with_expression,")
            
            if 'clusters' in table_names:
                query_parts.append("COUNT(DISTINCT c.cluster_id) as total_clusters,")
            
            # Add sequence statistics if sequences table exists
            if 'sequences' in table_names:
                query_parts.append("ROUND(AVG(s.length), 2) as avg_sequence_length,")
                query_parts.append("MIN(s.length) as min_sequence_length,")
                query_parts.append("MAX(s.length) as max_sequence_length")
            else:
                # Remove trailing comma from the last item if we're not adding sequence stats
                last_part = query_parts[-1]
                query_parts[-1] = last_part[:-1]  # Remove the trailing comma
            
            # Build the FROM and JOIN parts
            from_join_parts = [
                "FROM sequences s"
            ]
            
            # Add JOINs for tables that exist
            if 'annotations' in table_names:
                from_join_parts.append("LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id")
            
            if 'go_terms' in table_names:
                from_join_parts.append("LEFT JOIN go_terms g ON s.seqhash_id = g.seqhash_id")
            
            if 'ec_numbers' in table_names:
                from_join_parts.append("LEFT JOIN ec_numbers e ON s.seqhash_id = e.seqhash_id")
            
            if 'gene_protein_map' in table_names:
                from_join_parts.append("LEFT JOIN gene_protein_map gpm ON s.seqhash_id = gpm.protein_seqhash_id")
            
            if 'expression' in table_names and 'gene_protein_map' in table_names:
                from_join_parts.append("LEFT JOIN expression expr ON gpm.gene_seqhash_id = expr.gene_seqhash_id")
            
            if 'cluster_members' in table_names:
                from_join_parts.append("LEFT JOIN cluster_members cm ON s.seqhash_id = cm.seqhash_id")
            
            if 'clusters' in table_names:
                from_join_parts.append("LEFT JOIN clusters c ON cm.cluster_id = c.cluster_id")
            
            # Combine all parts
            summary_query = "\n".join(query_parts) + "\n" + "\n".join(from_join_parts)
            
            self.logger.debug(f"Summary query: {summary_query}")
            
            return self.con.execute(summary_query).df()
        
        except Exception as e:
            self.logger.error(f"Error generating database summary: {str(e)}")
            return pd.DataFrame([{
                'status': 'error',
                'message': str(e)
            }])
        
    def update_database(self, sample_ids: List[str]) -> pd.DataFrame:
        """Update existing database with new samples."""
        if not self.con.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        ).fetchall():
            raise ValueError(
                "Database not initialized. Use build_database() for new databases."
            )

        results = []
        for sample_id in sample_ids:
            result = self.process_sample(sample_id)
            results.append(result)

            if result.status == "success":
                self.logger.info(
                    f"Added {sample_id}: "
                    f"{result.sequences_loaded} sequences, "
                    f"{result.annotations_loaded} annotations"
                )
            else:
                self.logger.error(f"Failed to add {sample_id}: {result.error}")

        return pd.DataFrame([vars(r) for r in results])
