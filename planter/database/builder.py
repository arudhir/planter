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
        """Drop all existing tables"""
        self.logger.info("Cleaning up any existing tables...")
        cleanup_sql = """
        DROP TABLE IF EXISTS schema_version;
        DROP TABLE IF EXISTS cluster_members;
        DROP TABLE IF EXISTS clusters;
        DROP TABLE IF EXISTS kegg_info;
        DROP TABLE IF EXISTS ec_numbers;
        DROP TABLE IF EXISTS go_terms;
        DROP TABLE IF EXISTS annotations;
        DROP TABLE IF EXISTS expression;
        DROP TABLE IF EXISTS sequences;
        DROP TABLE IF EXISTS sra_metadata;
        DROP TABLE IF EXISTS temp_annotations;
        DROP TABLE IF EXISTS temp_clusters;
        """
        with self.transaction():
            self.con.execute(cleanup_sql)

    def init_database(self):
        """Initialize database with clean schema from migration files"""
        # Clean database first
        self.clean_database()

        with self.transaction():
            # Create version tracking
            self.con.execute(
                """
                CREATE TABLE IF NOT EXISTS schema_version (
                    version INTEGER PRIMARY KEY,
                    migration_name VARCHAR NOT NULL,
                    applied_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            """
            )

            # Apply migrations
            applied = set(
                row[0]
                for row in self.con.execute(
                    """
                SELECT migration_name 
                FROM schema_version 
                ORDER BY version
            """
                ).fetchall()
            )

            for migration_file in sorted(self.schema_dir.glob("*.sql")):
                if migration_file.name not in applied:
                    self.logger.info(f"Applying migration: {migration_file.name}")
                    with open(migration_file) as f:
                        self.con.execute(f.read())

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

    def _get_sample_paths(self, sample_id: str) -> SamplePaths:
        """Get paths to sample data files."""
        sample_dir = self.output_dir / sample_id
        return SamplePaths(
            sequences=sample_dir / f"transdecoder/{sample_id}.pep",
            annotations=sample_dir / f"eggnog/{sample_id}.emapper.annotations",
            expression=sample_dir / f"quants/{sample_id}.quant.json",
        )

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
        df = pd.DataFrame(sequences)
        self.con.execute(
            "INSERT INTO sequences SELECT * FROM df"
        )  # Remove transaction here

    def _load_sequences(
        self, sequence_path: Path, sample_id: str, duplicates: list
    ) -> int:
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
                gene_seqhash_id = protein_seqhash_id.split(".p")[0]

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

        # Since gene_seqhash_id is now a primary key, we need to ensure uniqueness
        # Group by gene_seqhash_id and take the first protein for each gene
        unique_maps = {}
        for mapping in gene_protein_maps:
            gene_id = mapping["gene_seqhash_id"]
            if gene_id not in unique_maps:
                unique_maps[gene_id] = mapping

        # Convert unique maps back to a list
        unique_gene_protein_maps = list(unique_maps.values())

        # Create dataframe and insert
        df = pd.DataFrame(unique_gene_protein_maps)
        self.con.execute(
            """
            INSERT INTO gene_protein_map
            SELECT * FROM df
            ON CONFLICT DO NOTHING
        """
        )

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

        return annotations_loaded

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
            # First, try to validate the JSON file by loading it with pandas
            try:
                import json

                import pandas as pd

                with open(expression_path, "r") as f:
                    data = json.load(f)
                    self.logger.info(f"Verified JSON file format: {len(data)} records")

                # Debug log the first few entries to understand structure
                if len(data) > 0:
                    self.logger.info(f"First entry structure: {data[0].keys()}")
            except Exception as e:
                self.logger.error(f"Error validating JSON file: {e}")

            # Load JSON file and insert into expression table
            with self.transaction():
                # Try a different approach: use pandas to load the JSON, then insert to DuckDB
                import json

                import pandas as pd

                # Read JSON from file
                with open(expression_path, "r") as f:
                    json_data = json.load(f)

                # Convert to pandas DataFrame
                df = pd.DataFrame(json_data)

                # Rename columns to match our schema
                df_renamed = df.rename(
                    columns={
                        "Name": "gene_seqhash_id",  # This is the gene ID without the .p1, .p2 suffix
                        "TPM": "tpm",
                        "NumReads": "num_reads",
                        "EffectiveLength": "effective_length",
                    }
                )

                # Add sample_id column
                df_renamed["sample_id"] = sample_id

                # Register dataframe with duckdb
                self.con.register("temp_expression_df", df_renamed)

                # Insert from dataframe to expression table
                # First make sure all gene_seqhash_ids in the expression data exist in the gene_protein_map
                # If they don't exist, create new entries in gene_protein_map for them
                # This ensures referential integrity for expression data
                self.con.execute(
                    """
                    INSERT INTO gene_protein_map (gene_seqhash_id, protein_seqhash_id)
                    SELECT DISTINCT 
                        temp_expression_df.gene_seqhash_id,
                        temp_expression_df.gene_seqhash_id || '.p1' as protein_seqhash_id
                    FROM temp_expression_df
                    WHERE NOT EXISTS (
                        SELECT 1 FROM gene_protein_map
                        WHERE gene_protein_map.gene_seqhash_id = temp_expression_df.gene_seqhash_id
                    )
                    AND temp_expression_df.gene_seqhash_id NOT IN (
                        SELECT gene_seqhash_id FROM gene_protein_map
                    )
                    ON CONFLICT DO NOTHING
                """
                )

                # Then insert the expression data
                self.con.execute(
                    """
                    INSERT INTO expression (gene_seqhash_id, sample_id, tpm, num_reads, effective_length)
                    SELECT 
                        gene_seqhash_id,
                        sample_id,
                        tpm,
                        num_reads,
                        effective_length
                    FROM temp_expression_df
                    WHERE NOT EXISTS (
                        SELECT 1 FROM expression 
                        WHERE expression.gene_seqhash_id = temp_expression_df.gene_seqhash_id 
                        AND expression.sample_id = temp_expression_df.sample_id
                    )
                """
                )

            expression_loaded = self.con.execute(
                "SELECT COUNT(*) FROM expression WHERE sample_id = ?", [sample_id]
            ).fetchone()[0]

            self.logger.info(
                f"Loaded {expression_loaded} expression records for {sample_id}"
            )
            return expression_loaded

        except Exception as e:
            self.logger.error(f"Error loading expression data: {e}")
            raise

    def process_sample(self, sample_id: str) -> ProcessingResult:
        """Process a single sample: fetch metadata and load sequence data."""
        try:
            # Check if sample already exists
            if self._check_sample_exists(sample_id):
                return ProcessingResult(
                    sample_id=sample_id,
                    status="error",
                    error=f"Sample {sample_id} already exists in the database",
                )

            paths = self._get_sample_paths(sample_id)

            if not paths.sequences.exists():
                raise FileNotFoundError(f"Sequence file not found: {paths.sequences}")

            metadata = self._fetch_sra_metadata(sample_id)

            def _process():
                # Insert metadata
                if metadata:
                    self.con.execute(
                        "DELETE FROM sra_metadata WHERE sample_id = ?", [sample_id]
                    )
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

                # Load sequences - this must happen first for foreign key constraints
                duplicates = []
                sequences_loaded = self._load_sequences(
                    paths.sequences, sample_id, duplicates
                )

                # Load annotations if available
                annotations_loaded = 0
                if paths.annotations.exists():
                    annotations_loaded = self._load_annotations(
                        paths.annotations, sample_id
                    )
                else:
                    self.logger.warning(
                        f"No annotation file found: {paths.annotations}"
                    )

                # Check if the sequences were properly loaded
                seq_count = self.con.execute(
                    "SELECT COUNT(*) FROM sequences WHERE sample_id = ?", [sample_id]
                ).fetchone()[0]
                self.logger.info(
                    f"Verified sequence loading: {seq_count} sequences found in DB"
                )

                # Load expression data if available - only after sequences are loaded
                expression_loaded = 0
                if paths.expression.exists():
                    # Make sure we have sequences to link expression data to
                    if seq_count > 0:
                        expression_loaded = self._load_expression(
                            paths.expression, sample_id
                        )
                    else:
                        self.logger.warning(
                            "Cannot load expression data: no sequences loaded"
                        )
                else:
                    self.logger.warning(f"No expression file found: {paths.expression}")

                return (
                    sequences_loaded,
                    annotations_loaded,
                    expression_loaded,
                    duplicates,
                )

            # Execute everything in a single transaction
            sequences_loaded, annotations_loaded, expression_loaded, duplicates = (
                self.execute_transaction(_process)
            )

            return ProcessingResult(
                sample_id=sample_id,
                status="success",
                sequences_loaded=sequences_loaded,
                annotations_loaded=annotations_loaded,
                duplicates=len(duplicates),
            )

        except Exception as e:
            self.logger.error(f"Error processing {sample_id}: {str(e)}")
            return ProcessingResult(sample_id=sample_id, status="error", error=str(e))

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
        summary_query = """
        SELECT 
            COUNT(DISTINCT s.seqhash_id) as total_sequences,
            COUNT(DISTINCT s.sample_id) as total_samples,
            COUNT(DISTINCT CASE WHEN s.is_representative THEN s.seqhash_id END) as representative_sequences,
            COUNT(DISTINCT CASE WHEN a.seqhash_id IS NOT NULL THEN s.seqhash_id END) as annotated_sequences,
            COUNT(DISTINCT CASE WHEN g.seqhash_id IS NOT NULL THEN s.seqhash_id END) as sequences_with_go,
            COUNT(DISTINCT CASE WHEN e.seqhash_id IS NOT NULL THEN s.seqhash_id END) as sequences_with_ec,
            COUNT(DISTINCT CASE WHEN expr.gene_seqhash_id IS NOT NULL THEN s.seqhash_id END) as sequences_with_expression,
            COUNT(DISTINCT c.cluster_id) as total_clusters,
            ROUND(AVG(s.length), 2) as avg_sequence_length,
            MIN(s.length) as min_sequence_length,
            MAX(s.length) as max_sequence_length
        FROM sequences s
        LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
        LEFT JOIN go_terms g ON s.seqhash_id = g.seqhash_id
        LEFT JOIN ec_numbers e ON s.seqhash_id = e.seqhash_id
        LEFT JOIN gene_protein_map gpm ON s.seqhash_id = gpm.protein_seqhash_id
        LEFT JOIN expression expr ON gpm.gene_seqhash_id = expr.gene_seqhash_id
        LEFT JOIN cluster_members cm ON s.seqhash_id = cm.seqhash_id
        LEFT JOIN clusters c ON cm.cluster_id = c.cluster_id
        """
        return self.con.execute(summary_query).df()

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
