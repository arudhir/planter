#!/usr/bin/env python3
"""
Clustering Manager - handles database operations for clustering.

This module manages the storage and retrieval of clustering results,
keeping the pure clustering logic separate from database concerns.
"""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Optional, Tuple

import duckdb
import pandas as pd

from .cluster import ClusteringParams, ClusteringResult, run_clustering

logger = logging.getLogger(__name__)


class ClusteringManager:
    """
    Manages clustering operations and their storage in the database.

    This class handles:
    - Creating and tracking clustering runs
    - Storing immutable clustering results
    - Retrieving current/historical clustering data
    - Exporting sequences for clustering
    """

    def __init__(self, db_path: str):
        """Initialize with path to DuckDB database."""
        self.db_path = str(db_path)

    def _ensure_schema(self, con: duckdb.DuckDBPyConnection) -> None:
        """Ensure the clustering schema exists."""
        schema_path = Path(__file__).parent.parent / "database" / "schema" / "migrations" / "005_immutable_clustering.sql"

        # First check if tables already exist
        existing_tables = {row[0] for row in con.execute(
            "SELECT table_name FROM information_schema.tables WHERE table_schema = 'main'"
        ).fetchall()}

        if "clustering_runs" in existing_tables:
            # Schema already applied
            return

        if schema_path.exists():
            schema_sql = schema_path.read_text()
            # Try to execute migration statements
            migration_failed = False
            for statement in schema_sql.split(";"):
                statement = statement.strip()
                if statement and not statement.startswith("--"):
                    try:
                        con.execute(statement)
                    except Exception as e:
                        error_msg = str(e).lower()
                        # Ignore "already exists" errors
                        if "already exists" in error_msg:
                            continue
                        # FK constraint failures or table-not-found errors mean we should use fallback
                        if any(x in error_msg for x in ["foreign key", "not found", "does not exist"]):
                            logger.debug(f"Migration statement failed (using fallback): {e}")
                            migration_failed = True
                            break
                        else:
                            # Other errors should be raised
                            raise

            if migration_failed:
                # Use fallback without FK constraints
                self._create_clustering_tables(con)
        else:
            # Fallback: create tables directly
            self._create_clustering_tables(con)

    def _create_clustering_tables(self, con: duckdb.DuckDBPyConnection) -> None:
        """Create clustering tables if they don't exist."""
        con.execute("""
            CREATE TABLE IF NOT EXISTS clustering_runs (
                run_id INTEGER PRIMARY KEY,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                completed_at TIMESTAMP,
                status VARCHAR DEFAULT 'pending',
                parameters VARCHAR,
                sequence_count INTEGER,
                cluster_count INTEGER,
                notes VARCHAR
            )
        """)

        con.execute("""
            CREATE TABLE IF NOT EXISTS versioned_clusters (
                run_id INTEGER NOT NULL,
                cluster_id VARCHAR NOT NULL,
                representative_seqhash_id VARCHAR NOT NULL,
                size INTEGER NOT NULL,
                PRIMARY KEY (run_id, cluster_id)
            )
        """)

        con.execute("""
            CREATE TABLE IF NOT EXISTS versioned_cluster_members (
                run_id INTEGER NOT NULL,
                seqhash_id VARCHAR NOT NULL,
                cluster_id VARCHAR NOT NULL,
                PRIMARY KEY (run_id, seqhash_id)
            )
        """)

    def create_run(
        self,
        params: Optional[ClusteringParams] = None,
        notes: Optional[str] = None
    ) -> int:
        """
        Create a new clustering run entry.

        Returns:
            run_id for the new clustering run
        """
        with duckdb.connect(self.db_path) as con:
            self._ensure_schema(con)

            # Get next run_id
            result = con.execute("SELECT COALESCE(MAX(run_id), 0) + 1 FROM clustering_runs").fetchone()
            run_id = result[0]

            params_json = json.dumps(params.to_dict()) if params else None

            con.execute("""
                INSERT INTO clustering_runs (run_id, status, parameters, notes)
                VALUES (?, 'pending', ?, ?)
            """, [run_id, params_json, notes])

            logger.info(f"Created clustering run {run_id}")
            return run_id

    def update_run_status(
        self,
        run_id: int,
        status: str,
        sequence_count: Optional[int] = None,
        cluster_count: Optional[int] = None
    ) -> None:
        """Update the status of a clustering run."""
        with duckdb.connect(self.db_path) as con:
            if status == "completed":
                con.execute("""
                    UPDATE clustering_runs
                    SET status = ?,
                        completed_at = CURRENT_TIMESTAMP,
                        sequence_count = ?,
                        cluster_count = ?
                    WHERE run_id = ?
                """, [status, sequence_count, cluster_count, run_id])
            else:
                con.execute("""
                    UPDATE clustering_runs
                    SET status = ?
                    WHERE run_id = ?
                """, [status, run_id])

    def store_results(self, run_id: int, result: ClusteringResult) -> None:
        """
        Store clustering results in the database.

        This is an atomic operation - either all results are stored or none.
        """
        if not result.success:
            raise ValueError(f"Cannot store failed clustering result: {result.error_message}")

        with duckdb.connect(self.db_path) as con:
            self._ensure_schema(con)

            try:
                con.execute("BEGIN TRANSACTION")

                # Update run status to running
                con.execute("""
                    UPDATE clustering_runs SET status = 'running' WHERE run_id = ?
                """, [run_id])

                # Calculate cluster sizes
                cluster_sizes = result.assignments.groupby("cluster_id").size().reset_index(name="size")
                cluster_reps = result.assignments[result.assignments["is_representative"]][
                    ["cluster_id", "representative_seqhash_id"]
                ].drop_duplicates()

                clusters_df = cluster_reps.merge(cluster_sizes, on="cluster_id")
                clusters_df["run_id"] = run_id

                # Insert clusters
                con.register("clusters_to_insert", clusters_df)
                con.execute("""
                    INSERT INTO versioned_clusters (run_id, cluster_id, representative_seqhash_id, size)
                    SELECT run_id, cluster_id, representative_seqhash_id, size
                    FROM clusters_to_insert
                """)

                # Insert cluster members
                members_df = result.assignments[["seqhash_id", "cluster_id"]].copy()
                members_df["run_id"] = run_id

                con.register("members_to_insert", members_df)
                con.execute("""
                    INSERT INTO versioned_cluster_members (run_id, seqhash_id, cluster_id)
                    SELECT run_id, seqhash_id, cluster_id
                    FROM members_to_insert
                """)

                # Update run as completed
                con.execute("""
                    UPDATE clustering_runs
                    SET status = 'completed',
                        completed_at = CURRENT_TIMESTAMP,
                        sequence_count = ?,
                        cluster_count = ?
                    WHERE run_id = ?
                """, [result.sequence_count, result.cluster_count, run_id])

                con.execute("COMMIT")
                logger.info(f"Stored clustering results for run {run_id}: "
                           f"{result.sequence_count} sequences, {result.cluster_count} clusters")

            except Exception as e:
                con.execute("ROLLBACK")
                # Mark run as failed
                con.execute("""
                    UPDATE clustering_runs SET status = 'failed' WHERE run_id = ?
                """, [run_id])
                raise RuntimeError(f"Failed to store clustering results: {e}")

    def export_sequences_for_clustering(self, output_path: Path) -> int:
        """
        Export all sequences to a FASTA file for clustering.

        Returns:
            Number of sequences exported
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with duckdb.connect(self.db_path) as con:
            sequences = con.execute("""
                SELECT seqhash_id, sequence
                FROM sequences
                ORDER BY seqhash_id
            """).fetchdf()

        count = 0
        with open(output_path, "w") as f:
            for _, row in sequences.iterrows():
                f.write(f">{row['seqhash_id']}\n")
                f.write(f"{row['sequence']}\n")
                count += 1

        logger.info(f"Exported {count} sequences to {output_path}")
        return count

    def get_current_representatives(self) -> pd.DataFrame:
        """
        Get representative sequences from the latest completed clustering run.

        Returns:
            DataFrame with seqhash_id, sequence for representatives
        """
        with duckdb.connect(self.db_path) as con:
            self._ensure_schema(con)

            # Check if we have any completed runs
            has_runs = con.execute("""
                SELECT COUNT(*) FROM clustering_runs WHERE status = 'completed'
            """).fetchone()[0]

            if has_runs == 0:
                logger.warning("No completed clustering runs found")
                return pd.DataFrame(columns=["seqhash_id", "sequence"])

            return con.execute("""
                SELECT s.seqhash_id, s.sequence
                FROM versioned_clusters vc
                JOIN sequences s ON vc.representative_seqhash_id = s.seqhash_id
                WHERE vc.run_id = (
                    SELECT MAX(run_id)
                    FROM clustering_runs
                    WHERE status = 'completed'
                )
            """).fetchdf()

    def export_representatives_fasta(self, output_path: Path) -> int:
        """
        Export representative sequences to FASTA file.

        Returns:
            Number of sequences exported
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        reps = self.get_current_representatives()

        if reps.empty:
            # No previous clustering - return empty file
            output_path.touch()
            return 0

        count = 0
        with open(output_path, "w") as f:
            for _, row in reps.iterrows():
                f.write(f">{row['seqhash_id']}\n")
                f.write(f"{row['sequence']}\n")
                count += 1

        logger.info(f"Exported {count} representative sequences to {output_path}")
        return count

    def get_current_clustering(self) -> pd.DataFrame:
        """
        Get the current (latest completed) clustering assignments.

        Returns:
            DataFrame with seqhash_id, cluster_id, representative_seqhash_id, cluster_size
        """
        with duckdb.connect(self.db_path) as con:
            self._ensure_schema(con)

            return con.execute("""
                SELECT
                    vcm.seqhash_id,
                    vcm.cluster_id,
                    vc.representative_seqhash_id,
                    vc.size as cluster_size
                FROM versioned_cluster_members vcm
                JOIN versioned_clusters vc
                    ON vcm.run_id = vc.run_id AND vcm.cluster_id = vc.cluster_id
                WHERE vcm.run_id = (
                    SELECT MAX(run_id)
                    FROM clustering_runs
                    WHERE status = 'completed'
                )
            """).fetchdf()

    def get_clustering_history(self) -> pd.DataFrame:
        """Get history of all clustering runs."""
        with duckdb.connect(self.db_path) as con:
            self._ensure_schema(con)

            return con.execute("""
                SELECT
                    run_id,
                    created_at,
                    completed_at,
                    status,
                    sequence_count,
                    cluster_count,
                    notes
                FROM clustering_runs
                ORDER BY run_id DESC
            """).fetchdf()

    def run_full_clustering(
        self,
        params: Optional[ClusteringParams] = None,
        work_dir: Optional[Path] = None,
        notes: Optional[str] = None
    ) -> Tuple[int, ClusteringResult]:
        """
        Run a complete clustering workflow.

        This is the main entry point for clustering:
        1. Creates a new run
        2. Exports all sequences
        3. Runs clustering
        4. Stores results

        Returns:
            Tuple of (run_id, ClusteringResult)
        """
        if params is None:
            params = ClusteringParams()

        # Create run entry
        run_id = self.create_run(params, notes)

        try:
            # Update status to running
            self.update_run_status(run_id, "running")

            # Create temp directory for clustering
            if work_dir is None:
                import tempfile
                temp_dir = tempfile.TemporaryDirectory()
                work_dir = Path(temp_dir.name)
            else:
                work_dir = Path(work_dir)
                work_dir.mkdir(parents=True, exist_ok=True)
                temp_dir = None

            try:
                # Export sequences
                fasta_path = work_dir / "all_sequences.fasta"
                seq_count = self.export_sequences_for_clustering(fasta_path)

                if seq_count == 0:
                    raise ValueError("No sequences to cluster")

                # Run clustering (pure function)
                result = run_clustering(
                    fasta_path=fasta_path,
                    params=params,
                    work_dir=work_dir
                )

                if not result.success:
                    self.update_run_status(run_id, "failed")
                    return run_id, result

                # Store results
                self.store_results(run_id, result)

                return run_id, result

            finally:
                if temp_dir:
                    temp_dir.cleanup()

        except Exception as e:
            logger.error(f"Clustering run {run_id} failed: {e}")
            self.update_run_status(run_id, "failed")
            raise

    def get_cluster_stats(self) -> pd.DataFrame:
        """Get statistics about the current clustering."""
        with duckdb.connect(self.db_path) as con:
            self._ensure_schema(con)

            return con.execute("""
                SELECT
                    COUNT(DISTINCT cluster_id) as total_clusters,
                    COUNT(DISTINCT seqhash_id) as total_sequences,
                    AVG(cluster_size) as avg_cluster_size,
                    MIN(cluster_size) as min_cluster_size,
                    MAX(cluster_size) as max_cluster_size,
                    MEDIAN(cluster_size) as median_cluster_size
                FROM (
                    SELECT
                        vcm.seqhash_id,
                        vcm.cluster_id,
                        vc.size as cluster_size
                    FROM versioned_cluster_members vcm
                    JOIN versioned_clusters vc
                        ON vcm.run_id = vc.run_id AND vcm.cluster_id = vc.cluster_id
                    WHERE vcm.run_id = (
                        SELECT MAX(run_id)
                        FROM clustering_runs
                        WHERE status = 'completed'
                    )
                )
            """).fetchdf()
