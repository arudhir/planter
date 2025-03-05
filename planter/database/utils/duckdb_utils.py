#!/usr/bin/env python3
"""
DuckDB utility functions for merging and updating databases.

This module provides functions for working with DuckDB databases,
including merging multiple databases and updating cluster information.
It includes schema compatibility features to ensure forward and backward
compatibility between different schema versions.
"""
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import duckdb
import pandas as pd

from planter.database.schema.schema_version import (SCHEMA_VERSIONS,
                                                    ensure_compatibility,
                                                    get_db_schema_version)

logger = logging.getLogger(__name__)


def merge_duckdbs(
    duckdb_paths: List[Union[str, Path]],
    master_db_path: Union[str, Path],
    schema_sql_path: Union[str, Path],
    upgrade_schema: bool = True,
    target_schema_version: Optional[int] = None,
) -> str:
    """
    Merge multiple DuckDB databases into a master DuckDB with schema compatibility.

    Args:
        duckdb_paths: List of paths to source DuckDB files.
        master_db_path: Path to the master (merged) DuckDB.
        schema_sql_path: Path to the SQL file defining the schema.
        upgrade_schema: Whether to automatically upgrade schema (default: True).
        target_schema_version: Target schema version (default: latest version).

    Returns:
        Path to the master (merged) DuckDB.

    The function:
      - Creates (or opens) the master database.
      - Executes the schema SQL to create tables if they don't exist.
      - Ensures schema compatibility across different versions.
      - Iterates through each source database, attaches it,
        and inserts data into the master tables in dependency order.
      - Uses INSERT OR IGNORE to avoid duplicate key errors.
      - Detaches each source database after merging.
    """
    master_db_path = str(master_db_path)
    schema_sql_path = Path(schema_sql_path)

    # Read the schema SQL
    schema_sql = schema_sql_path.read_text()

    # First create the master database with the initial schema
    with duckdb.connect(master_db_path) as master_conn:
        # Set up the schema in the master database
        master_conn.execute(schema_sql)

    # Now check the schema version and ensure compatibility
    if upgrade_schema:
        schema_version, was_upgraded = ensure_compatibility(
            master_db_path, required_version=target_schema_version
        )
        if was_upgraded:
            logger.info(f"Master database schema upgraded to version {schema_version}")
    else:
        schema_version = get_db_schema_version(master_db_path)

    logger.info(f"Merging databases using schema version {schema_version}")

    # Define known tables by schema version
    tables_by_version = {
        1: [
            "sra_metadata",
            "sequences",
            "annotations",
            "go_terms",
            "ec_numbers",
            "kegg_info",
            "clusters",
            "cluster_members",
        ],
        2: [
            "sra_metadata",
            "sequences",
            "annotations",
            "go_terms",
            "ec_numbers",
            "kegg_info",
            "clusters",
            "cluster_members",
            "expression",
            "gene_protein_map",
        ],
    }

    # Get the list of tables for the current schema version
    tables = tables_by_version.get(
        schema_version, tables_by_version[max(tables_by_version.keys())]
    )

    # Now merge the source databases
    with duckdb.connect(master_db_path) as master_conn:
        # Process each source DuckDB
        for i, source_db in enumerate(duckdb_paths):
            alias = f"db{i}"
            source_db_str = str(source_db)
            logger.info(f"Attaching {source_db_str} as {alias}...")

            # Check the source database schema version
            source_schema_version = get_db_schema_version(source_db_str)
            logger.info(f"Source database schema version: {source_schema_version}")

            master_conn.execute(f"ATTACH '{source_db_str}' AS {alias};")

            # Insert data for each table, but only if it exists in the source
            for table in tables:
                # Check if the table exists in the source database
                table_exists_in_source = master_conn.execute(
                    f"""
                    SELECT COUNT(*) FROM {alias}.sqlite_master 
                    WHERE type='table' AND name='{table}'
                """
                ).fetchone()[0]

                if table_exists_in_source:
                    # Get column information for both source and target tables
                    target_columns = master_conn.execute(
                        f"""
                        PRAGMA table_info({table})
                    """
                    ).fetchall()
                    target_col_names = [col[1] for col in target_columns]

                    source_columns = master_conn.execute(
                        f"""
                        PRAGMA table_info({alias}.{table})
                    """
                    ).fetchall()
                    source_col_names = [col[1] for col in source_columns]

                    # Find common columns for safe copy
                    common_cols = [
                        col for col in target_col_names if col in source_col_names
                    ]

                    if common_cols:
                        # Log the column adaptation
                        if set(target_col_names) != set(source_col_names):
                            logger.info(
                                f"Adapting {table} schema: matching {len(common_cols)} columns"
                            )
                            logger.info(f"Target columns: {target_col_names}")
                            logger.info(f"Source columns: {source_col_names}")

                        # Build a dynamic query with only common columns
                        col_str = ", ".join(common_cols)

                        # Special case for sequences table in v2 schema (add is_representative)
                        if (
                            table == "sequences"
                            and "is_representative" in target_col_names
                            and "is_representative" not in source_col_names
                        ):
                            # Add is_representative column with FALSE default
                            all_cols = list(common_cols)
                            all_cols.append("is_representative")
                            target_col_str = ", ".join(all_cols)

                            source_col_str = ", ".join(common_cols)
                            source_col_str += ", FALSE as is_representative"

                            master_conn.execute(
                                f"""
                                INSERT OR IGNORE INTO {table} ({target_col_str})
                                SELECT {source_col_str}
                                FROM {alias}.{table}
                            """
                            )
                        else:
                            # Standard copy of matching columns
                            master_conn.execute(
                                f"""
                                INSERT OR IGNORE INTO {table} ({col_str})
                                SELECT {col_str}
                                FROM {alias}.{table}
                            """
                            )
                    else:
                        logger.warning(
                            f"No common columns found for table {table}, skipping import"
                        )

                    # Log the number of inserted rows
                    row_count = master_conn.execute(
                        f"""
                        SELECT COUNT(*) FROM {table}
                    """
                    ).fetchone()[0]
                    logger.info(f"  - Inserted rows in {table}: {row_count}")

            master_conn.execute(f"DETACH {alias};")
            logger.info(f"Finished merging {source_db_str}")

        # Optional commit; DuckDB auto-commits by default.
        master_conn.commit()

    logger.info(f"All databases have been merged into: {master_db_path}")
    return master_db_path


def update_duckdb_with_cluster_info(
    db_path: Union[str, Path], tsv_path: Union[str, Path], upgrade_schema: bool = True
) -> None:
    """
    Updates a DuckDB database with clustering information from an MMSeqs2 TSV file.
    Handles schema compatibility across different versions.

    Args:
        db_path: Path to the DuckDB database file.
        tsv_path: Path to the MMSeqs2 TSV file containing clustering data.
        upgrade_schema: Whether to automatically upgrade the schema if needed (default: True)
    """
    db_path = str(db_path)
    tsv_path = str(tsv_path)

    # Check and ensure schema compatibility
    schema_version, was_upgraded = (
        ensure_compatibility(db_path)
        if upgrade_schema
        else (get_db_schema_version(db_path), False)
    )
    if was_upgraded:
        logger.info(
            f"Database schema was automatically upgraded to version {schema_version}"
        )

    # Load the MMSeqs2 TSV into a DataFrame
    df = pd.read_csv(
        tsv_path, sep="\t", names=["representative_seqhash_id", "seqhash_id"]
    )
    df = df.drop_duplicates()

    # Connect to DuckDB
    con = duckdb.connect(db_path)

    try:
        logger.info(
            f"Updating database '{db_path}' with MMSeqs2 clustering information from '{tsv_path}'"
        )
        logger.info(f"Database schema version: {schema_version}")
        logger.info(f"Found {len(df)} cluster assignments")

        # Begin transaction
        con.execute("BEGIN TRANSACTION")

        # Create a temporary table for the MMSeqs2 data
        con.execute(
            """
            CREATE TEMPORARY TABLE mmseqs2_clusters (
                representative_seqhash_id VARCHAR,
                seqhash_id VARCHAR
            )
        """
        )

        # Insert data into the temporary table
        con.executemany(
            "INSERT INTO mmseqs2_clusters VALUES (?, ?)", df.values.tolist()
        )

        # First, verify that all seqhash_ids exist in the sequences table
        missing_seqs = con.execute(
            """
            SELECT DISTINCT mc.seqhash_id 
            FROM mmseqs2_clusters mc
            LEFT JOIN sequences s ON mc.seqhash_id = s.seqhash_id
            WHERE s.seqhash_id IS NULL
        """
        ).fetchall()

        if missing_seqs:
            missing_count = len(missing_seqs)
            logger.warning(
                f"{missing_count} sequence IDs from cluster file not found in the database"
            )
            if missing_count > 10:
                logger.warning(
                    f"First 10 missing IDs: {[m[0] for m in missing_seqs[:10]]}"
                )
            else:
                logger.warning(f"Missing IDs: {[m[0] for m in missing_seqs]}")

        # Handle different schema versions
        if schema_version >= 2:
            # Schema version 2+ has is_representative column
            logger.info("Using schema v2+ with is_representative column")

            # Mark representative sequences in the sequences table
            con.execute(
                """
                UPDATE sequences 
                SET is_representative = FALSE
            """
            )

            con.execute(
                """
                UPDATE sequences 
                SET is_representative = TRUE
                WHERE seqhash_id IN (
                    SELECT DISTINCT representative_seqhash_id FROM mmseqs2_clusters
                )
            """
            )
        else:
            # Schema version 1 doesn't have is_representative column
            logger.info("Using schema v1 without is_representative column")

        # Update the sequences table to assign the representative sequence - works in all versions
        con.execute(
            """
            UPDATE sequences 
            SET repseq_id = mm.representative_seqhash_id
            FROM mmseqs2_clusters mm
            WHERE sequences.seqhash_id = mm.seqhash_id
        """
        )

        # Get the count of updated sequences for reporting
        updated_count = con.execute(
            """
            SELECT COUNT(DISTINCT seqhash_id) FROM mmseqs2_clusters
        """
        ).fetchone()[0]

        # Clear existing cluster information to rebuild it
        con.execute("DELETE FROM cluster_members")
        con.execute("DELETE FROM clusters")

        # Insert new clusters into the clusters table
        # Use proper cluster_id format: 'CLUSTER_nnn'
        con.execute(
            """
            INSERT INTO clusters (cluster_id, representative_seqhash_id, size)
            SELECT 
                'CLUSTER_' || ROW_NUMBER() OVER (ORDER BY representative_seqhash_id) as cluster_id,
                representative_seqhash_id, 
                COUNT(*) as size
            FROM mmseqs2_clusters
            GROUP BY representative_seqhash_id
        """
        )

        # Insert cluster memberships into the cluster_members table
        # Use ON CONFLICT to handle potential duplicates
        con.execute(
            """
            INSERT INTO cluster_members (seqhash_id, cluster_id)
            SELECT m.seqhash_id, c.cluster_id 
            FROM mmseqs2_clusters m
            JOIN clusters c ON m.representative_seqhash_id = c.representative_seqhash_id
            ON CONFLICT (seqhash_id) DO UPDATE SET cluster_id = EXCLUDED.cluster_id
        """
        )

        # Get cluster statistics for reporting
        cluster_stats = con.execute(
            """
            SELECT COUNT(*) as total_clusters, AVG(size) as avg_size, MAX(size) as max_size 
            FROM clusters
        """
        ).fetchone()

        # Commit transaction
        con.execute("COMMIT")

        # Report success statistics
        logger.info(
            f"Successfully updated {updated_count} sequences with new cluster assignments"
        )
        logger.info(
            f"Created {cluster_stats[0]} clusters with average size {cluster_stats[1]:.1f} (max: {cluster_stats[2]})"
        )

    except Exception as e:
        # Rollback on error
        con.execute("ROLLBACK")
        logger.error(f"Error updating database: {str(e)}")
        raise
    finally:
        # Drop the temporary table and close connection
        con.execute("DROP TABLE IF EXISTS mmseqs2_clusters")
        con.close()

    logger.info(
        f"Database '{db_path}' successfully updated with MMSeqs2 clustering information."
    )
