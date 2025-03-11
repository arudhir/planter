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
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from planter.database.schema.schema_version import (SCHEMA_VERSIONS,
                                                    ensure_compatibility,
                                                    get_db_schema_version)
from planter.database.builder import SequenceDBBuilder
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
    db_path: Union[str, Path], 
    tsv_path: Union[str, Path],
    upgrade_schema: bool = True
) -> None:
    """
    Updates cluster information in a DuckDB database, working within the 
    constraints of DuckDB's foreign key implementation.
    """
    db_path = str(db_path)
    tsv_path = str(tsv_path)
    
    con = duckdb.connect(db_path)
    
    try:
        con.execute("BEGIN TRANSACTION")
        
        # 1. Load the cluster file into a temp table
        con.execute(f"""
            CREATE TEMP TABLE new_clustering AS
            SELECT 
                column0 AS rep_id,
                column1 AS seq_id
            FROM read_csv_auto('{tsv_path}', sep='\t', header=FALSE)
        """)
        
        # 2. Update sequences table with representative info
        con.execute("""
            UPDATE sequences
            SET repseq_id = nc.rep_id,
                is_representative = (sequences.seqhash_id = nc.rep_id)
            FROM new_clustering nc
            WHERE sequences.seqhash_id = nc.seq_id
        """)
        
        # 3. Create mapping between sequences and new cluster IDs
        con.execute("""
            CREATE TEMP TABLE cluster_assignments AS
            SELECT 
                rep_id,
                'CLUSTER_' || ROW_NUMBER() OVER (ORDER BY rep_id) AS cluster_id,
                COUNT(*) OVER (PARTITION BY rep_id) AS size
            FROM new_clustering
            GROUP BY rep_id
        """)
        
        # 4. Update cluster_members (handle foreign key constraints)
        # First, clear cluster_members for each sequence that will be reassigned
        # (This avoids violating the PK constraint on seqhash_id)
        con.execute("""
            DELETE FROM cluster_members
            WHERE seqhash_id IN (SELECT seq_id FROM new_clustering)
        """)
        
        # Now safely insert new cluster memberships
        con.execute("""
            INSERT INTO cluster_members (seqhash_id, cluster_id)
            SELECT nc.seq_id, ca.cluster_id
            FROM new_clustering nc
            JOIN cluster_assignments ca ON nc.rep_id = ca.rep_id
        """)
        
        # 5. Update clusters table
        # First, identify clusters that aren't used anymore
        con.execute("""
            CREATE TEMP TABLE unused_clusters AS
            SELECT cluster_id 
            FROM clusters
            WHERE cluster_id NOT IN (SELECT cluster_id FROM cluster_members)
        """)
        
        # Delete unused clusters (safe because no members reference them)
        con.execute("""
            DELETE FROM clusters
            WHERE cluster_id IN (SELECT cluster_id FROM unused_clusters)
        """)
        
        # Now insert new clusters
        con.execute("""
            INSERT INTO clusters (cluster_id, representative_seqhash_id, size)
            SELECT 
                ca.cluster_id,
                ca.rep_id,
                ca.size
            FROM cluster_assignments ca
            WHERE NOT EXISTS (
                SELECT 1 FROM clusters c WHERE c.representative_seqhash_id = ca.rep_id
            )
        """)
        
        # Update sizes for existing clusters
        con.execute("""
            UPDATE clusters
            SET size = ca.size
            FROM cluster_assignments ca
            WHERE clusters.representative_seqhash_id = ca.rep_id
        """)
        
        con.execute("COMMIT")
        
    except Exception as e:
        con.execute("ROLLBACK")
        raise e
    finally:
        # Clean up temp tables
        for table in ["new_clustering", "cluster_assignments", "unused_clusters"]:
            try:
                con.execute(f"DROP TABLE IF EXISTS {table}")
            except:
                pass
        con.close()
        
def extract_representative_sequences(
    db_path: Union[str, Path],
    output_path: Union[str, Path]
) -> None:
    """
    Extracts representative sequences from a DuckDB database and saves them to a FASTA file.

    Args:
        db_path (Union[str, Path]): Path to the DuckDB database.
        output_path (Union[str, Path]): Path to save the extracted representative sequences in FASTA format.

    Returns:
        None
    """
    # Ensure db_path and output_path are strings before converting to Path
    db_path = Path(str(db_path)).resolve()
    output_path = Path(str(output_path)).resolve()

    logger.info(f"Extracting representative sequences from {db_path}")

    query = """
        SELECT seqhash_id, sequence
        FROM sequences
        WHERE repseq_id = seqhash_id;
    """

    try:
        # Connect to DuckDB
        con = duckdb.connect(str(db_path))

        # Execute query and fetch results
        results = con.execute(query).fetchall()

        # Convert results to SeqRecords
        seq_records = [
            SeqRecord(Seq(seq), id=seq_id, description="")
            for seq_id, seq in results
        ]

        # Write to FASTA file using SeqIO
        with output_path.open("w") as fasta_file:
            SeqIO.write(seq_records, fasta_file, "fasta")

        logger.info(f"Successfully extracted {len(seq_records)} sequences to {output_path}")

    except Exception as e:
        logger.error(f"Error extracting sequences: {e}")

    finally:
        con.close()
    return output_path

def create_duckdb(sample_id: str, outdir: str, duckdb_out: str):
    """Create DuckDB database for the given sample ID."""

    logger.info(f"Processing sample {sample_id}")
    with SequenceDBBuilder(duckdb_out, output_dir=outdir) as builder:
        results = builder.build_database([sample_id])
        logger.info(f"Build results: {results}")

        summary = builder.get_database_summary()
        logger.info(f"Database summary: {summary}")
