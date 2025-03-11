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
    """
    master_db_path = str(master_db_path)
    schema_sql_path = Path(schema_sql_path)

    # Read the schema SQL
    schema_sql = schema_sql_path.read_text()
    
    # Check if master database already exists
    master_exists = Path(master_db_path).exists()

    # First create the master database with the initial schema
    with duckdb.connect(master_db_path) as master_conn:
        # Set up the schema in the master database only if it doesn't exist
        if not master_exists:
            try:
                # Process the schema SQL to make it safer for execution
                safe_schema_sql = _process_schema_sql_for_safety(schema_sql)
                master_conn.execute(safe_schema_sql)
                logger.info(f"Created new master database with schema from {schema_sql_path}")
            except Exception as e:
                logger.error(f"Error creating schema: {str(e)}")
                # Create a minimal schema with essential tables if the full schema fails
                minimal_schema = _create_minimal_schema()
                master_conn.execute(minimal_schema)
                logger.info(f"Created new master database with minimal schema")
        else:
            logger.info(f"Using existing master database at {master_db_path}")

    # Now check the schema version and ensure compatibility
    if upgrade_schema:
        try:
            schema_version, was_upgraded = ensure_compatibility(
                master_db_path, required_version=target_schema_version
            )
            if was_upgraded:
                logger.info(f"Master database schema upgraded to version {schema_version}")
        except Exception as e:
            logger.error(f"Schema upgrade failed: {str(e)}")
            schema_version = 2  # Assume latest schema version if upgrade fails
    else:
        try:
            schema_version = get_db_schema_version(master_db_path)
        except Exception as e:
            logger.error(f"Could not get schema version: {str(e)}")
            schema_version = 2  # Assume latest schema version

    logger.info(f"Merging databases using schema version {schema_version}")

    # Define known tables by schema version and their dependency order (for foreign keys)
    tables_by_version = {
        1: [
            "sra_metadata",
            "sequences",          # No dependencies
            "annotations",        # Depends on sequences
            "go_terms",           # Depends on sequences
            "ec_numbers",         # Depends on sequences
            "clusters",           # Depends on sequences
            "cluster_members",    # Depends on clusters and sequences
        ],
        2: [
            "sra_metadata",       # No dependencies
            "sequences",          # No dependencies
            "annotations",        # Depends on sequences
            "go_terms",           # Depends on sequences
            "ec_numbers",         # Depends on sequences
            "clusters",           # Depends on sequences
            "cluster_members",    # Depends on clusters and sequences
            "gene_protein_map",   # Depends on sequences (must be before expression)
            "expression",         # Depends on gene_protein_map and sequences
        ],
    }

    # Get the list of tables for the current schema version
    tables = tables_by_version.get(
        schema_version, tables_by_version[max(tables_by_version.keys())]
    )
    
    # Note that kegg_info is missing from the tables list - it may not exist in the schema
    # We'll check for it dynamically

    # Now merge the source databases
    with duckdb.connect(master_db_path) as master_conn:
        # Get list of existing tables in master database
        existing_tables = set()
        try:
            result = master_conn.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
            existing_tables = {row[0].lower() for row in result}
            logger.info(f"Existing tables in master database: {', '.join(existing_tables)}")
        except Exception as e:
            logger.warning(f"Could not query existing tables: {str(e)}")
            
        # Check if kegg_info exists in the master schema
        has_kegg_info = "kegg_info" in existing_tables
        if has_kegg_info and "kegg_info" not in tables:
            # Add kegg_info to the tables list if it exists in the schema
            tables.append("kegg_info")
            
        # Process each source DuckDB
        for i, source_db in enumerate(duckdb_paths):
            source_db_str = str(source_db)
            if not Path(source_db_str).exists():
                logger.warning(f"Source database {source_db_str} does not exist. Skipping.")
                continue
                
            alias = f"db{i}"
            logger.info(f"Attaching {source_db_str} as {alias}...")

            try:
                # Check the source database schema version
                try:
                    source_schema_version = get_db_schema_version(source_db_str)
                    logger.info(f"Source database schema version: {source_schema_version}")
                except Exception as e:
                    logger.warning(f"Could not get schema version for {source_db_str}: {str(e)}")
                    source_schema_version = None

                # Attach the database
                master_conn.execute(f"ATTACH '{source_db_str}' AS {alias};")

                # Get list of tables in source database
                source_tables = set()
                try:
                    result = master_conn.execute(f"SELECT name FROM {alias}.sqlite_master WHERE type='table'").fetchall()
                    source_tables = {row[0].lower() for row in result}
                    logger.info(f"Tables in source database: {', '.join(source_tables)}")
                except Exception as e:
                    logger.warning(f"Could not query tables in source database: {str(e)}")

                # Insert data for each table, but only if it exists in both source and master
                for table in tables:
                    table_lower = table.lower()
                    
                    # Skip tables that don't exist in master
                    if table_lower not in existing_tables and len(existing_tables) > 0:
                        logger.warning(f"Table {table} doesn't exist in master database. Skipping.")
                        continue
                        
                    # Skip tables that don't exist in source
                    if table_lower not in source_tables and len(source_tables) > 0:
                        logger.warning(f"Table {table} doesn't exist in source database. Skipping.")
                        continue

                    try:
                        # Get column information for both source and target tables
                        target_columns = master_conn.execute(
                            f"PRAGMA table_info({table})"
                        ).fetchall()
                        target_col_names = [col[1] for col in target_columns]

                        source_columns = master_conn.execute(
                            f"PRAGMA table_info({alias}.{table})"
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
                            # Handle special cases for tables with common constraint violations
                            elif table == "gene_protein_map":
                                # Process gene_protein_map in batches to avoid duplicate key errors
                                try:
                                    # Get all rows from source
                                    rows = master_conn.execute(
                                        f"SELECT {col_str} FROM {alias}.{table}"
                                    ).fetchall()
                                    
                                    # Process rows one by one to avoid batch duplicates
                                    for row in rows:
                                        cols = []
                                        vals = []
                                        for j, col in enumerate(common_cols):
                                            cols.append(col)
                                            if isinstance(row[j], str):
                                                vals.append(f"'{row[j]}'")
                                            elif row[j] is None:
                                                vals.append("NULL")
                                            else:
                                                vals.append(str(row[j]))
                                                
                                        col_list = ", ".join(cols)
                                        val_list = ", ".join(vals)
                                        
                                        try:
                                            master_conn.execute(
                                                f"INSERT OR IGNORE INTO {table} ({col_list}) VALUES ({val_list})"
                                            )
                                        except Exception as e:
                                            if "duplicate key" not in str(e) and "constraint violated" not in str(e):
                                                logger.warning(f"Error inserting row into {table}: {str(e)}")
                                    
                                    logger.info(f"Processed {len(rows)} rows for {table} individually")
                                except Exception as e:
                                    logger.error(f"Error during batch processing of {table}: {str(e)}")
                            elif table == "expression":
                                # Ensure we only insert expression records with valid gene references
                                master_conn.execute(
                                    f"""
                                    INSERT OR IGNORE INTO {table} ({col_str})
                                    SELECT e.{col_str} 
                                    FROM {alias}.{table} e
                                    WHERE EXISTS (
                                        SELECT 1 FROM sequences s
                                        WHERE s.seqhash_id = e.gene_seqhash_id
                                    )
                                    """
                                )
                            else:
                                # Standard insert for other tables
                                master_conn.execute(
                                    f"""
                                    INSERT OR IGNORE INTO {table} ({col_str})
                                    SELECT {col_str}
                                    FROM {alias}.{table}
                                    """
                                )

                            # Log the number of inserted rows
                            try:
                                row_count = master_conn.execute(
                                    f"SELECT COUNT(*) FROM {table}"
                                ).fetchone()[0]
                                logger.info(f"  - Inserted rows in {table}: {row_count}")
                            except Exception as e:
                                logger.warning(f"Could not get row count for {table}: {str(e)}")
                        else:
                            logger.warning(
                                f"No common columns found for table {table}, skipping import"
                            )

                    except Exception as e:
                        logger.error(f"Error processing table {table}: {str(e)}")
                        continue

                master_conn.execute(f"DETACH {alias};")
                logger.info(f"Finished merging {source_db_str}")
            except Exception as e:
                logger.error(f"Error processing database {source_db_str}: {str(e)}")
                try:
                    master_conn.execute(f"DETACH {alias};")
                except:
                    pass
                continue

        # Optional commit; DuckDB auto-commits by default.
        master_conn.commit()

    logger.info(f"All databases have been merged into: {master_db_path}")
    return master_db_path


def _process_schema_sql_for_safety(schema_sql: str) -> str:
    """
    Process schema SQL to make it safe for repeated execution:
    - Add IF NOT EXISTS to CREATE TABLE statements
    - Remove or modify problematic statements like CREATE TABLE AS
    - Handle foreign key constraints safely
    - Skip backup table creation
    
    Args:
        schema_sql: Original schema SQL
        
    Returns:
        Modified schema SQL that's safe for repeated execution
    """
    lines = schema_sql.split('\n')
    safe_lines = []
    
    # State tracking
    in_create_table = False
    current_table = None
    skip_until_semicolon = False
    in_alter_table = False
    
    # List of tables to skip entirely
    skip_tables = ["expression_backup", "temp_", "tmp_"]
    
    # Process line by line
    for line in lines:
        # Skip blocks of code until we reach a semicolon
        if skip_until_semicolon:
            if ';' in line:
                skip_until_semicolon = False
            continue
        
        # Skip CREATE TABLE AS statements entirely
        if re.search(r'CREATE\s+TABLE\s+.*\s+AS\s+SELECT', line, re.IGNORECASE):
            skip_until_semicolon = True
            continue
            
        # Check for ALTER TABLE statements (which we want to handle carefully)
        if re.search(r'ALTER\s+TABLE', line, re.IGNORECASE):
            in_alter_table = True
            
            # Skip ALTER statements as they may be problematic
            skip_until_semicolon = True
            safe_lines.append(f"-- Skipped: {line}")
            continue
            
        # Check for CREATE TABLE statements
        create_match = re.search(r'CREATE\s+TABLE\s+(?:IF\s+NOT\s+EXISTS\s+)?(\w+)', line, re.IGNORECASE)
        if create_match:
            current_table = create_match.group(1)
            in_create_table = True
            
            # Check if this is a table we want to skip entirely
            should_skip = False
            for skip_pattern in skip_tables:
                if skip_pattern in current_table:
                    should_skip = True
                    skip_until_semicolon = True
                    break
                    
            if should_skip:
                continue
                
            # Add IF NOT EXISTS if not present
            if 'IF NOT EXISTS' not in line.upper():
                line = line.replace('CREATE TABLE', 'CREATE TABLE IF NOT EXISTS')
        
        # Add the processed line
        safe_lines.append(line)
        
        # Check if we're ending a CREATE TABLE statement
        if in_create_table and ');' in line:
            in_create_table = False
            current_table = None
            
    return '\n'.join(safe_lines)

def _create_minimal_schema() -> str:
    """
    Creates a minimal schema for the master database with the essential tables.
    Used as a fallback when the provided schema SQL fails.
    
    Returns:
        SQL string for creating a minimal schema
    """
    return """
    CREATE TABLE IF NOT EXISTS schema_version (
        version INTEGER PRIMARY KEY,
        migration_name VARCHAR,
        applied_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    );
    
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
        run_bases INTEGER,
        run_published TIMESTAMP
    );
    
    CREATE TABLE IF NOT EXISTS sequences (
        seqhash_id VARCHAR PRIMARY KEY,
        sequence VARCHAR,
        sample_id VARCHAR,
        assembly_date TIMESTAMP,
        is_representative BOOLEAN DEFAULT FALSE,
        repseq_id VARCHAR,
        length INTEGER
    );
    
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
        sample_id VARCHAR
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
        ko_number VARCHAR,
        pathway VARCHAR,
        PRIMARY KEY (seqhash_id, ko_number, pathway)
    );
    
    CREATE TABLE IF NOT EXISTS clusters (
        cluster_id VARCHAR PRIMARY KEY,
        representative_seqhash_id VARCHAR,
        size INTEGER
    );
    
    CREATE TABLE IF NOT EXISTS cluster_members (
        seqhash_id VARCHAR,
        cluster_id VARCHAR,
        PRIMARY KEY (seqhash_id, cluster_id)
    );
    
    CREATE TABLE IF NOT EXISTS gene_protein_map (
        gene_seqhash_id VARCHAR,
        protein_seqhash_id VARCHAR,
        PRIMARY KEY (gene_seqhash_id, protein_seqhash_id)
    );
    
    CREATE TABLE IF NOT EXISTS expression (
        gene_seqhash_id VARCHAR,
        sample_id VARCHAR,
        tpm DOUBLE,
        num_reads DOUBLE,
        effective_length DOUBLE,
        PRIMARY KEY (gene_seqhash_id, sample_id)
    );
    
    INSERT OR IGNORE INTO schema_version (version, migration_name, applied_at) VALUES (2, 'minimal_schema', CURRENT_TIMESTAMP);
    """


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

