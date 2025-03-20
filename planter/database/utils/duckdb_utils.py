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
import tempfile
import shutil
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Union
import os
import re
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
            # Master DB should exist but handle the case if it doesn't
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
                                # First make sure gene_protein_map entries exist for these genes
                                try:
                                    # Get all genes from the expression table that need to be in gene_protein_map
                                    gene_ids = master_conn.execute(
                                        f"""
                                        SELECT DISTINCT gene_seqhash_id 
                                        FROM {alias}.{table} 
                                        WHERE gene_seqhash_id NOT IN (
                                            SELECT gene_seqhash_id FROM gene_protein_map
                                        )
                                        """
                                    ).fetchall()
                                    
                                    # For each missing gene, try to find a matching protein and add an entry
                                    for (gene_id,) in gene_ids:
                                        # Try to find a protein with matching ID pattern (removing .p1, etc.)
                                        base_id = gene_id.split('.')[0] if '.' in gene_id else gene_id
                                        protein_query = f"{base_id}.p1"
                                        
                                        # Add gene-protein mapping if the protein exists
                                        master_conn.execute(
                                            f"""
                                            INSERT OR IGNORE INTO gene_protein_map (gene_seqhash_id, protein_seqhash_id)
                                            VALUES ('{gene_id}', '{protein_query}')
                                            """
                                        )
                                except Exception as e:
                                    logger.warning(f"Error ensuring gene_protein_map entries: {str(e)}")
                                
                                # Now insert expression data that has valid gene mappings
                                master_conn.execute(
                                    f"""
                                    INSERT OR IGNORE INTO {table} ({col_str})
                                    SELECT e.{col_str} 
                                    FROM {alias}.{table} e
                                    WHERE EXISTS (
                                        SELECT 1 FROM gene_protein_map gpm
                                        WHERE gpm.gene_seqhash_id = e.gene_seqhash_id
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

def update_clusters(
    db_path: Union[str, Path], 
    tsv_path: Union[str, Path],
    backup_first: bool = True,
    log_path: Union[str, Path] = None,
    handle_duplicates: str = "replace"  # "replace", "ignore", or "error"
) -> None:
    """
    Update cluster information, skipping missing sequences and logging them.
    
    Args:
        db_path: Path to the DuckDB database file
        tsv_path: Path to the TSV file containing cluster information
        backup_first: Whether to make a backup of the database first
        log_path: Path to save the log file (default: generated based on db_path)
        handle_duplicates: How to handle duplicate keys ("replace", "ignore", or "error")
    """
    db_path = Path(str(db_path))
    tsv_path = Path(str(tsv_path))
    
    # Validate handle_duplicates parameter
    if handle_duplicates not in ["replace", "ignore", "error"]:
        raise ValueError(f"Invalid value for handle_duplicates: {handle_duplicates}")
    
    # Set up logging
    if log_path is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_path = db_path.parent / f"cluster_update_{timestamp}.log"
    
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )
    
    # Create console handler to display logs as well
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter("%(message)s")
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)
    
    logging.info(f"Starting cluster update for database: {db_path}")
    logging.info(f"Using clustering data from: {tsv_path}")
    logging.info(f"Duplicate handling strategy: {handle_duplicates}")
    
    # Create a backup if requested
    if backup_first:
        backup_path = f"{db_path}.backup"
        logging.info(f"Creating backup at: {backup_path}")
        shutil.copy2(db_path, backup_path)
    
    con = duckdb.connect(str(db_path))
    
    try:
        logging.info("Beginning database transaction")
        con.execute("BEGIN TRANSACTION")
        
        # Step 1: Load the TSV file into a temp table
        logging.info("Loading cluster data from TSV...")
        con.execute(f"""
            CREATE TEMP TABLE new_clustering AS
            SELECT 
                column0 AS representative_seqhash_id,
                column1 AS seqhash_id
            FROM read_csv_auto('{tsv_path}', sep='\t', header=FALSE)
        """)
        
        total_entries = con.execute("SELECT COUNT(*) FROM new_clustering").fetchone()[0]
        logging.info(f"Loaded {total_entries} entries from clustering TSV")
        
        # Step 2: Check if we should handle duplicates preemptively
        if handle_duplicates == "replace":
            # Drop existing cluster tables
            logging.info("Strategy is 'replace': Dropping existing cluster tables...")
            try:
                con.execute("DROP TABLE IF EXISTS cluster_members")
                con.execute("DROP TABLE IF EXISTS clusters")
            except Exception as e:
                logging.warning(f"Error dropping tables: {str(e)}")
                
            # Recreate the cluster tables
            logging.info("Recreating cluster tables...")
            con.execute("""
                CREATE TABLE clusters (
                    cluster_id VARCHAR PRIMARY KEY,
                    representative_seqhash_id VARCHAR NOT NULL,
                    size INTEGER NOT NULL
                )
            """)
            
            con.execute("""
                CREATE TABLE cluster_members (
                    seqhash_id VARCHAR NOT NULL,
                    cluster_id VARCHAR NOT NULL,
                    PRIMARY KEY (seqhash_id, cluster_id),
                    FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id),
                    FOREIGN KEY (cluster_id) REFERENCES clusters(cluster_id)
                )
            """)
            
        else:
            # Check for potential duplicates if we're not replacing
            cluster_members_count = 0
            try:
                cluster_members_count = con.execute("SELECT COUNT(*) FROM cluster_members").fetchone()[0]
            except:
                # Table might not exist
                logging.info("cluster_members table doesn't exist yet, will create it")
                con.execute("""
                    CREATE TABLE IF NOT EXISTS clusters (
                        cluster_id VARCHAR PRIMARY KEY,
                        representative_seqhash_id VARCHAR NOT NULL,
                        size INTEGER NOT NULL
                    )
                """)
                
                con.execute("""
                    CREATE TABLE IF NOT EXISTS cluster_members (
                        seqhash_id VARCHAR PRIMARY KEY,
                        cluster_id VARCHAR NOT NULL,
                        FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id),
                        FOREIGN KEY (cluster_id) REFERENCES clusters(cluster_id)
                    )
                """)
            
            if cluster_members_count > 0 and handle_duplicates == "error":
                # Look for potential duplicates that would cause errors
                potential_duplicates = con.execute("""
                    SELECT COUNT(*) FROM new_clustering nc
                    JOIN cluster_members cm ON nc.seqhash_id = cm.seqhash_id
                """).fetchone()[0]
                
                if potential_duplicates > 0:
                    logging.warning(f"Found {potential_duplicates} potential duplicate entries that could cause errors")
                    # We'll let the error happen naturally during insert
        
        # Step 3: Identify missing sequences
        logging.info("Identifying missing sequences...")
        
        # Find missing sequences
        con.execute("""
            CREATE TEMP TABLE missing_sequences AS
            SELECT DISTINCT seqhash_id
            FROM new_clustering
            WHERE NOT EXISTS (
                SELECT 1 FROM sequences s WHERE s.seqhash_id = new_clustering.seqhash_id
            )
        """)
        
        # Find missing representatives
        con.execute("""
            CREATE TEMP TABLE missing_representatives AS
            SELECT DISTINCT representative_seqhash_id
            FROM new_clustering
            WHERE NOT EXISTS (
                SELECT 1 FROM sequences s WHERE s.seqhash_id = new_clustering.representative_seqhash_id
            )
        """)
        
        missing_seq_count = con.execute("SELECT COUNT(*) FROM missing_sequences").fetchone()[0]
        missing_rep_count = con.execute("SELECT COUNT(*) FROM missing_representatives").fetchone()[0]
        
        logging.info(f"Found {missing_seq_count} unique sequences missing from database")
        logging.info(f"Found {missing_rep_count} unique representatives missing from database")
        
        # Log some examples of missing sequences
        if missing_seq_count > 0:
            missing_examples = con.execute(
                "SELECT seqhash_id FROM missing_sequences LIMIT 10"
            ).fetchall()
            logging.info("Examples of missing sequences:")
            for i, (seq_id,) in enumerate(missing_examples, 1):
                logging.info(f"  {i}. {seq_id}")
        
        # Step 4: Filter to only valid entries
        logging.info("Filtering to valid entries only...")
        con.execute("""
            CREATE TEMP TABLE valid_clustering AS
            SELECT nc.*
            FROM new_clustering nc
            WHERE EXISTS (SELECT 1 FROM sequences s1 WHERE s1.seqhash_id = nc.seqhash_id)
              AND EXISTS (SELECT 1 FROM sequences s2 WHERE s2.seqhash_id = nc.representative_seqhash_id)
        """)
        
        valid_count = con.execute("SELECT COUNT(*) FROM valid_clustering").fetchone()[0]
        logging.info(f"Retained {valid_count}/{total_entries} entries after filtering")
        
        # Step 5: Update sequences with representative info
        logging.info("Updating sequences with representative information...")
        
        # Check if is_representative column exists
        column_exists = con.execute("""
            SELECT COUNT(*) FROM pragma_table_info('sequences') 
            WHERE name = 'is_representative'
        """).fetchone()[0] > 0
        
        if column_exists:
            con.execute("""
                UPDATE sequences
                SET repseq_id = vc.representative_seqhash_id,
                    is_representative = (sequences.seqhash_id = vc.representative_seqhash_id)
                FROM valid_clustering vc
                WHERE sequences.seqhash_id = vc.seqhash_id
            """)
        else:
            # Update without is_representative column
            con.execute("""
                UPDATE sequences
                SET repseq_id = vc.representative_seqhash_id
                FROM valid_clustering vc
                WHERE sequences.seqhash_id = vc.seqhash_id
            """)
            logging.info("Note: is_representative column not found in sequences table, skipping that update")
        
        updated_count = con.execute(
            "SELECT COUNT(*) FROM sequences WHERE repseq_id IS NOT NULL"
        ).fetchone()[0]
        logging.info(f"Updated {updated_count} sequences with representative info")
        
        # Step 6: Insert clusters
        logging.info("Creating clusters...")
        if handle_duplicates == "ignore":
            # Use INSERT OR IGNORE syntax
            con.execute("""
                INSERT OR IGNORE INTO clusters (cluster_id, representative_seqhash_id, size)
                SELECT 
                    representative_seqhash_id AS cluster_id,
                    representative_seqhash_id,
                    COUNT(*) AS size
                FROM valid_clustering
                GROUP BY representative_seqhash_id
            """)
        else:
            # Normal insert (will error on duplicates if they exist and handle_duplicates="error")
            con.execute("""
                INSERT INTO clusters (cluster_id, representative_seqhash_id, size)
                SELECT 
                    representative_seqhash_id AS cluster_id,
                    representative_seqhash_id,
                    COUNT(*) AS size
                FROM valid_clustering
                GROUP BY representative_seqhash_id
            """)
        
        cluster_count = con.execute("SELECT COUNT(*) FROM clusters").fetchone()[0]
        logging.info(f"Created {cluster_count} clusters")
        
        # Step 7: Insert cluster memberships with explicit duplicate handling
        logging.info("Adding cluster memberships with careful duplicate handling...")
        
        # Create a temporary table of what we want to insert
        con.execute("""
            CREATE TEMP TABLE mem_to_insert AS
            SELECT DISTINCT 
                seqhash_id, 
                representative_seqhash_id AS cluster_id
            FROM valid_clustering
        """)
        
        # For ignore or error strategies, we need to check for existing entries
        if handle_duplicates != "replace" and con.execute("SELECT COUNT(*) FROM cluster_members").fetchone()[0] > 0:
            # Remove entries that would cause duplicates
            logging.info("Checking for duplicate memberships...")
            con.execute("""
                CREATE TEMP TABLE existing_members AS
                SELECT seqhash_id, cluster_id FROM cluster_members
            """)
            
            # Count how many would cause duplicates
            duplicates = con.execute("""
                SELECT COUNT(*) FROM mem_to_insert m
                WHERE EXISTS (
                    SELECT 1 FROM existing_members e
                    WHERE e.seqhash_id = m.seqhash_id AND e.cluster_id = m.cluster_id
                )
            """).fetchone()[0]
            logging.info(f"Found {duplicates} entries that would cause duplicates")
            
            if handle_duplicates == "ignore":
                # Remove entries that would cause duplicates
                con.execute("""
                    DELETE FROM mem_to_insert m
                    WHERE EXISTS (
                        SELECT 1 FROM existing_members e
                        WHERE e.seqhash_id = m.seqhash_id AND e.cluster_id = m.cluster_id
                    )
                """)
                logging.info(f"Filtered out {duplicates} duplicate entries")
        
        # Now insert only the non-duplicate entries
        insert_count = 0
        if handle_duplicates == "ignore" or handle_duplicates == "replace":
            # Safe insert with no errors
            con.execute("""
                INSERT OR IGNORE INTO cluster_members (seqhash_id, cluster_id)
                SELECT seqhash_id, cluster_id 
                FROM mem_to_insert
            """)
            insert_count = con.execute("SELECT COUNT(*) FROM mem_to_insert").fetchone()[0]
        else:
            # With "error" strategy, we'll just let it fail if there are duplicates
            con.execute("""
                INSERT INTO cluster_members (seqhash_id, cluster_id)
                SELECT seqhash_id, cluster_id
                FROM mem_to_insert
            """)
            insert_count = con.execute("SELECT COUNT(*) FROM mem_to_insert").fetchone()[0]
        
        logging.info(f"Inserted {insert_count} cluster memberships")
        
        # Log final counts
        member_count = con.execute("SELECT COUNT(*) FROM cluster_members").fetchone()[0]
        logging.info(f"Final cluster_members count: {member_count}")
        
        # Log cluster statistics
        avg_size = con.execute("SELECT AVG(size) FROM clusters").fetchone()[0]
        max_size = con.execute("SELECT MAX(size) FROM clusters").fetchone()[0]
        logging.info(f"Average cluster size: {avg_size:.2f}")
        logging.info(f"Largest cluster size: {max_size}")
        
        # Commit the transaction
        con.execute("COMMIT")
        logging.info("Successfully updated cluster information")
        logging.info(f"Log file saved to: {log_path}")
        
    except Exception as e:
        con.execute("ROLLBACK")
        logging.error(f"Error updating cluster information: {str(e)}")
        raise
    finally:
        # Clean up temp tables
        try:
            for table in ["new_clustering", "missing_sequences", "missing_representatives", 
                         "valid_clustering", "mem_to_insert", "existing_members"]:
                con.execute(f"DROP TABLE IF EXISTS {table}")
        except:
            pass
        con.close()
    if backup_first:
        return backup_path
    return db_path

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
        run_bases BIGINT,
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

def verify_database(db_path: Union[str, Path]):
    """Verify key tables and their row counts"""
    db_path = str(db_path)
    
    con = duckdb.connect(db_path)
    try:
        # Get all tables
        tables = con.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        ).fetchall()
        
        print("Database contains the following tables:")
        for table_row in tables:
            table_name = table_row[0]
            row_count = con.execute(f"SELECT COUNT(*) FROM {table_name}").fetchone()[0]
            print(f"  - {table_name}: {row_count} rows")
    finally:
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

def validate_duckdb_schema(db_path):
    """Validate schema and relationships in the database."""
    print(f"Validating database: {db_path}")
    
    with duckdb.connect(db_path) as con:
        # 1. Get list of all tables
        print("\n=== Tables in database ===")
        tables = con.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
        for table in tables:
            table_name = table[0]
            count = con.execute(f"SELECT COUNT(*) FROM {table_name}").fetchone()[0]
            print(f"Table: {table_name} - {count} rows")
        
        # 2. Check schema of each table
        print("\n=== Schema for each table ===")
        for table in tables:
            table_name = table[0]
            print(f"\nSchema for {table_name}:")
            schema = con.execute(f"PRAGMA table_info({table_name})").fetchall()
            for col in schema:
                print(f"  {col[1]} ({col[2]}){' PRIMARY KEY' if col[5] > 0 else ''}")
        
        # 3. Check the gene-protein relationships
        print("\n=== Gene-Protein Relationships ===")
        try:
            gene_protein_stats = con.execute("""
                SELECT 
                    COUNT(DISTINCT gene_seqhash_id) AS total_genes,
                    COUNT(DISTINCT protein_seqhash_id) AS total_proteins,
                    COUNT(*) AS total_relationships
                FROM gene_protein_map
            """).fetchone()
            
            print(f"Total genes: {gene_protein_stats[0]}")
            print(f"Total proteins: {gene_protein_stats[1]}")
            print(f"Total gene-protein relationships: {gene_protein_stats[2]}")
            
            # Check for genes with multiple proteins
            multi_protein_genes = con.execute("""
                SELECT gene_seqhash_id, COUNT(protein_seqhash_id) as protein_count
                FROM gene_protein_map
                GROUP BY gene_seqhash_id
                HAVING COUNT(protein_seqhash_id) > 1
                ORDER BY COUNT(protein_seqhash_id) DESC
                LIMIT 5
            """).fetchall()
            
            if multi_protein_genes:
                print("\nTop 5 genes with multiple proteins:")
                for gene in multi_protein_genes:
                    print(f"  Gene {gene[0]} has {gene[1]} proteins")
                    # Sample of proteins for this gene
                    proteins = con.execute(f"""
                        SELECT protein_seqhash_id 
                        FROM gene_protein_map 
                        WHERE gene_seqhash_id = '{gene[0]}'
                        LIMIT 3
                    """).fetchall()
                    for protein in proteins:
                        print(f"    - {protein[0]}")
        except Exception as e:
            print(f"Error checking gene-protein relationships: {e}")
        
        # 4. Check expression data and linkage to genes
        print("\n=== Expression Data ===")
        try:
            expr_stats = con.execute("""
                SELECT 
                    COUNT(DISTINCT gene_seqhash_id) AS genes_with_expression,
                    COUNT(DISTINCT sample_id) AS samples_with_expression,
                    AVG(tpm) AS avg_tpm
                FROM expression
            """).fetchone()
            
            print(f"Genes with expression data: {expr_stats[0]}")
            print(f"Samples with expression data: {expr_stats[1]}")
            print(f"Average TPM: {expr_stats[2]:.2f}")
            
            # Check gene-expression-protein linkage
            gene_expr_protein = con.execute("""
                SELECT 
                    COUNT(DISTINCT e.gene_seqhash_id) AS genes_with_expr_and_protein,
                    COUNT(DISTINCT gpm.protein_seqhash_id) AS proteins_linked_to_expr
                FROM expression e
                JOIN gene_protein_map gpm ON e.gene_seqhash_id = gpm.gene_seqhash_id
            """).fetchone()
            
            print(f"Genes with both expression and protein mappings: {gene_expr_protein[0]}")
            print(f"Proteins linked to genes with expression: {gene_expr_protein[1]}")
        except Exception as e:
            print(f"Error checking expression data: {e}")
        
        # 5. Check annotations
        print("\n=== Annotation Data ===")
        try:
            anno_stats = con.execute("""
                SELECT 
                    COUNT(*) AS total_annotations,
                    COUNT(DISTINCT sample_id) AS samples_with_annotations
                FROM annotations
            """).fetchone()
            
            print(f"Total annotated sequences: {anno_stats[0]}")
            print(f"Samples with annotations: {anno_stats[1]}")
            
            # Check annotation-sequence-gene-expression linkage
            anno_gene_expr = con.execute("""
                SELECT 
                    COUNT(DISTINCT a.seqhash_id) AS annotated_proteins,
                    COUNT(DISTINCT gpm.gene_seqhash_id) AS genes_of_annotated_proteins,
                    COUNT(DISTINCT e.gene_seqhash_id) AS genes_with_annotation_and_expression
                FROM annotations a
                JOIN sequences s ON a.seqhash_id = s.seqhash_id
                JOIN gene_protein_map gpm ON s.seqhash_id = gpm.protein_seqhash_id
                LEFT JOIN expression e ON gpm.gene_seqhash_id = e.gene_seqhash_id
            """).fetchone()
            
            print(f"Annotated proteins: {anno_gene_expr[0]}")
            print(f"Genes of annotated proteins: {anno_gene_expr[1]}")
            print(f"Genes with both annotation and expression: {anno_gene_expr[2]}")
        except Exception as e:
            print(f"Error checking annotation data: {e}")
        
        # 6. Check clusters
        print("\n=== Cluster Data ===")
        try:
            cluster_stats = con.execute("""
                SELECT 
                    COUNT(DISTINCT cluster_id) AS total_clusters,
                    AVG(size) AS avg_cluster_size,
                    MAX(size) AS largest_cluster
                FROM clusters
            """).fetchone()
            
            print(f"Total clusters: {cluster_stats[0]}")
            print(f"Average cluster size: {cluster_stats[1]:.2f}")
            print(f"Largest cluster size: {cluster_stats[2]}")
            
            # Check cluster members
            member_stats = con.execute("""
                SELECT 
                    COUNT(DISTINCT seqhash_id) AS unique_members,
                    COUNT(*) AS total_membership_records
                FROM cluster_members
            """).fetchone()
            
            print(f"Unique sequences in clusters: {member_stats[0]}")
            print(f"Total cluster membership records: {member_stats[1]}")
            
            # Get top clusters
            top_clusters = con.execute("""
                SELECT cluster_id, size
                FROM clusters
                ORDER BY size DESC
                LIMIT 3
            """).fetchall()
            
            if top_clusters:
                print("\nTop 3 largest clusters:")
                for cluster in top_clusters:
                    print(f"  Cluster {cluster[0]}: {cluster[1]} members")
        except Exception as e:
            print(f"Error checking cluster data: {e}")
