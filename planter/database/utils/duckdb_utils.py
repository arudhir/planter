#!/usr/bin/env python3
"""
DuckDB utility functions for merging and updating databases.
"""
import duckdb
import pandas as pd
from pathlib import Path
from typing import List, Union


def merge_duckdbs(
    duckdb_paths: List[Union[str, Path]],
    master_db_path: Union[str, Path],
    schema_sql_path: Union[str, Path]
) -> str:
    """
    Merge multiple DuckDB databases into a master DuckDB.
    
    Args:
        duckdb_paths: List of paths to source DuckDB files.
        master_db_path: Path to the master (merged) DuckDB.
        schema_sql_path: Path to the SQL file defining the schema.
    
    Returns:
        Path to the master (merged) DuckDB.
    
    The function:
      - Creates (or opens) the master database.
      - Executes the schema SQL to create tables if they don't exist.
      - Iterates through each source database, attaches it,
        and inserts data into the master tables in dependency order.
      - Uses INSERT OR IGNORE to avoid duplicate key errors.
      - Detaches each source database after merging.
    """
    master_db_path = str(master_db_path)
    schema_sql_path = Path(schema_sql_path)
    
    # Read the schema SQL
    schema_sql = schema_sql_path.read_text()
    
    with duckdb.connect(master_db_path) as master_conn:
        # Set up the schema in the master database
        master_conn.execute(schema_sql)
        
        # Process each source DuckDB
        for i, source_db in enumerate(duckdb_paths):
            alias = f"db{i}"
            source_db_str = str(source_db)
            print(f"Attaching {source_db_str} as {alias}...")
            master_conn.execute(f"ATTACH '{source_db_str}' AS {alias};")
            
            # Insert data in dependency order
            master_conn.execute(f"""
                INSERT OR IGNORE INTO sra_metadata
                SELECT * FROM {alias}.sra_metadata;
            """)
            master_conn.execute(f"""
                INSERT OR IGNORE INTO sequences
                SELECT * FROM {alias}.sequences;
            """)
            master_conn.execute(f"""
                INSERT OR IGNORE INTO annotations
                SELECT * FROM {alias}.annotations;
            """)
            master_conn.execute(f"""
                INSERT OR IGNORE INTO go_terms
                SELECT * FROM {alias}.go_terms;
            """)
            master_conn.execute(f"""
                INSERT OR IGNORE INTO ec_numbers
                SELECT * FROM {alias}.ec_numbers;
            """)
            master_conn.execute(f"""
                INSERT OR IGNORE INTO kegg_info
                SELECT * FROM {alias}.kegg_info;
            """)
            master_conn.execute(f"""
                INSERT OR IGNORE INTO clusters
                SELECT * FROM {alias}.clusters;
            """)
            master_conn.execute(f"""
                INSERT OR IGNORE INTO cluster_members
                SELECT * FROM {alias}.cluster_members;
            """)
            
            # Check if expression table exists in both databases
            expression_exists = master_conn.execute("""
                SELECT COUNT(*) FROM sqlite_master 
                WHERE type='table' AND name='expression'
            """).fetchone()[0]
            
            expression_exists_in_source = master_conn.execute(f"""
                SELECT COUNT(*) FROM {alias}.sqlite_master 
                WHERE type='table' AND name='expression'
            """).fetchone()[0]
            
            if expression_exists and expression_exists_in_source:
                master_conn.execute(f"""
                    INSERT OR IGNORE INTO expression
                    SELECT * FROM {alias}.expression;
                """)
            
            # Check if gene_protein_map table exists in both databases
            gene_protein_map_exists = master_conn.execute("""
                SELECT COUNT(*) FROM sqlite_master 
                WHERE type='table' AND name='gene_protein_map'
            """).fetchone()[0]
            
            gene_protein_map_exists_in_source = master_conn.execute(f"""
                SELECT COUNT(*) FROM {alias}.sqlite_master 
                WHERE type='table' AND name='gene_protein_map'
            """).fetchone()[0]
            
            if gene_protein_map_exists and gene_protein_map_exists_in_source:
                master_conn.execute(f"""
                    INSERT OR IGNORE INTO gene_protein_map
                    SELECT * FROM {alias}.gene_protein_map;
                """)
            
            master_conn.execute(f"DETACH {alias};")
            print(f"Finished merging {source_db_str}\n")
        
        # Optional commit; DuckDB auto-commits by default.
        master_conn.commit()
    
    print("All databases have been merged into:", master_db_path)
    return master_db_path


def update_duckdb_with_cluster_info(db_path: Union[str, Path], tsv_path: Union[str, Path]) -> None:
    """
    Updates a DuckDB database with clustering information from an MMSeqs2 TSV file.

    Args:
        db_path: Path to the DuckDB database file.
        tsv_path: Path to the MMSeqs2 TSV file containing clustering data.
    """
    db_path = str(db_path)
    tsv_path = str(tsv_path)
    
    # Load the MMSeqs2 TSV into a DataFrame
    df = pd.read_csv(tsv_path, sep="\t", names=["representative_seqhash_id", "seqhash_id"])
    df = df.drop_duplicates()
    
    # Connect to DuckDB
    con = duckdb.connect(db_path)

    try:
        print(f"Updating database '{db_path}' with MMSeqs2 clustering information from '{tsv_path}'")
        print(f"Found {len(df)} cluster assignments")
        
        # Begin transaction
        con.execute("BEGIN TRANSACTION")
        
        # Create a temporary table for the MMSeqs2 data
        con.execute("""
            CREATE TEMPORARY TABLE mmseqs2_clusters (
                representative_seqhash_id VARCHAR,
                seqhash_id VARCHAR
            )
        """)

        # Insert data into the temporary table
        con.executemany("INSERT INTO mmseqs2_clusters VALUES (?, ?)", df.values.tolist())

        # First, verify that all seqhash_ids exist in the sequences table
        missing_seqs = con.execute("""
            SELECT DISTINCT mc.seqhash_id 
            FROM mmseqs2_clusters mc
            LEFT JOIN sequences s ON mc.seqhash_id = s.seqhash_id
            WHERE s.seqhash_id IS NULL
        """).fetchall()
        
        if missing_seqs:
            missing_count = len(missing_seqs)
            print(f"Warning: {missing_count} sequence IDs from cluster file not found in the database")
            if missing_count > 10:
                print(f"First 10 missing IDs: {[m[0] for m in missing_seqs[:10]]}")
            else:
                print(f"Missing IDs: {[m[0] for m in missing_seqs]}")
        
        # Check if is_representative column exists in the sequences table
        has_representative_column = con.execute("""
            SELECT COUNT(*) FROM pragma_table_info('sequences') 
            WHERE name = 'is_representative'
        """).fetchone()[0] > 0
        
        if has_representative_column:
            # Mark representative sequences in the sequences table if the column exists
            con.execute("""
                UPDATE sequences 
                SET is_representative = FALSE
            """)
            
            con.execute("""
                UPDATE sequences 
                SET is_representative = TRUE
                WHERE seqhash_id IN (
                    SELECT DISTINCT representative_seqhash_id FROM mmseqs2_clusters
                )
            """)

        # Update the sequences table to assign the representative sequence
        con.execute("""
            UPDATE sequences 
            SET repseq_id = mm.representative_seqhash_id
            FROM mmseqs2_clusters mm
            WHERE sequences.seqhash_id = mm.seqhash_id
        """)

        # Get the count of updated sequences for reporting
        updated_count = con.execute("""
            SELECT COUNT(DISTINCT seqhash_id) FROM mmseqs2_clusters
        """).fetchone()[0]
        
        # Clear existing cluster information to rebuild it
        con.execute("DELETE FROM cluster_members")
        con.execute("DELETE FROM clusters")
        
        # Insert new clusters into the clusters table
        # Use proper cluster_id format: 'CLUSTER_nnn'
        con.execute("""
            INSERT INTO clusters (cluster_id, representative_seqhash_id, size)
            SELECT 
                'CLUSTER_' || ROW_NUMBER() OVER (ORDER BY representative_seqhash_id) as cluster_id,
                representative_seqhash_id, 
                COUNT(*) as size
            FROM mmseqs2_clusters
            GROUP BY representative_seqhash_id
        """)

        # Insert cluster memberships into the cluster_members table
        # Use ON CONFLICT to handle potential duplicates
        con.execute("""
            INSERT INTO cluster_members (seqhash_id, cluster_id)
            SELECT m.seqhash_id, c.cluster_id 
            FROM mmseqs2_clusters m
            JOIN clusters c ON m.representative_seqhash_id = c.representative_seqhash_id
            ON CONFLICT (seqhash_id) DO UPDATE SET cluster_id = EXCLUDED.cluster_id
        """)
        
        # Get cluster statistics for reporting
        cluster_stats = con.execute("""
            SELECT COUNT(*) as total_clusters, AVG(size) as avg_size, MAX(size) as max_size 
            FROM clusters
        """).fetchone()
        
        # Commit transaction
        con.execute("COMMIT")
        
        # Report success statistics
        print(f"Successfully updated {updated_count} sequences with new cluster assignments")
        print(f"Created {cluster_stats[0]} clusters with average size {cluster_stats[1]:.1f} (max: {cluster_stats[2]})")
        
    except Exception as e:
        # Rollback on error
        con.execute("ROLLBACK")
        print(f"Error updating database: {str(e)}")
        raise
    finally:
        # Drop the temporary table and close connection
        con.execute("DROP TABLE IF EXISTS mmseqs2_clusters")
        con.close()

    print(f"Database '{db_path}' successfully updated with MMSeqs2 clustering information.")