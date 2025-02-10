#!/usr/bin/env python3
import argparse
import os
import sys
import duckdb
import shutil
import tempfile

def merge_databases(add_db: str, output_db: str):
    """
    Attach the add_db (new sample) as "new_sample" to the output_db,
    check if its sample_id already exists in the master, and if not,
    merge its tables into the output_db.
    """
    con = duckdb.connect(output_db)
    # Attach the new sample database as "new_sample"
    con.execute(f"ATTACH '{add_db}' AS new_sample;")
    
    # Check if the sample already exists:
    sample_ids = con.execute("SELECT DISTINCT sample_id FROM new_sample.sra_metadata").fetchall()
    for sid_tuple in sample_ids:
        sample_id = sid_tuple[0]
        existing_count = con.execute(
            "SELECT COUNT(*) FROM sra_metadata WHERE sample_id = ?", [sample_id]
        ).fetchone()[0]
        if existing_count > 0:
            print(f"Sample '{sample_id}' already exists in the master DB. Merge aborted.")
            con.execute("DETACH new_sample;")
            con.close()
            sys.exit(0)
    
    try:
        con.execute("BEGIN")
        
        print("Merging sra_metadata...")
        con.execute("""
            INSERT INTO sra_metadata
            SELECT * FROM new_sample.sra_metadata;
        """)
        
        print("Merging sequences...")
        con.execute("""
            INSERT INTO sequences
            SELECT * FROM new_sample.sequences;
        """)
        
        print("Merging annotations...")
        con.execute("""
            INSERT INTO annotations
            SELECT * FROM new_sample.annotations;
        """)
        
        print("Merging go_terms...")
        con.execute("""
            INSERT INTO go_terms
            SELECT * FROM new_sample.go_terms;
        """)
        
        print("Merging ec_numbers...")
        con.execute("""
            INSERT INTO ec_numbers
            SELECT * FROM new_sample.ec_numbers;
        """)
        
        print("Merging clusters...")
        con.execute("""
            INSERT INTO clusters
            SELECT * FROM new_sample.clusters;
        """)
        
        print("Merging cluster_members...")
        con.execute("""
            INSERT INTO cluster_members
            SELECT * FROM new_sample.cluster_members;
        """)
        
        con.execute("COMMIT")
        print("Merge completed successfully.")
    except Exception as e:
        con.execute("ROLLBACK")
        print("Error during merge:", e)
    finally:
        con.execute("DETACH new_sample;")
        con.close()

def main():
    parser = argparse.ArgumentParser(
        description="Merge a current (master) DuckDB with an additional sample DuckDB."
    )
    parser.add_argument(
        "--current",
        required=True,
        help="Path to the current (master) DuckDB file."
    )
    parser.add_argument(
        "--add",
        required=True,
        help="Path to the DuckDB file to add/merge into the master."
    )
    parser.add_argument(
        "--output",
        default=None,
        help=("Path for the output merged DuckDB file. "
              "If not provided, a temporary file will be created and its path printed.")
    )
    
    args = parser.parse_args()
    
    # Check that the input files exist.
    if not os.path.exists(args.current):
        parser.error(f"Current DB file '{args.current}' does not exist.")
    if not os.path.exists(args.add):
        parser.error(f"Add DB file '{args.add}' does not exist.")
    
    # Determine the output file path.
    # If --output is not provided, create a temporary file (which won't be auto-deleted).
    output_path = args.output
    if output_path is None:
        tmp = tempfile.NamedTemporaryFile(prefix="merged_", suffix=".duckdb", delete=False)
        output_path = tmp.name
        tmp.close()
        print(f"No output file provided. Using temporary file: {output_path}")
    else:
        print(f"Using output file: {output_path}")
    
    # Copy the current (master) DB to the output location so we do not modify the original.
    try:
        shutil.copyfile(args.current, output_path)
        print(f"Copied current DB '{args.current}' to '{output_path}'.")
    except Exception as e:
        print(f"Error copying current DB: {e}")
        return
    
    # Merge the add database into the output database.
    merge_databases(args.add, output_path)
    
    print(f"Merged database is available at: {output_path}")

if __name__ == "__main__":
    main()

