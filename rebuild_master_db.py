#!/usr/bin/env python3
"""
Rebuild master database from samples that match the clustering.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from planter.database.utils.duckdb_utils import merge_duckdbs

def main():
    # Read list of samples that have DuckDB files
    samples_file = Path('/tmp/samples_with_duckdb.txt')

    if not samples_file.exists():
        print("Error: /tmp/samples_with_duckdb.txt not found")
        print("Run this first:")
        print("  find /mnt/data4/recombia.planter -name '*.pep' -path '*/transdecoder/*.pep' | grep -v transdecoder_dir | cut -d'/' -f5 | sort > /tmp/completed_samples.txt")
        print("  while read sample; do if [ -f \"/mnt/data4/recombia.planter/$sample/$sample.duckdb\" ]; then echo \"$sample\"; fi; done < /tmp/completed_samples.txt > /tmp/samples_with_duckdb.txt")
        sys.exit(1)

    with open(samples_file) as f:
        samples = [line.strip() for line in f if line.strip()]

    print(f"Found {len(samples)} samples with DuckDB files")

    # Build paths
    base_dir = Path('/mnt/data4/recombia.planter')
    duckdb_paths = []

    for sample in samples:
        db_path = base_dir / sample / f"{sample}.duckdb"
        if db_path.exists():
            duckdb_paths.append(str(db_path))
        else:
            print(f"Warning: {db_path} not found, skipping")

    print(f"\nMerging {len(duckdb_paths)} sample databases into master_rebuilt.duckdb...")

    # Merge into new database
    output_db = Path('/mnt/data4/master_rebuilt.duckdb')

    if output_db.exists():
        print(f"\nWarning: {output_db} already exists. Remove it first or choose different name.")
        response = input("Delete existing file and continue? (yes/no): ")
        if response.lower() != 'yes':
            print("Aborted.")
            sys.exit(1)
        output_db.unlink()

    # Get schema path - use the initial schema and let merge_duckdbs upgrade it
    schema_path = Path('/home/ubuntu/planter/planter/database/schema/migrations/001_initial_schema.sql')

    merge_duckdbs(
        duckdb_paths=duckdb_paths,
        master_db_path=str(output_db),
        schema_sql_path=str(schema_path),
        upgrade_schema=True
    )

    print(f"\n✓ Master database created: {output_db}")
    print(f"\nNext steps:")
    print(f"1. Load cluster data:")
    print(f"   python load_clusters_v6.py")
    print(f"\n2. If successful, replace old database:")
    print(f"   mv /mnt/data4/master.duckdb /mnt/data4/master.duckdb.backup")
    print(f"   mv {output_db} /mnt/data4/master.duckdb")

if __name__ == '__main__':
    main()
