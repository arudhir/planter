#!/usr/bin/env python3
"""
Load cluster data from repseq_v6 into the master database.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from planter.database.utils.duckdb_utils import update_clusters

def main():
    db_path = '/mnt/data4/master.duckdb'
    tsv_path = '/mnt/data4/planter_outputs/repseq_v6/output130/newClusterDB.tsv'

    print(f"Loading cluster data into master database...")
    print(f"Database: {db_path}")
    print(f"Cluster TSV: {tsv_path}")
    print()

    update_clusters(
        db_path=db_path,
        tsv_path=tsv_path,
        backup_first=True,
        handle_duplicates='replace'  # Replace existing cluster data
    )

    print("\n✓ Cluster data loaded successfully!")
    print(f"\nBackup created at: {db_path}.backup")

if __name__ == '__main__':
    main()
