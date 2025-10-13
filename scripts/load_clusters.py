#!/usr/bin/env python3
"""
Load cluster information from MMSeqs2 TSV into DuckDB database.

This script imports cluster membership data that was computed during clustering
but never loaded into the database.
"""
import argparse
import sys
from pathlib import Path

# Add planter to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from planter.database.builder import SequenceDBBuilder


def main():
    parser = argparse.ArgumentParser(
        description="Load MMSeqs2 cluster data into DuckDB"
    )
    parser.add_argument(
        "-d", "--database",
        required=True,
        help="Path to DuckDB database"
    )
    parser.add_argument(
        "-t", "--tsv",
        required=True,
        help="Path to cluster TSV file (e.g., output130/newClusterDB.tsv)"
    )
    parser.add_argument(
        "-o", "--output-dir",
        default="/mnt/data4/planter_outputs",
        help="Output directory (default: /mnt/data4/planter_outputs)"
    )

    args = parser.parse_args()

    db_path = Path(args.database)
    tsv_path = Path(args.tsv)

    if not db_path.exists():
        print(f"Error: Database not found at {db_path}")
        sys.exit(1)

    if not tsv_path.exists():
        print(f"Error: TSV file not found at {tsv_path}")
        sys.exit(1)

    print(f"Database: {db_path}")
    print(f"Cluster TSV: {tsv_path}")
    print()

    # Load clusters
    print("Loading cluster data...")
    with SequenceDBBuilder(str(db_path), output_dir=args.output_dir) as builder:
        # First, check current state
        result = builder.con.execute("""
            SELECT
                (SELECT COUNT(*) FROM sequences) as total_sequences,
                (SELECT COUNT(*) FROM cluster_members) as clustered_sequences,
                (SELECT COUNT(*) FROM clusters) as total_clusters
        """).fetchone()

        print(f"Before loading:")
        print(f"  Total sequences: {result[0]:,}")
        print(f"  Sequences with cluster assignments: {result[1]:,}")
        print(f"  Total clusters: {result[2]:,}")
        print()

        # Clear existing cluster data
        print("Clearing existing cluster data...")
        builder.con.execute("DELETE FROM cluster_members")
        builder.con.execute("DELETE FROM clusters")
        builder.con.execute("UPDATE sequences SET is_representative = FALSE")
        print("✓ Cleared")
        print()

        # Load new cluster data
        print("Loading cluster data from TSV...")
        builder.load_clusters_from_tsv(str(tsv_path))
        print("✓ Loaded")
        print()

        # Check new state
        result = builder.con.execute("""
            SELECT
                (SELECT COUNT(*) FROM sequences) as total_sequences,
                (SELECT COUNT(*) FROM cluster_members) as clustered_sequences,
                (SELECT COUNT(*) FROM clusters) as total_clusters,
                (SELECT COUNT(*) FROM sequences WHERE is_representative = TRUE) as representatives
        """).fetchone()

        print(f"After loading:")
        print(f"  Total sequences: {result[0]:,}")
        print(f"  Sequences with cluster assignments: {result[1]:,}")
        print(f"  Total clusters: {result[2]:,}")
        print(f"  Representative sequences: {result[3]:,}")
        print()

        coverage = (result[1] / result[0] * 100) if result[0] > 0 else 0
        print(f"  Coverage: {coverage:.1f}% of sequences have cluster assignments")
        print()
        print("✓ Cluster data loaded successfully!")


if __name__ == "__main__":
    main()
