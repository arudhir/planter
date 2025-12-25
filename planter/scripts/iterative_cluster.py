#!/usr/bin/env python3

import argparse
import glob
import os
import subprocess
import sys
from pathlib import Path

from tqdm import tqdm


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Batch update script for MMseqs clustering"
    )
    parser.add_argument(
        "-g",
        "--glob-pattern",
        required=True,
        help='Glob pattern for input files (e.g., "*.pep")',
    )
    parser.add_argument(
        "-o", "--output-dir", required=True, help="Base output directory"
    )
    parser.add_argument(
        "-i",
        "--initial-reps",
        default="",
        help="Initial old representative sequences file",
    )
    parser.add_argument(
        "-d",
        "--database",
        default="",
        help="DuckDB database path (if provided, cluster data will be loaded after clustering)",
    )
    return parser.parse_args()


def setup_output_directory(base_dir):
    """Create the base output directory if it doesn't exist."""
    os.makedirs(base_dir, exist_ok=True)
    return base_dir


def get_sorted_files(glob_pattern):
    """Get and sort the list of files matching the glob pattern."""
    files = glob.glob(glob_pattern)
    if not files:
        print(f"No files found matching the pattern {glob_pattern}")
        sys.exit(1)
    return sorted(files)


def run_mmseqs_update(script_path, old_reps, new_file, output_dir):
    """Run the mmseqs_cluster_update.py script with the given parameters."""
    print("=" * 39)
    print(f"Processing file: {new_file}")
    print(f"Old representative sequences: {old_reps}")
    print(f"Output directory: {output_dir}")
    print("=" * 39)

    # Get the directory containing the current script
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # Construct path to mmseqs_cluster_update.py in the same directory
    mmseqs_script = os.path.join(current_dir, "mmseqs_cluster_update.py")

    try:
        subprocess.run(
            [
                sys.executable,
                mmseqs_script,
                "-i",
                old_reps,
                "-o",
                output_dir,
                new_file,
            ],
            check=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"Error running mmseqs_cluster_update.py: {e}")
        sys.exit(1)


def load_clusters_to_database(database_path, cluster_tsv, output_dir):
    """Load cluster data into DuckDB database."""
    print("\n" + "=" * 80)
    print("Loading cluster data into database")
    print("=" * 80)
    print(f"Database: {database_path}")
    print(f"Cluster TSV: {cluster_tsv}")

    # Import here to avoid requiring database dependencies if not using this feature
    sys.path.insert(0, str(Path(__file__).parent.parent.parent))
    from planter.database.builder import SequenceDBBuilder

    try:
        with SequenceDBBuilder(database_path, output_dir=output_dir) as builder:
            # Check current state
            result = builder.con.execute("""
                SELECT
                    (SELECT COUNT(*) FROM sequences) as total_sequences,
                    (SELECT COUNT(*) FROM cluster_members) as before_clustered
            """).fetchone()

            print(f"\nBefore loading:")
            print(f"  Total sequences: {result[0]:,}")
            print(f"  Sequences with cluster assignments: {result[1]:,}")

            # Clear existing cluster data to avoid conflicts
            print("\nClearing existing cluster data...")
            builder.con.execute("DELETE FROM cluster_members")
            builder.con.execute("DELETE FROM clusters")
            builder.con.execute("UPDATE sequences SET is_representative = FALSE")

            # Load new cluster data
            print("Loading cluster data from TSV...")
            builder.load_clusters_from_tsv(cluster_tsv)

            # Check new state
            result = builder.con.execute("""
                SELECT
                    (SELECT COUNT(*) FROM cluster_members) as after_clustered,
                    (SELECT COUNT(*) FROM clusters) as total_clusters,
                    (SELECT COUNT(*) FROM sequences WHERE is_representative = TRUE) as representatives
            """).fetchone()

            print(f"\nAfter loading:")
            print(f"  Sequences with cluster assignments: {result[0]:,}")
            print(f"  Total clusters: {result[1]:,}")
            print(f"  Representative sequences: {result[2]:,}")
            print("\n✓ Cluster data loaded successfully!")
            print("=" * 80)

    except Exception as e:
        print(f"\n✗ Error loading cluster data: {e}")
        print("Clustering completed but cluster data was not loaded into database.")
        print("You can load it manually later using:")
        print(f"  python scripts/load_clusters.py -d {database_path} -t {cluster_tsv}")
        raise


def main():
    args = parse_arguments()
    base_dir = setup_output_directory(args.output_dir)
    files = get_sorted_files(args.glob_pattern)

    # Initialize old_rep_seqs with the provided initial reps or the first file
    old_rep_seqs = args.initial_reps
    if not old_rep_seqs:
        old_rep_seqs = files[0]
        print(
            f"No initial representative sequences provided. Using {files[0]} as initial --old input."
        )
        files = files[1:]  # Skip the first file if we're using it as initial input

    # Process each file
    for i, file in enumerate(
        tqdm(files, desc="Processing files", unit="file"), start=1
    ):
        output_dir = os.path.join(base_dir, f"output{i}")
        os.makedirs(output_dir, exist_ok=True)

        run_mmseqs_update(
            script_path="mmseqs_cluster_update.py",
            old_reps=old_rep_seqs,
            new_file=file,
            output_dir=output_dir,
        )

        # Update old_rep_seqs for the next iteration
        old_rep_seqs = os.path.join(output_dir, "newRepSeqDB.fasta")

    # After all clustering is complete, load cluster data into database if requested
    if args.database:
        final_output_num = len(files)
        final_output_dir = os.path.join(base_dir, f"output{final_output_num}")
        cluster_tsv = os.path.join(final_output_dir, "newClusterDB.tsv")

        if not os.path.exists(cluster_tsv):
            print(f"\nWarning: Cluster TSV not found at {cluster_tsv}")
            print("Cluster data will not be loaded into database.")
        else:
            load_clusters_to_database(args.database, cluster_tsv, base_dir)


if __name__ == "__main__":
    main()
