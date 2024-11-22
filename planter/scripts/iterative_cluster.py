#!/usr/bin/env python3

import argparse
import glob
import os
import subprocess
import sys
from pathlib import Path

def parse_arguments():
    parser = argparse.ArgumentParser(description='Batch update script for MMseqs clustering')
    parser.add_argument('-g', '--glob-pattern', required=True,
                        help='Glob pattern for input files (e.g., "*.pep")')
    parser.add_argument('-o', '--output-dir', required=True,
                        help='Base output directory')
    parser.add_argument('-i', '--initial-reps', default='',
                        help='Initial old representative sequences file')
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
        subprocess.run([
            sys.executable, mmseqs_script,
            "--old", old_reps,
            "--new", new_file,
            "--output", output_dir
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running mmseqs_cluster_update.py: {e}")
        sys.exit(1)

def main():
    args = parse_arguments()
    base_dir = setup_output_directory(args.output_dir)
    files = get_sorted_files(args.glob_pattern)

    # Initialize old_rep_seqs with the provided initial reps or the first file
    old_rep_seqs = args.initial_reps
    if not old_rep_seqs:
        old_rep_seqs = files[0]
        print(f"No initial representative sequences provided. Using {files[0]} as initial --old input.")
        files = files[1:]  # Skip the first file if we're using it as initial input

    # Process each file
    for i, file in enumerate(files, start=1):
        output_dir = os.path.join(base_dir, f"output{i}")
        os.makedirs(output_dir, exist_ok=True)

        run_mmseqs_update(
            script_path="mmseqs_cluster_update.py",
            old_reps=old_rep_seqs,
            new_file=file,
            output_dir=output_dir
        )

        # Update old_rep_seqs for the next iteration
        old_rep_seqs = os.path.join(output_dir, "newRepSeqDB.fasta")

if __name__ == "__main__":
    main()