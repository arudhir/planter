#!/usr/bin/env python3
"""
Full re-clustering script for MMSeqs2.

Unlike iterative_cluster.py which uses clusterupdate (preserving existing representatives),
this script performs a complete re-clustering from scratch. This ensures the longest
sequences are selected as representatives.

Use this periodically or when you need to ensure optimal representative selection.
"""
import argparse
import glob
import logging
import os
import subprocess
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Full MMSeqs2 re-clustering (not incremental update)"
    )
    parser.add_argument(
        "-g",
        "--glob-pattern",
        required=True,
        help='Glob pattern for input files (e.g., "*.pep")',
    )
    parser.add_argument(
        "-o", "--output-dir", required=True, help="Output directory"
    )
    parser.add_argument(
        "--threads", default=16, type=int, help="Number of threads (default: 16)"
    )
    parser.add_argument(
        "--min-seq-id", default=0.3, type=float, help="Minimum sequence identity (default: 0.3)"
    )
    parser.add_argument(
        "-c", "--coverage", default=0.8, type=float, help="Coverage threshold (default: 0.8)"
    )
    return parser.parse_args()


def combine_fasta_files(input_files, output_file):
    """Combine all input FASTA files into one."""
    logger.info(f"Combining {len(input_files)} files into {output_file}")

    with open(output_file, 'w') as outf:
        for i, input_file in enumerate(input_files, 1):
            if i % 10 == 0:
                logger.info(f"  Processing file {i}/{len(input_files)}")
            with open(input_file, 'r') as inf:
                outf.write(inf.read())

    logger.info(f"Combined FASTA written to {output_file}")


def run_full_clustering(combined_fasta, output_dir, threads, min_seq_id, coverage):
    """Run complete MMSeqs2 clustering from scratch."""
    os.makedirs(output_dir, exist_ok=True)

    seqdb = os.path.join(output_dir, "sequenceDB")
    clusterdb = os.path.join(output_dir, "clusterDB")
    repseqdb = os.path.join(output_dir, "repSeqDB")
    tmp_dir = os.path.join(output_dir, "tmp")

    os.makedirs(tmp_dir, exist_ok=True)

    # Step 1: Create sequence database
    logger.info("Step 1: Creating sequence database")
    subprocess.run(
        ["mmseqs", "createdb", combined_fasta, seqdb],
        check=True
    )

    # Step 2: Cluster with correct parameters
    logger.info("Step 2: Clustering sequences")
    logger.info(f"  Parameters: cluster-mode=2, cov-mode=1, min-seq-id={min_seq_id}, coverage={coverage}")
    subprocess.run([
        "mmseqs", "cluster",
        seqdb, clusterdb, tmp_dir,
        "--cluster-mode", "2",  # Greedy clustering by sequence length
        "--cov-mode", "1",      # Coverage of target (allows shorter to join longer)
        "--min-seq-id", str(min_seq_id),
        "-c", str(coverage),
        "--threads", str(threads)
    ], check=True)

    # Step 3: Create TSV output
    logger.info("Step 3: Creating cluster TSV")
    cluster_tsv = os.path.join(output_dir, "clusters.tsv")
    subprocess.run([
        "mmseqs", "createtsv",
        seqdb, seqdb, clusterdb, cluster_tsv
    ], check=True)

    # Step 4: Extract representative sequences
    logger.info("Step 4: Extracting representative sequences")
    subprocess.run([
        "mmseqs", "result2repseq",
        seqdb, clusterdb, repseqdb
    ], check=True)

    # Step 5: Convert to FASTA
    logger.info("Step 5: Converting representatives to FASTA")
    repseq_fasta = os.path.join(output_dir, "repseq.fasta")
    subprocess.run([
        "mmseqs", "convert2fasta",
        repseqdb, repseq_fasta
    ], check=True)

    # Count clusters
    cluster_count = subprocess.run(
        f"cut -f1 {cluster_tsv} | sort | uniq | wc -l",
        shell=True,
        capture_output=True,
        text=True
    ).stdout.strip()

    logger.info(f"Clustering complete!")
    logger.info(f"  Total clusters: {cluster_count}")
    logger.info(f"  Cluster TSV: {cluster_tsv}")
    logger.info(f"  Representative sequences: {repseq_fasta}")

    return cluster_tsv, repseq_fasta


def main():
    args = parse_arguments()

    # Get input files
    input_files = sorted(glob.glob(args.glob_pattern))
    if not input_files:
        logger.error(f"No files found matching pattern: {args.glob_pattern}")
        sys.exit(1)

    logger.info(f"Found {len(input_files)} input files")

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Combine all FASTA files
    combined_fasta = os.path.join(args.output_dir, "all_sequences.fasta")
    combine_fasta_files(input_files, combined_fasta)

    # Run full clustering
    cluster_tsv, repseq_fasta = run_full_clustering(
        combined_fasta,
        args.output_dir,
        args.threads,
        args.min_seq_id,
        args.coverage
    )

    logger.info("\n" + "="*80)
    logger.info("FULL RE-CLUSTERING COMPLETE")
    logger.info("="*80)
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Cluster assignments: {cluster_tsv}")
    logger.info(f"Representative sequences: {repseq_fasta}")
    logger.info("="*80)


if __name__ == "__main__":
    main()
