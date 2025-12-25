#!/usr/bin/env python3
"""
Standalone script to re-cluster existing sequences with --cluster-mode 2
to prefer longer sequences as representatives.

This script is useful for fixing clustering issues where truncated sequences
were selected as representatives instead of full-length sequences.

Usage:
    python recluster_repseq.py --input repseq.faa --output output_dir --db master.duckdb

On EC2:
    python /usr/src/planter/planter/scripts/recluster_repseq.py \
        --input /mnt/data4/repseq.faa \
        --output /mnt/data4/reclustered \
        --db /mnt/data4/master.duckdb
"""

import argparse
import os
import sys
import subprocess
import logging
import shutil
from pathlib import Path

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def run_command(cmd, description):
    """Run a command and log its execution."""
    logger.info(f"{description}")
    logger.info(f"Running: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error(f"Command failed with return code {result.returncode}")
        logger.error(f"STDERR: {result.stderr}")
        sys.exit(1)

    if result.stdout:
        logger.info(f"STDOUT: {result.stdout}")

    return result


def recluster_sequences(input_fasta, output_dir, db_path=None):
    """
    Re-cluster sequences with --cluster-mode 2 to prefer longest sequences.

    Args:
        input_fasta: Path to input FASTA file (e.g., repseq.faa)
        output_dir: Directory for output files
        db_path: Optional path to DuckDB to update with new cluster info
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    tmp_dir = os.path.join(output_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    logger.info(f"Starting re-clustering with cluster-mode 2 (longest sequence as rep)")
    logger.info(f"Input: {input_fasta}")
    logger.info(f"Output: {output_dir}")

    # Define paths
    sequence_db = os.path.join(output_dir, "sequenceDB")
    cluster_db = os.path.join(output_dir, "clusterDB")
    rep_seq_db = os.path.join(output_dir, "repSeqDB")
    cluster_tsv = os.path.join(output_dir, "clusterDB.tsv")
    new_repseq_fasta = os.path.join(output_dir, "repseq_longest.faa")

    try:
        # Step 1: Create MMSeqs2 database from FASTA
        run_command(
            ["mmseqs", "createdb", input_fasta, sequence_db],
            "Creating MMSeqs2 sequence database"
        )

        # Step 2: Run clustering with --cluster-mode 2 (longest sequence as representative)
        # Also adding common parameters for better control
        run_command(
            [
                "mmseqs", "cluster",
                sequence_db,
                cluster_db,
                tmp_dir,
                "--cluster-mode", "2",      # Longest sequence as representative
                "-v", "3"                     # Verbose output
            ],
            "Running MMSeqs2 clustering (cluster-mode 2: longest sequence as rep)"
        )

        # Step 3: Create TSV output
        run_command(
            [
                "mmseqs", "createtsv",
                sequence_db,
                sequence_db,
                cluster_db,
                cluster_tsv
            ],
            "Creating cluster TSV file"
        )

        # Step 4: Extract representative sequences
        run_command(
            [
                "mmseqs", "result2repseq",
                sequence_db,
                cluster_db,
                rep_seq_db
            ],
            "Extracting representative sequences"
        )

        # Step 5: Convert to FASTA
        run_command(
            [
                "mmseqs", "convert2fasta",
                rep_seq_db,
                new_repseq_fasta
            ],
            "Converting representatives to FASTA"
        )

        # Step 6: Report statistics
        logger.info("\n" + "="*80)
        logger.info("CLUSTERING STATISTICS")
        logger.info("="*80)

        # Count original sequences
        orig_count = int(subprocess.check_output(
            f"grep -c '^>' {input_fasta}",
            shell=True
        ).decode().strip())

        # Count clusters (representatives)
        cluster_count = int(subprocess.check_output(
            f"cut -f1 {cluster_tsv} | sort | uniq | wc -l",
            shell=True
        ).decode().strip())

        # Count new representatives
        new_rep_count = int(subprocess.check_output(
            f"grep -c '^>' {new_repseq_fasta}",
            shell=True
        ).decode().strip())

        logger.info(f"Original sequences: {orig_count}")
        logger.info(f"Number of clusters: {cluster_count}")
        logger.info(f"New representatives: {new_rep_count}")
        logger.info(f"Reduction: {orig_count - new_rep_count} sequences ({100*(1-new_rep_count/orig_count):.1f}%)")
        logger.info("="*80 + "\n")

        # Step 7: Update DuckDB if provided
        if db_path:
            logger.info(f"Updating DuckDB: {db_path}")
            try:
                from planter.database.utils.duckdb_utils import update_clusters

                # Backup first
                backup_path = f"{db_path}.backup_before_recluster"
                shutil.copy2(db_path, backup_path)
                logger.info(f"Created backup: {backup_path}")

                # Update clusters
                update_clusters(
                    db_path=db_path,
                    tsv_path=cluster_tsv,
                    backup_first=False,  # We already made a backup
                    handle_duplicates="ignore"
                )

                logger.info("Successfully updated DuckDB with new cluster information")

            except Exception as e:
                logger.error(f"Failed to update DuckDB: {e}")
                logger.error("You can manually update later using the TSV file at:")
                logger.error(f"  {cluster_tsv}")

        # Step 8: Cleanup temporary files
        logger.info("Cleaning up temporary files")
        shutil.rmtree(tmp_dir)

        logger.info("\n" + "="*80)
        logger.info("RE-CLUSTERING COMPLETE!")
        logger.info("="*80)
        logger.info(f"New representative sequences: {new_repseq_fasta}")
        logger.info(f"Cluster TSV file: {cluster_tsv}")
        logger.info("="*80 + "\n")

    except Exception as e:
        logger.error(f"Error during re-clustering: {e}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Re-cluster sequences with --cluster-mode 2 to prefer longest sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Just re-cluster the FASTA
  python recluster_repseq.py --input repseq.faa --output ./reclustered

  # Re-cluster and update DuckDB
  python recluster_repseq.py \\
    --input /mnt/data4/repseq.faa \\
    --output /mnt/data4/reclustered \\
    --db /mnt/data4/master.duckdb
        """
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input FASTA file (e.g., repseq.faa)"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output directory for re-clustered results"
    )

    parser.add_argument(
        "-d", "--db",
        help="Optional: DuckDB path to update with new cluster information"
    )

    args = parser.parse_args()

    # Validate input file exists
    if not os.path.exists(args.input):
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)

    # Validate database exists if provided
    if args.db and not os.path.exists(args.db):
        logger.error(f"Database file not found: {args.db}")
        sys.exit(1)

    recluster_sequences(args.input, args.output, args.db)


if __name__ == "__main__":
    main()
