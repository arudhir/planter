#!/usr/bin/env python3
"""
Simple script to create a DuckDB database for a
single sample using Planter.
"""

import argparse
import logging
import sys
from pathlib import Path

from planter.database.builder import SequenceDBBuilder

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Create a DuckDB database for a genomic sample"
    )
    parser.add_argument(
        "--sample_id", help="Transcriptomic SRA ID to process", required=True
    )

    parser.add_argument("--outdir", help="Output directory", required=True)

    parser.add_argument("--duckdb_out", help="duckdb output file", required=True)

    return parser.parse_args()


def create_duckdb(sample_id: str, outdir: str, duckdb_out: str):
    """Create DuckDB database for the given sample ID."""

    logger.info(f"Processing sample {sample_id}")
    with SequenceDBBuilder(duckdb_out, output_dir=outdir) as builder:
        results = builder.build_database([sample_id])
        logger.info(f"Build results: {results}")

        summary = builder.get_database_summary()
        logger.info(f"Database summary: {summary}")


def main():
    args = parse_args()
    create_duckdb(
        sample_id=args.sample_id, outdir=args.outdir, duckdb_out=args.duckdb_out
    )


if __name__ == "__main__":
    main()
