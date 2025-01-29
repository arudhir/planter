#!/usr/bin/env python3
"""
Simple script to create a DuckDB database for a 
single sample using Planter.
"""

import argparse
import logging
from pathlib import Path
import sys
sys.path.append(str(Path(__file__).parent.parent))
from planter.database.builder import SequenceDBBuilder

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Create a DuckDB database for a genomic sample')
    parser.add_argument(
        '--sample_id', 
        help='Transcriptomic SRA ID to process', 
        required=True
    )
    return parser.parse_args()

def create_duckdb(sample_id: str):
    """Create DuckDB database for the given sample ID."""
    duckdb_path = Path(f"/mnt/data3/planter_outputs/{sample_id}/{sample_id}.duckdb")
    base_output_dir = Path(f"/mnt/data3/planter_outputs/")

    logger.info(f"Processing sample {sample_id}")
    with SequenceDBBuilder(duckdb_path, output_dir=base_output_dir) as builder:
        results = builder.build_database([sample_id])
        logger.info(f"Build results: {results}")
        
        summary = builder.get_database_summary()
        logger.info(f"Database summary: {summary}")

def main():
    args = parse_args()
    create_duckdb(args.sample_id)

if __name__ == "__main__":
    main()