#!/usr/bin/env python3
"""
Rebuild all sample DuckDB databases with the correct schema.

This will delete and recreate each sample's DuckDB file using the
correct schema where repseq_id is nullable.
"""
import sys
from pathlib import Path
from tqdm import tqdm
import logging

from planter.database.builder import SequenceDBBuilder

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def rebuild_sample_duckdb(sample_id: str, base_dir: Path) -> bool:
    """
    Rebuild DuckDB for a single sample.

    Returns True if successful, False otherwise.
    """
    sample_dir = base_dir / sample_id
    duckdb_path = sample_dir / f"{sample_id}.duckdb"

    # Check if sample directory exists
    if not sample_dir.exists():
        print(f"  ⚠ {sample_id}: Sample directory not found, skipping")
        return True

    # Delete old DuckDB and WAL files if they exist
    try:
        if duckdb_path.exists():
            duckdb_path.unlink()

        wal_path = Path(f"{duckdb_path}.wal")
        if wal_path.exists():
            wal_path.unlink()

    except Exception as e:
        print(f"  ✗ {sample_id}: Error deleting old database: {e}")
        return False

    # Rebuild the database
    try:
        with SequenceDBBuilder(str(duckdb_path), output_dir=str(sample_dir)) as builder:
            results = builder.build_database([sample_id])

            if results and sample_id in results:
                return True
            else:
                print(f"  ✗ {sample_id}: Build returned no results")
                return False

    except Exception as e:
        print(f"  ✗ {sample_id}: Error building database: {e}")
        return False


def main():
    # Read sample list
    samples_file = Path('/tmp/samples_with_duckdb.txt')

    if not samples_file.exists():
        print(f"Error: {samples_file} not found")
        print("This should contain the list of samples to rebuild.")
        sys.exit(1)

    with open(samples_file) as f:
        samples = [line.strip() for line in f if line.strip()]

    print(f"Rebuilding DuckDB databases for {len(samples)} samples...")
    print("This will delete and recreate each sample's DuckDB file.\n")

    response = input("Continue? (yes/no): ")
    if response.lower() != 'yes':
        print("Aborted.")
        sys.exit(0)

    base_dir = Path('/mnt/data4/recombia.planter')
    success_count = 0
    fail_count = 0

    for sample in tqdm(samples, desc="Rebuilding databases"):
        if rebuild_sample_duckdb(sample, base_dir):
            success_count += 1
        else:
            fail_count += 1

    print(f"\n{'='*60}")
    print(f"✓ Successfully rebuilt: {success_count} databases")
    if fail_count > 0:
        print(f"✗ Failed: {fail_count} databases")
    print(f"{'='*60}\n")

    if fail_count > 0:
        print("⚠ Some databases failed to rebuild. Check the errors above.")
        sys.exit(1)
    else:
        print("✓ All databases rebuilt successfully!")
        print("\nYou can now run: python rebuild_master_db.py")


if __name__ == '__main__':
    main()
