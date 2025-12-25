#!/usr/bin/env python3
"""
Create DuckDB files for samples that have .pep files but no .duckdb
"""
import sys
from pathlib import Path
import logging

sys.path.insert(0, str(Path(__file__).parent))

from planter.database.builder import SequenceDBBuilder

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def create_sample_db(sample_id: str, base_dir: Path):
    """Create DuckDB for a single sample"""
    sample_dir = base_dir / sample_id

    # Check required files exist
    pep_file = sample_dir / "transdecoder" / f"{sample_id}.pep"
    emapper_file = sample_dir / "eggnog" / f"{sample_id}.emapper.annotations"

    if not pep_file.exists():
        logging.warning(f"Skipping {sample_id}: .pep file not found at {pep_file}")
        return False

    if not emapper_file.exists():
        logging.warning(f"Skipping {sample_id}: emapper file not found at {emapper_file}")
        return False

    # Output database path
    db_path = sample_dir / f"{sample_id}.duckdb"

    if db_path.exists():
        logging.info(f"Skipping {sample_id}: database already exists")
        return True

    logging.info(f"Creating database for {sample_id}...")

    try:
        # Create database builder
        with SequenceDBBuilder(str(db_path), output_dir=str(base_dir)) as builder:
            # Build database for this sample (loads sequences, annotations, expression)
            results = builder.build_database([sample_id])
            logging.info(f"  Build results: {results}")

        logging.info(f"✓ Created {db_path}")
        return True

    except Exception as e:
        logging.error(f"✗ Failed to create database for {sample_id}: {e}")
        # Clean up partial database
        if db_path.exists():
            db_path.unlink()
        return False


def main():
    # Read list of samples
    samples_file = Path('/tmp/missing_duckdb_samples.txt')

    if not samples_file.exists():
        print("Error: /tmp/missing_duckdb_samples.txt not found")
        sys.exit(1)

    with open(samples_file) as f:
        samples = [line.strip() for line in f if line.strip()]

    base_dir = Path('/mnt/data4/recombia.planter')

    logging.info(f"Creating DuckDB files for {len(samples)} samples...")

    success_count = 0
    failed_samples = []

    for i, sample in enumerate(samples, 1):
        logging.info(f"\n[{i}/{len(samples)}] Processing {sample}")

        if create_sample_db(sample, base_dir):
            success_count += 1
        else:
            failed_samples.append(sample)

    logging.info(f"\n{'='*80}")
    logging.info(f"Summary:")
    logging.info(f"  Successfully created: {success_count}/{len(samples)}")

    if failed_samples:
        logging.warning(f"  Failed samples ({len(failed_samples)}):")
        for sample in failed_samples:
            logging.warning(f"    - {sample}")

    logging.info(f"{'='*80}")

    # Update the samples list
    if success_count > 0:
        logging.info("\nNext steps:")
        logging.info("1. Find all samples with DuckDB files:")
        logging.info("   find /mnt/data4/recombia.planter -name '*.pep' -path '*/transdecoder/*.pep' | grep -v transdecoder_dir | cut -d'/' -f5 | sort > /tmp/completed_samples.txt")
        logging.info("   while read sample; do if [ -f \"/mnt/data4/recombia.planter/$sample/$sample.duckdb\" ]; then echo \"$sample\"; fi; done < /tmp/completed_samples.txt > /tmp/all_samples_with_duckdb.txt")
        logging.info("\n2. Re-merge the database:")
        logging.info("   python smart_merge.py  # (after updating it to use /tmp/all_samples_with_duckdb.txt)")
        logging.info("\n3. Load clusters:")
        logging.info("   python scripts/load_clusters.py -d /mnt/data4/master_complete.duckdb -t /mnt/data4/planter_outputs/repseq_v3/output130/newClusterDB.tsv")


if __name__ == '__main__':
    main()
