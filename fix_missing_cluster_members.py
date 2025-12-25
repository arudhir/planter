#!/usr/bin/env python3
"""
Fix missing cluster_members entries for sequences that have repseq_id but are not in cluster_members table.
This handles the case where sequences were merged after clustering was done.
"""
import duckdb
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

DB_PATH = "/mnt/data4/recombia.planter/master.duckdb"

def fix_missing_cluster_members():
    """Add missing cluster_members entries for all sequences with repseq_id."""

    logger.info(f"Connecting to database: {DB_PATH}")
    con = duckdb.connect(DB_PATH)

    try:
        logger.info("Starting transaction...")
        con.execute("BEGIN TRANSACTION")

        # Count sequences that need to be added
        missing_count = con.execute("""
            SELECT COUNT(*)
            FROM sequences s
            WHERE s.repseq_id IS NOT NULL
              AND s.seqhash_id NOT IN (SELECT seqhash_id FROM cluster_members)
        """).fetchone()[0]

        logger.info(f"Found {missing_count} sequences with repseq_id missing from cluster_members")

        if missing_count == 0:
            logger.info("No missing entries found. Nothing to do!")
            con.execute("ROLLBACK")
            return

        # Create missing cluster entries first
        logger.info("Creating missing cluster entries...")
        con.execute("""
            INSERT OR IGNORE INTO clusters (cluster_id, representative_seqhash_id, size)
            SELECT DISTINCT
                repseq_id AS cluster_id,
                repseq_id AS representative_seqhash_id,
                0 AS size
            FROM sequences
            WHERE repseq_id IS NOT NULL
              AND repseq_id NOT IN (SELECT cluster_id FROM clusters)
        """)

        new_clusters = con.execute("""
            SELECT COUNT(*)
            FROM clusters
            WHERE size = 0
        """).fetchone()[0]
        logger.info(f"Created {new_clusters} new cluster entries")

        # Insert missing cluster_members
        logger.info("Inserting missing cluster_members...")
        con.execute("""
            INSERT OR IGNORE INTO cluster_members (seqhash_id, cluster_id)
            SELECT
                seqhash_id,
                repseq_id AS cluster_id
            FROM sequences
            WHERE repseq_id IS NOT NULL
              AND seqhash_id NOT IN (SELECT seqhash_id FROM cluster_members)
        """)

        # Update all cluster sizes to reflect true membership
        logger.info("Updating cluster sizes...")
        con.execute("""
            UPDATE clusters
            SET size = (
                SELECT COUNT(*)
                FROM cluster_members cm
                WHERE cm.cluster_id = clusters.cluster_id
            )
        """)

        # Get final statistics
        total_members = con.execute("SELECT COUNT(*) FROM cluster_members").fetchone()[0]
        total_clusters = con.execute("SELECT COUNT(*) FROM clusters").fetchone()[0]

        logger.info(f"Final totals: {total_clusters} clusters, {total_members} cluster members")

        # Verify the fix
        still_missing = con.execute("""
            SELECT COUNT(*)
            FROM sequences s
            WHERE s.repseq_id IS NOT NULL
              AND s.seqhash_id NOT IN (SELECT seqhash_id FROM cluster_members)
        """).fetchone()[0]

        if still_missing > 0:
            logger.warning(f"Warning: Still have {still_missing} sequences missing from cluster_members!")
        else:
            logger.info("✅ All sequences with repseq_id are now in cluster_members!")

        # Commit the transaction
        con.execute("COMMIT")
        logger.info("Successfully fixed missing cluster_members")

    except Exception as e:
        logger.error(f"Error: {e}")
        con.execute("ROLLBACK")
        raise
    finally:
        con.close()

if __name__ == "__main__":
    fix_missing_cluster_members()
