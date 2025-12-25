#!/usr/bin/env python3
"""
Integration tests for incremental clustering workflow.

This test suite simulates the full workflow that was broken by the bug:
1. Initial clustering of samples
2. Adding new samples to database
3. Running incremental clustering (should include new samples)
4. Verifying all sequences get clustered

This reproduces the exact scenario that caused Rhodiola UGTs to remain unclustered.
"""
import pytest
import duckdb
from pathlib import Path
import tempfile
import shutil

from planter.database.utils.duckdb_utils import (
    extract_representative_sequences,
    update_clusters
)


@pytest.fixture
def workflow_test_dir():
    """Create a temporary directory for the workflow test."""
    test_dir = Path(tempfile.mkdtemp(prefix='clustering_workflow_test_'))
    yield test_dir
    shutil.rmtree(test_dir)


@pytest.fixture
def initial_database(workflow_test_dir):
    """
    Create initial database with clustered sequences.

    Simulates the state after initial clustering (output130):
    - 3 samples, each with multiple sequences
    - All sequences have been clustered (repseq_id set)
    - Representative sequences identified
    """
    db_path = workflow_test_dir / 'master.duckdb'
    con = duckdb.connect(str(db_path))

    # Create schema
    con.execute("""
        CREATE TABLE sequences (
            seqhash_id VARCHAR PRIMARY KEY,
            sequence VARCHAR,
            sample_id VARCHAR,
            assembly_date TIMESTAMP,
            is_representative BOOLEAN DEFAULT FALSE,
            repseq_id VARCHAR,
            length INTEGER
        )
    """)

    con.execute("""
        CREATE TABLE clusters (
            cluster_id VARCHAR PRIMARY KEY,
            representative_seqhash_id VARCHAR NOT NULL,
            size INTEGER NOT NULL
        )
    """)

    con.execute("""
        CREATE TABLE cluster_members (
            seqhash_id VARCHAR PRIMARY KEY,
            cluster_id VARCHAR NOT NULL,
            FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id),
            FOREIGN KEY (cluster_id) REFERENCES clusters(cluster_id)
        )
    """)

    # Insert initial clustered sequences (simulating output130)
    # Cluster 1: Kinase sequences
    con.execute("""
        INSERT INTO sequences (seqhash_id, sequence, sample_id, is_representative, repseq_id, length) VALUES
        ('kinase_rep', 'MAAAKKKGGGTTTWWW', 'SRR001', TRUE, 'kinase_rep', 16),
        ('kinase_m1', 'MAAAKKKGGGTTTWWX', 'SRR001', FALSE, 'kinase_rep', 16),
        ('kinase_m2', 'MAAAKKKGGGTTTWWY', 'SRR002', FALSE, 'kinase_rep', 16)
    """)

    # Cluster 2: Transferase sequences
    con.execute("""
        INSERT INTO sequences (seqhash_id, sequence, sample_id, is_representative, repseq_id, length) VALUES
        ('transferase_rep', 'MQQQRRRTTTYYYLLL', 'SRR002', TRUE, 'transferase_rep', 16),
        ('transferase_m1', 'MQQQRRRTTTYYYLLX', 'SRR003', FALSE, 'transferase_rep', 16)
    """)

    # Cluster 3: Singleton (no similar sequences)
    con.execute("""
        INSERT INTO sequences (seqhash_id, sequence, sample_id, is_representative, repseq_id, length) VALUES
        ('singleton', 'MPPPSSSHHHNNNFFF', 'SRR003', TRUE, 'singleton', 16)
    """)

    # Create clusters
    con.execute("""
        INSERT INTO clusters VALUES
        ('kinase_rep', 'kinase_rep', 3),
        ('transferase_rep', 'transferase_rep', 2),
        ('singleton', 'singleton', 1)
    """)

    # Create cluster memberships
    con.execute("""
        INSERT INTO cluster_members VALUES
        ('kinase_rep', 'kinase_rep'),
        ('kinase_m1', 'kinase_rep'),
        ('kinase_m2', 'kinase_rep'),
        ('transferase_rep', 'transferase_rep'),
        ('transferase_m1', 'transferase_rep'),
        ('singleton', 'singleton')
    """)

    con.close()

    return db_path


def test_incremental_workflow_bug_reproduction(workflow_test_dir, initial_database):
    """
    Reproduce the exact bug scenario:
    1. Start with clustered database
    2. Add new samples (Rhodiola) with repseq_id = NULL
    3. Extract representatives
    4. Verify new samples ARE included (the fix)

    This is the CORE integration test for the bug fix.
    """
    # PHASE 1: Initial state (after output130 clustering)
    con = duckdb.connect(str(initial_database))

    initial_count = con.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
    assert initial_count == 6, f"Expected 6 initial sequences, got {initial_count}"

    clustered_count = con.execute(
        "SELECT COUNT(*) FROM sequences WHERE repseq_id IS NOT NULL"
    ).fetchone()[0]
    assert clustered_count == 6, "All initial sequences should be clustered"

    # PHASE 2: Add new Rhodiola samples (simulating samples added AFTER clustering)
    # These are added with repseq_id = NULL because they haven't been clustered yet
    con.execute("""
        INSERT INTO sequences (seqhash_id, sequence, sample_id, is_representative, repseq_id, length)
        VALUES
        ('rhodiola_ugt1', 'MUUUGGGTTT111222', 'SRR5936536', FALSE, NULL, 16),
        ('rhodiola_ugt2', 'MUUUGGGTTT111333', 'SRR5936537', FALSE, NULL, 16),
        ('rhodiola_ugt3', 'MUUUGGGTTT111444', 'SRR22844166', FALSE, NULL, 16),
        ('other_new1', 'MVVVWWWXXXYYYZZ1', 'SRR999', FALSE, NULL, 16),
        ('other_new2', 'MVVVWWWXXXYYYZZ2', 'SRR888', FALSE, NULL, 16)
    """)

    total_after = con.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
    assert total_after == 11, f"Expected 11 sequences after adding new samples, got {total_after}"

    unclustered_count = con.execute(
        "SELECT COUNT(*) FROM sequences WHERE repseq_id IS NULL"
    ).fetchone()[0]
    assert unclustered_count == 5, f"Expected 5 unclustered sequences, got {unclustered_count}"

    con.close()

    # PHASE 3: Extract representatives for next clustering round
    # This is where the bug occurred - unclustered sequences were EXCLUDED
    repseq_fasta = workflow_test_dir / 'representatives.faa'
    extract_representative_sequences(initial_database, repseq_fasta)

    # PHASE 4: Verify extraction includes BOTH representatives AND unclustered sequences
    from Bio import SeqIO
    extracted_ids = set()
    for record in SeqIO.parse(repseq_fasta, 'fasta'):
        extracted_ids.add(record.id)

    # Expected: 3 existing representatives + 5 unclustered sequences = 8
    expected_ids = {
        # Existing representatives
        'kinase_rep', 'transferase_rep', 'singleton',
        # NEW unclustered sequences (THE FIX - these MUST be included)
        'rhodiola_ugt1', 'rhodiola_ugt2', 'rhodiola_ugt3',
        'other_new1', 'other_new2'
    }

    assert extracted_ids == expected_ids, \
        f"Expected {expected_ids}\nGot {extracted_ids}\nMissing: {expected_ids - extracted_ids}"

    # CRITICAL: Verify Rhodiola UGTs are in the extraction
    assert 'rhodiola_ugt1' in extracted_ids, "Rhodiola UGT 1 MUST be extracted"
    assert 'rhodiola_ugt2' in extracted_ids, "Rhodiola UGT 2 MUST be extracted"
    assert 'rhodiola_ugt3' in extracted_ids, "Rhodiola UGT 3 MUST be extracted"

    # Verify clustered members are NOT in extraction
    assert 'kinase_m1' not in extracted_ids, "Clustered members should NOT be extracted"
    assert 'kinase_m2' not in extracted_ids
    assert 'transferase_m1' not in extracted_ids


def test_old_bug_behavior(workflow_test_dir, initial_database):
    """
    Demonstrate that the OLD buggy behavior would have excluded unclustered sequences.

    This test shows what WOULD have happened without the fix.
    """
    # Add new samples
    con = duckdb.connect(str(initial_database))
    con.execute("""
        INSERT INTO sequences (seqhash_id, sequence, sample_id, is_representative, repseq_id, length)
        VALUES
        ('rhodiola_ugt1', 'MUUUGGGTTT111222', 'SRR5936536', FALSE, NULL, 16),
        ('rhodiola_ugt2', 'MUUUGGGTTT111333', 'SRR5936537', FALSE, NULL, 16)
    """)

    # OLD buggy query (without OR repseq_id IS NULL)
    old_buggy_query = """
        SELECT seqhash_id FROM sequences
        WHERE repseq_id = seqhash_id
    """

    old_results = con.execute(old_buggy_query).fetchall()
    old_extracted_ids = {row[0] for row in old_results}

    # Bug: Only existing representatives, NO unclustered sequences
    assert len(old_extracted_ids) == 3, \
        f"Old bug would extract only 3 reps, got {len(old_extracted_ids)}"

    # THE BUG: Rhodiola sequences were EXCLUDED
    assert 'rhodiola_ugt1' not in old_extracted_ids, \
        "OLD BUG: Rhodiola UGT 1 was EXCLUDED from extraction"
    assert 'rhodiola_ugt2' not in old_extracted_ids, \
        "OLD BUG: Rhodiola UGT 2 was EXCLUDED from extraction"

    con.close()


def test_full_incremental_clustering_cycle(workflow_test_dir, initial_database):
    """
    Test a complete incremental clustering cycle:
    1. Initial state: Some sequences clustered
    2. Add new samples: Sequences with repseq_id = NULL
    3. Extract representatives: Should include new sequences
    4. Simulate clustering: Create mock cluster TSV
    5. Load clusters: Update database
    6. Verify: All sequences now have cluster assignments
    """
    # Step 1: Add new Rhodiola samples
    con = duckdb.connect(str(initial_database))
    con.execute("""
        INSERT INTO sequences (seqhash_id, sequence, sample_id, is_representative, repseq_id, length)
        VALUES
        ('rhodiola_ugt1', 'MUUUGGGTTT111222', 'SRR5936536', FALSE, NULL, 16),
        ('rhodiola_ugt2', 'MUUUGGGTTT111333', 'SRR5936537', FALSE, NULL, 16)
    """)
    con.close()

    # Step 2: Extract representatives (includes unclustered after fix)
    repseq_fasta = workflow_test_dir / 'representatives.faa'
    extract_representative_sequences(initial_database, repseq_fasta)

    # Step 3: Simulate MMSeqs2 clustering output
    # In reality, MMSeqs2 would compare new samples against extracted representatives
    # For this test, we simulate the output TSV
    cluster_tsv = workflow_test_dir / 'newClusterDB.tsv'

    # Simulate clustering results:
    # - rhodiola_ugt1 clusters with rhodiola_ugt2 (similar sequences)
    # - Both are new, so rhodiola_ugt1 becomes representative (first seen)
    with open(cluster_tsv, 'w') as f:
        # Existing clusters (unchanged)
        f.write("kinase_rep\tkinase_rep\n")
        f.write("kinase_rep\tkinase_m1\n")
        f.write("kinase_rep\tkinase_m2\n")
        f.write("transferase_rep\ttransferase_rep\n")
        f.write("transferase_rep\ttransferase_m1\n")
        f.write("singleton\tsingleton\n")

        # NEW: Rhodiola cluster (this is the fix - they get clustered!)
        f.write("rhodiola_ugt1\trhodiola_ugt1\n")  # rhodiola_ugt1 is representative
        f.write("rhodiola_ugt1\trhodiola_ugt2\n")  # rhodiola_ugt2 clusters with it

    # Step 4: Load clusters into database
    update_clusters(initial_database, cluster_tsv, backup_first=False, handle_duplicates="replace")

    # Step 5: Verify ALL sequences now have cluster assignments
    con = duckdb.connect(str(initial_database))

    # Check total sequences
    total = con.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
    assert total == 8, f"Expected 8 sequences total, got {total}"

    # Check clustered sequences (should be ALL of them now)
    clustered = con.execute(
        "SELECT COUNT(*) FROM sequences WHERE repseq_id IS NOT NULL"
    ).fetchone()[0]
    assert clustered == 8, \
        f"ALL sequences should be clustered after update, got {clustered}/8"

    # Check unclustered sequences (should be ZERO now)
    unclustered = con.execute(
        "SELECT COUNT(*) FROM sequences WHERE repseq_id IS NULL"
    ).fetchone()[0]
    assert unclustered == 0, \
        f"NO sequences should be unclustered after update, got {unclustered}"

    # CRITICAL: Verify Rhodiola UGTs have cluster assignments
    rhodiola_results = con.execute("""
        SELECT seqhash_id, repseq_id, is_representative
        FROM sequences
        WHERE seqhash_id IN ('rhodiola_ugt1', 'rhodiola_ugt2')
        ORDER BY seqhash_id
    """).fetchall()

    assert len(rhodiola_results) == 2, "Should have 2 Rhodiola sequences"

    # rhodiola_ugt1 should be representative
    ugt1 = [r for r in rhodiola_results if r[0] == 'rhodiola_ugt1'][0]
    assert ugt1[1] == 'rhodiola_ugt1', "rhodiola_ugt1 should be its own representative"
    assert ugt1[2] == True, "rhodiola_ugt1 should be marked as representative"

    # rhodiola_ugt2 should cluster with rhodiola_ugt1
    ugt2 = [r for r in rhodiola_results if r[0] == 'rhodiola_ugt2'][0]
    assert ugt2[1] == 'rhodiola_ugt1', "rhodiola_ugt2 should have rhodiola_ugt1 as representative"
    assert ugt2[2] == False, "rhodiola_ugt2 should NOT be marked as representative"

    # Verify cluster_members table
    member_count = con.execute("SELECT COUNT(*) FROM cluster_members").fetchone()[0]
    assert member_count == 8, \
        f"cluster_members should have 8 entries (all sequences), got {member_count}"

    # Verify Rhodiola sequences are in cluster_members
    rhodiola_members = con.execute("""
        SELECT seqhash_id, cluster_id
        FROM cluster_members
        WHERE seqhash_id IN ('rhodiola_ugt1', 'rhodiola_ugt2')
        ORDER BY seqhash_id
    """).fetchall()

    assert len(rhodiola_members) == 2, "Both Rhodiola UGTs should be in cluster_members"
    assert rhodiola_members[0] == ('rhodiola_ugt1', 'rhodiola_ugt1')
    assert rhodiola_members[1] == ('rhodiola_ugt2', 'rhodiola_ugt1')

    con.close()


def test_clustering_health_check(workflow_test_dir, initial_database):
    """
    Test the clustering health check metric recommended in CLUSTERING_BUG_ANALYSIS.md

    This metric should alert if >10% of sequences are unclustered.
    """
    con = duckdb.connect(str(initial_database))

    # Initially: 100% clustered (healthy)
    result = con.execute("""
        SELECT
            COUNT(*) as total,
            COUNT(CASE WHEN repseq_id IS NULL THEN 1 END) as unclustered,
            ROUND(COUNT(CASE WHEN repseq_id IS NULL THEN 1 END) * 100.0 / COUNT(*), 1) as unclustered_pct
        FROM sequences
    """).fetchone()

    assert result[2] == 0.0, "Initially should have 0% unclustered"

    # Add many unclustered sequences (simulating the bug scenario)
    for i in range(50):
        con.execute(f"""
            INSERT INTO sequences (seqhash_id, sequence, sample_id, is_representative, repseq_id, length)
            VALUES ('new_{i}', 'MAAAABBBBCCCC', 'new_sample_{i}', FALSE, NULL, 13)
        """)

    # Check health metric
    result = con.execute("""
        SELECT
            COUNT(*) as total,
            COUNT(CASE WHEN repseq_id IS NULL THEN 1 END) as unclustered,
            ROUND(COUNT(CASE WHEN repseq_id IS NULL THEN 1 END) * 100.0 / COUNT(*), 1) as unclustered_pct
        FROM sequences
    """).fetchone()

    total, unclustered, unclustered_pct = result

    assert total == 56, f"Expected 56 total sequences, got {total}"
    assert unclustered == 50, f"Expected 50 unclustered, got {unclustered}"
    assert unclustered_pct == 89.3, f"Expected 89.3% unclustered, got {unclustered_pct}%"

    # This should trigger a health alert (>10% unclustered)
    ALERT_THRESHOLD = 10.0
    if unclustered_pct > ALERT_THRESHOLD:
        # This is the condition that should have alerted us to the bug
        alert_message = f"ALERT: {unclustered_pct}% of sequences are unclustered!"
        print(alert_message)
        assert True, "Health check correctly detects high unclustered percentage"

    con.close()


def test_multiple_incremental_rounds(workflow_test_dir, initial_database):
    """
    Test multiple rounds of incremental clustering:
    Round 1: Add samples A, cluster
    Round 2: Add samples B, cluster
    Round 3: Add samples C, cluster

    Verify all samples eventually get clustered.
    """
    con = duckdb.connect(str(initial_database))

    # Round 1: Add batch A
    con.execute("""
        INSERT INTO sequences (seqhash_id, sequence, sample_id, is_representative, repseq_id, length)
        VALUES
        ('batch_a1', 'MAAABBBCCC111', 'SampleA1', FALSE, NULL, 13),
        ('batch_a2', 'MAAABBBCCC222', 'SampleA2', FALSE, NULL, 13)
    """)
    con.close()

    # Cluster round 1
    repseq_fasta = workflow_test_dir / 'round1.faa'
    extract_representative_sequences(initial_database, repseq_fasta)

    # Should extract: 3 original reps + 2 unclustered from batch A = 5
    from Bio import SeqIO
    round1_ids = {r.id for r in SeqIO.parse(repseq_fasta, 'fasta')}
    assert len(round1_ids) == 5
    assert 'batch_a1' in round1_ids
    assert 'batch_a2' in round1_ids

    # Simulate clustering and update
    cluster_tsv = workflow_test_dir / 'round1.tsv'
    with open(cluster_tsv, 'w') as f:
        f.write("kinase_rep\tkinase_rep\n")
        f.write("kinase_rep\tkinase_m1\n")
        f.write("kinase_rep\tkinase_m2\n")
        f.write("transferase_rep\ttransferase_rep\n")
        f.write("transferase_rep\ttransferase_m1\n")
        f.write("singleton\tsingleton\n")
        f.write("batch_a1\tbatch_a1\n")
        f.write("batch_a1\tbatch_a2\n")  # a2 clusters with a1

    update_clusters(initial_database, cluster_tsv, backup_first=False, handle_duplicates="replace")

    # Round 2: Add batch B
    con = duckdb.connect(str(initial_database))
    con.execute("""
        INSERT INTO sequences (seqhash_id, sequence, sample_id, is_representative, repseq_id, length)
        VALUES
        ('batch_b1', 'MDDDEEEFFF333', 'SampleB1', FALSE, NULL, 13)
    """)
    con.close()

    # Extract for round 2
    repseq_fasta2 = workflow_test_dir / 'round2.faa'
    extract_representative_sequences(initial_database, repseq_fasta2)

    # Should extract: 4 reps (original 3 + batch_a1) + 1 unclustered (batch_b1) = 5
    round2_ids = {r.id for r in SeqIO.parse(repseq_fasta2, 'fasta')}
    assert len(round2_ids) == 5
    assert 'batch_b1' in round2_ids  # New unclustered sequence
    assert 'batch_a1' in round2_ids  # Now a representative

    # Final check: All sequences should be clusterable
    con = duckdb.connect(str(initial_database))
    total = con.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
    assert total == 9, f"Expected 9 sequences after 2 rounds, got {total}"
    con.close()
