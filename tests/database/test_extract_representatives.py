#!/usr/bin/env python3
"""
Unit tests for extract_representative_sequences() fix.

This test suite verifies that the bug fix correctly includes unclustered sequences
(repseq_id IS NULL) in the extraction for incremental clustering.
"""
import pytest
import duckdb
from pathlib import Path
import tempfile
from Bio import SeqIO

from planter.database.utils.duckdb_utils import extract_representative_sequences


@pytest.fixture
def test_db():
    """Create a test database with a mix of clustered and unclustered sequences."""
    with tempfile.NamedTemporaryFile(suffix='.duckdb', delete=False) as f:
        db_path = Path(f.name)

    # Remove the empty file created by NamedTemporaryFile
    db_path.unlink()

    con = duckdb.connect(str(db_path))

    # Create sequences table
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

    # Insert test data:
    # - 3 existing representatives (repseq_id = seqhash_id) from initial clustering
    # - 2 clustered members (repseq_id points to a representative)
    # - 5 unclustered sequences (repseq_id IS NULL) from new samples

    test_sequences = [
        # Existing representatives (from initial clustering)
        ('rep1', 'MKKLLVVGGAAGGTT', 'sample1', True, 'rep1'),
        ('rep2', 'MAAQQWWERRTTYYUU', 'sample2', True, 'rep2'),
        ('rep3', 'MGGGHHHJJJKKKLLL', 'sample3', True, 'rep3'),

        # Clustered members (from initial clustering)
        ('member1', 'MKKLLVVGGAAGGTX', 'sample1', False, 'rep1'),  # clusters with rep1
        ('member2', 'MAAQQWWERRTTYYUX', 'sample2', False, 'rep2'),  # clusters with rep2

        # Unclustered sequences (NEW - added after initial clustering)
        ('new1', 'MPPPQQQQRRRRSSSS', 'rhodiola1', False, None),
        ('new2', 'MTTTTTTUUUUUVVVV', 'rhodiola2', False, None),
        ('new3', 'MWWWWWWXXXXXYYYY', 'rhodiola3', False, None),
        ('new4', 'MZZZZZAAAABBBBCC', 'sample4', False, None),
        ('new5', 'MDDDDDEEEEEFFFFF', 'sample5', False, None),
    ]

    for seq_id, sequence, sample_id, is_rep, repseq_id in test_sequences:
        length = len(sequence)
        repseq_sql = f"'{repseq_id}'" if repseq_id else 'NULL'
        con.execute(f"""
            INSERT INTO sequences (seqhash_id, sequence, sample_id, is_representative, repseq_id, length)
            VALUES ('{seq_id}', '{sequence}', '{sample_id}', {is_rep}, {repseq_sql}, {length})
        """)

    con.close()

    yield db_path

    # Cleanup
    db_path.unlink()


def test_extract_includes_unclustered_sequences(test_db):
    """
    Test that extract_representative_sequences() includes BOTH:
    1. Existing representatives (repseq_id = seqhash_id)
    2. Unclustered sequences (repseq_id IS NULL)

    This is the PRIMARY test for the bug fix.
    """
    with tempfile.NamedTemporaryFile(suffix='.faa', delete=False) as f:
        output_fasta = Path(f.name)

    try:
        # Extract representatives
        extract_representative_sequences(test_db, output_fasta)

        # Parse the output FASTA
        extracted_ids = set()
        for record in SeqIO.parse(output_fasta, 'fasta'):
            extracted_ids.add(record.id)

        # Expected: 3 existing representatives + 5 unclustered sequences = 8 total
        expected_ids = {
            # Existing representatives
            'rep1', 'rep2', 'rep3',
            # Unclustered sequences (THE FIX)
            'new1', 'new2', 'new3', 'new4', 'new5'
        }

        assert extracted_ids == expected_ids, \
            f"Expected {expected_ids}, but got {extracted_ids}"

        # Verify count
        assert len(extracted_ids) == 8, \
            f"Expected 8 sequences (3 reps + 5 unclustered), got {len(extracted_ids)}"

    finally:
        output_fasta.unlink()


def test_extract_excludes_clustered_members(test_db):
    """
    Test that extract_representative_sequences() EXCLUDES:
    - Non-representative clustered members (repseq_id != seqhash_id AND repseq_id IS NOT NULL)

    These should NOT be in the extraction because they're already represented by their cluster rep.
    """
    with tempfile.NamedTemporaryFile(suffix='.faa', delete=False) as f:
        output_fasta = Path(f.name)

    try:
        extract_representative_sequences(test_db, output_fasta)

        extracted_ids = set()
        for record in SeqIO.parse(output_fasta, 'fasta'):
            extracted_ids.add(record.id)

        # These should NOT be in the extraction
        excluded_ids = {'member1', 'member2'}

        for excluded in excluded_ids:
            assert excluded not in extracted_ids, \
                f"Clustered member '{excluded}' should NOT be extracted"

    finally:
        output_fasta.unlink()


def test_extract_sequences_are_valid(test_db):
    """Test that extracted sequences have valid sequence data."""
    with tempfile.NamedTemporaryFile(suffix='.faa', delete=False) as f:
        output_fasta = Path(f.name)

    try:
        extract_representative_sequences(test_db, output_fasta)

        # Verify all sequences are valid
        for record in SeqIO.parse(output_fasta, 'fasta'):
            assert len(record.seq) > 0, \
                f"Sequence {record.id} has no sequence data"
            assert str(record.seq).startswith('M'), \
                f"Sequence {record.id} doesn't start with M (methionine)"

    finally:
        output_fasta.unlink()


def test_old_behavior_would_fail(test_db):
    """
    Test that demonstrates the OLD buggy behavior would have failed.

    This test shows what WOULD have happened with the old query:
    WHERE repseq_id = seqhash_id (without the OR repseq_id IS NULL)
    """
    con = duckdb.connect(str(test_db))

    # Old buggy query
    old_query = """
        SELECT seqhash_id, sequence
        FROM sequences
        WHERE repseq_id = seqhash_id;
    """

    old_results = con.execute(old_query).fetchall()
    old_ids = {row[0] for row in old_results}

    # Old behavior: only 3 representatives, missing all 5 unclustered sequences
    assert len(old_ids) == 3, \
        f"Old query should return only 3 reps, got {len(old_ids)}"
    assert old_ids == {'rep1', 'rep2', 'rep3'}, \
        f"Old query should only return existing reps"

    # This is the BUG: unclustered sequences are NOT included
    unclustered_ids = {'new1', 'new2', 'new3', 'new4', 'new5'}
    for unclustered in unclustered_ids:
        assert unclustered not in old_ids, \
            f"Old buggy behavior EXCLUDED unclustered sequence {unclustered}"

    con.close()


def test_new_behavior_includes_unclustered(test_db):
    """
    Test that demonstrates the NEW fixed behavior includes unclustered sequences.
    """
    con = duckdb.connect(str(test_db))

    # New fixed query
    new_query = """
        SELECT seqhash_id, sequence
        FROM sequences
        WHERE repseq_id = seqhash_id           -- Existing representatives
           OR repseq_id IS NULL;                -- Unclustered sequences (new samples)
    """

    new_results = con.execute(new_query).fetchall()
    new_ids = {row[0] for row in new_results}

    # New behavior: 3 representatives + 5 unclustered = 8 total
    assert len(new_ids) == 8, \
        f"New query should return 8 sequences (3 reps + 5 unclustered), got {len(new_ids)}"

    # Verify representatives are included
    assert 'rep1' in new_ids
    assert 'rep2' in new_ids
    assert 'rep3' in new_ids

    # Verify unclustered sequences are included (THE FIX)
    assert 'new1' in new_ids
    assert 'new2' in new_ids
    assert 'new3' in new_ids
    assert 'new4' in new_ids
    assert 'new5' in new_ids

    # Verify clustered members are still excluded
    assert 'member1' not in new_ids
    assert 'member2' not in new_ids

    con.close()


@pytest.fixture
def empty_db():
    """Create an empty test database (no sequences)."""
    with tempfile.NamedTemporaryFile(suffix='.duckdb', delete=False) as f:
        db_path = Path(f.name)

    # Remove the empty file created by NamedTemporaryFile
    db_path.unlink()

    con = duckdb.connect(str(db_path))
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
    con.close()

    yield db_path
    db_path.unlink()


def test_extract_empty_database(empty_db):
    """Test that extracting from an empty database doesn't crash."""
    with tempfile.NamedTemporaryFile(suffix='.faa', delete=False) as f:
        output_fasta = Path(f.name)

    try:
        extract_representative_sequences(empty_db, output_fasta)

        # Should create an empty FASTA file
        assert output_fasta.exists()
        records = list(SeqIO.parse(output_fasta, 'fasta'))
        assert len(records) == 0

    finally:
        output_fasta.unlink()


@pytest.fixture
def all_unclustered_db():
    """Create a database where ALL sequences are unclustered (simulating first clustering run)."""
    with tempfile.NamedTemporaryFile(suffix='.duckdb', delete=False) as f:
        db_path = Path(f.name)

    # Remove the empty file created by NamedTemporaryFile
    db_path.unlink()

    con = duckdb.connect(str(db_path))
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

    # All sequences have repseq_id = NULL (never been clustered)
    for i in range(10):
        con.execute(f"""
            INSERT INTO sequences (seqhash_id, sequence, sample_id, is_representative, repseq_id, length)
            VALUES ('seq{i}', 'MAAAABBBBCCCC', 'sample{i}', FALSE, NULL, 13)
        """)

    con.close()
    yield db_path
    db_path.unlink()


def test_extract_all_unclustered(all_unclustered_db):
    """Test extraction when ALL sequences are unclustered (first clustering run scenario)."""
    with tempfile.NamedTemporaryFile(suffix='.faa', delete=False) as f:
        output_fasta = Path(f.name)

    try:
        extract_representative_sequences(all_unclustered_db, output_fasta)

        # Should extract all 10 sequences
        extracted_ids = set()
        for record in SeqIO.parse(output_fasta, 'fasta'):
            extracted_ids.add(record.id)

        expected_ids = {f'seq{i}' for i in range(10)}
        assert extracted_ids == expected_ids
        assert len(extracted_ids) == 10

    finally:
        output_fasta.unlink()


def test_clustering_coverage_metric(test_db):
    """
    Test that we can measure clustering coverage correctly.

    This helps verify the health check metric mentioned in CLUSTERING_BUG_ANALYSIS.md
    """
    con = duckdb.connect(str(test_db))

    result = con.execute("""
        SELECT
            COUNT(*) as total,
            COUNT(CASE WHEN repseq_id IS NOT NULL THEN 1 END) as clustered,
            COUNT(CASE WHEN repseq_id IS NULL THEN 1 END) as unclustered,
            ROUND(COUNT(CASE WHEN repseq_id IS NOT NULL THEN 1 END) * 100.0 / COUNT(*), 1) as coverage_pct
        FROM sequences
    """).fetchone()

    total, clustered, unclustered, coverage_pct = result

    assert total == 10, f"Expected 10 total sequences, got {total}"
    assert clustered == 5, f"Expected 5 clustered sequences (3 reps + 2 members), got {clustered}"
    assert unclustered == 5, f"Expected 5 unclustered sequences, got {unclustered}"
    assert coverage_pct == 50.0, f"Expected 50% coverage, got {coverage_pct}%"

    con.close()
