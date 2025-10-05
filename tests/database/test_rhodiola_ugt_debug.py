"""
Debug test for Rhodiola UGT search discrepancy between MMSeqs2 and NCBI BLAST.

Issue: MMSeqs2 returns truncated UGT (240 aa, AUI41142) instead of expected
full-length sequences (498 aa, AUI41147 and related hits).

Expected NCBI BLAST results: AUI41147, AUI41146, AUI41143, AUI41132, AUI41144
"""
import duckdb
import pytest
from pathlib import Path


# Test sequence from user
AUI41147_SEQUENCE = """MSLIEKPLTAIETREKPHAVCIPYPAQGHINPMMQLAKLLHHSGFHITFVHTEYNYDRLVKSQGSACVAGLPDFRFEAIPDGLPSTNGDVTQDIPLLSSSTSKTCLKPFKELLKRLQDKCKELPDDVPPLSCIVSDAAMSFTIDASEEFGVPIALLWTASACGFLGYTHYPYLIDRGVIPLKDESQLTNGYLDMSIDGIPCMEGIRLRDLPSFLRTTDLDDMMFSYILHEIKQVSRGSAIILNTFEALDHDVLDSLSKIYQNVILPVGPLHVSLNKIPKHYPLQSLSSNLWKDDTDCIPWLSSKASKSVIYVNFGSITTVSPKQIVEFAWGLANSKHPFLWIIRPDLVAGEASIIPQDFMDETKGRGLLAGWCDQELVLNHPSIGGFLTHCGWNSIIESISAGVPTVCWPFFAEQQTNCWFACKKWCIGMEMHTDVKRDEVDKLLRELMEGDKGEELKRKATNWKRLAEEAVSSTGLSTLNFRTLVNQVLLSKTKHIR"""

EXPECTED_NCBI_HITS = ['AUI41147', 'AUI41146', 'AUI41143', 'AUI41132', 'AUI41144']
UNEXPECTED_MMSEQS_HIT = 'AUI41142'


@pytest.fixture
def db_connection():
    """Connect to the master database."""
    db_path = Path("database/master.duckdb")
    if not db_path.exists():
        pytest.skip(f"Database not found at {db_path}")

    conn = duckdb.connect(str(db_path), read_only=True)
    yield conn
    conn.close()


def test_sequence_length():
    """Verify AUI41147 sequence is 498 amino acids."""
    seq = AUI41147_SEQUENCE.replace('\n', '').strip()
    assert len(seq) == 498, f"Expected 498 aa, got {len(seq)} aa"


def test_rhodiola_sequences_in_db(db_connection):
    """Query all Rhodiola rosea sequences in the database."""
    query = """
    SELECT
        seqhash,
        ncbi_accession,
        length(sequence) as seq_length,
        is_representative,
        cluster_id,
        organism
    FROM sequences
    WHERE organism LIKE '%Rhodiola rosea%'
    ORDER BY seq_length DESC
    """

    result = db_connection.execute(query).fetchdf()
    print(f"\n\nFound {len(result)} Rhodiola rosea sequences:")
    print(result.to_string())

    assert len(result) > 0, "No Rhodiola rosea sequences found in database"


def test_specific_accessions_in_db(db_connection):
    """Search for the specific accessions mentioned."""
    all_accessions = EXPECTED_NCBI_HITS + [UNEXPECTED_MMSEQS_HIT]

    query = """
    SELECT
        ncbi_accession,
        seqhash,
        length(sequence) as seq_length,
        is_representative,
        cluster_id,
        organism
    FROM sequences
    WHERE ncbi_accession IN ({})
    ORDER BY seq_length DESC
    """.format(','.join([f"'{acc}'" for acc in all_accessions]))

    result = db_connection.execute(query).fetchdf()
    print(f"\n\nSearching for specific accessions:")
    print(result.to_string())

    found_accessions = set(result['ncbi_accession'].tolist())
    print(f"\n\nExpected NCBI hits: {EXPECTED_NCBI_HITS}")
    print(f"Found in DB: {found_accessions}")
    print(f"Missing: {set(EXPECTED_NCBI_HITS) - found_accessions}")


def test_sequence_length_distribution(db_connection):
    """Check the length distribution of Rhodiola UGT sequences."""
    query = """
    SELECT
        ncbi_accession,
        length(sequence) as seq_length,
        is_representative,
        cluster_id
    FROM sequences
    WHERE organism LIKE '%Rhodiola rosea%'
    AND ncbi_accession IN ({})
    ORDER BY seq_length DESC
    """.format(','.join([f"'{acc}'" for acc in EXPECTED_NCBI_HITS + [UNEXPECTED_MMSEQS_HIT]]))

    result = db_connection.execute(query).fetchdf()
    print(f"\n\nSequence lengths:")
    print(result.to_string())

    # Check for truncated sequences
    truncated = result[result['seq_length'] < 400]
    if len(truncated) > 0:
        print(f"\n\nWARNING: Found truncated sequences (<400 aa):")
        print(truncated.to_string())


def test_representative_status(db_connection):
    """Check which sequences are marked as representative."""
    all_accessions = EXPECTED_NCBI_HITS + [UNEXPECTED_MMSEQS_HIT]

    query = """
    SELECT
        ncbi_accession,
        seqhash,
        is_representative,
        cluster_id,
        length(sequence) as seq_length
    FROM sequences
    WHERE ncbi_accession IN ({})
    ORDER BY cluster_id, is_representative DESC
    """.format(','.join([f"'{acc}'" for acc in all_accessions]))

    result = db_connection.execute(query).fetchdf()
    print(f"\n\nRepresentative status:")
    print(result.to_string())

    # Check if truncated sequence is marked as representative
    if UNEXPECTED_MMSEQS_HIT in result['ncbi_accession'].values:
        aui41142 = result[result['ncbi_accession'] == UNEXPECTED_MMSEQS_HIT]
        if not aui41142.empty and aui41142.iloc[0]['is_representative']:
            print(f"\n\nWARNING: {UNEXPECTED_MMSEQS_HIT} (truncated) is marked as representative!")


def test_cluster_membership(db_connection):
    """Check cluster membership for all related sequences."""
    all_accessions = EXPECTED_NCBI_HITS + [UNEXPECTED_MMSEQS_HIT]

    query = """
    SELECT
        cluster_id,
        COUNT(*) as cluster_size,
        GROUP_CONCAT(ncbi_accession, ', ') as members,
        SUM(CASE WHEN is_representative THEN 1 ELSE 0 END) as num_reps
    FROM sequences
    WHERE ncbi_accession IN ({})
    GROUP BY cluster_id
    ORDER BY cluster_size DESC
    """.format(','.join([f"'{acc}'" for acc in all_accessions]))

    result = db_connection.execute(query).fetchdf()
    print(f"\n\nCluster membership:")
    print(result.to_string())
