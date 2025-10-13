"""
Debug test for Rhodiola UGT search discrepancy between MMSeqs2 and NCBI BLAST.

Issue: MMSeqs2 returns truncated UGT (240 aa, AUI41142) instead of expected
full-length sequences (498 aa, AUI41147 and related hits).

Expected NCBI BLAST results: AUI41147, AUI41146, AUI41143, AUI41132, AUI41144
"""
import duckdb
import pytest
import subprocess
import tempfile
import os
import pandas as pd
from pathlib import Path
from planter.database.utils.duckdb_utils import get_fasta_path


# Test sequence from user
AUI41147_SEQUENCE = """MSLIEKPLTAIETREKPHAVCIPYPAQGHINPMMQLAKLLHHSGFHITFVHTEYNYDRLVKSQGSACVAGLPDFRFEAIPDGLPSTNGDVTQDIPLLSSSTSKTCLKPFKELLKRLQDKCKELPDDVPPLSCIVSDAAMSFTIDASEEFGVPIALLWTASACGFLGYTHYPYLIDRGVIPLKDESQLTNGYLDMSIDGIPCMEGIRLRDLPSFLRTTDLDDMMFSYILHEIKQVSRGSAIILNTFEALDHDVLDSLSKIYQNVILPVGPLHVSLNKIPKHYPLQSLSSNLWKDDTDCIPWLSSKASKSVIYVNFGSITTVSPKQIVEFAWGLANSKHPFLWIIRPDLVAGEASIIPQDFMDETKGRGLLAGWCDQELVLNHPSIGGFLTHCGWNSIIESISAGVPTVCWPFFAEQQTNCWFACKKWCIGMEMHTDVKRDEVDKLLRELMEGDKGEELKRKATNWKRLAEEAVSSTGLSTLNFRTLVNQVLLSKTKHIR"""

EXPECTED_NCBI_HITS = ['AUI41147', 'AUI41146', 'AUI41143', 'AUI41132', 'AUI41144']
UNEXPECTED_MMSEQS_HIT = 'AUI41142'


@pytest.fixture
def db_connection():
    """Connect to the master database."""
    db_path = Path("/mnt/data4/master.duckdb")
    if not db_path.exists():
        pytest.skip(f"Database not found at {db_path}")

    conn = duckdb.connect(str(db_path), read_only=True)
    yield conn
    conn.close()


@pytest.fixture
def repseq_fasta():
    """Path to the representative sequences FASTA file."""
    fasta_path = Path("/mnt/data4/repseq.faa")
    if not fasta_path.exists():
        pytest.skip(f"Representative sequences FASTA not found at {fasta_path}")
    return str(fasta_path)


@pytest.fixture
def repseq_v2_fasta():
    """Path to the v2 representative sequences FASTA file (cluster-mode 2 only)."""
    fasta_path = Path("/mnt/data4/planter_outputs/repseq_v2/output118/newRepSeqDB.fasta")
    if not fasta_path.exists():
        pytest.skip(f"Representative sequences v2 FASTA not found at {fasta_path}")
    return str(fasta_path)


@pytest.fixture
def repseq_v3_fasta():
    """Path to the v3 representative sequences FASTA file (cluster-mode 2 + cov-mode 1)."""
    fasta_path = Path("/mnt/data4/planter_outputs/repseq_v3/output130/newRepSeqDB.fasta")
    if not fasta_path.exists():
        pytest.skip(f"Representative sequences v3 FASTA not found at {fasta_path}")
    return str(fasta_path)


@pytest.fixture
def all_sequences_fasta(db_connection):
    """Path to ALL sequences FASTA file (generated from database)."""
    db_path = Path("/mnt/data4/master.duckdb")
    fasta_path = get_fasta_path(db_path, representatives_only=False)
    return str(fasta_path)


def test_sequence_length():
    """Verify AUI41147 sequence is 498 amino acids."""
    seq = AUI41147_SEQUENCE.replace('\n', '').strip()
    assert len(seq) == 498, f"Expected 498 aa, got {len(seq)} aa"


def test_show_config_paths(repseq_fasta, repseq_v2_fasta, db_connection):
    """Display the exact paths being used in tests vs what the app uses."""
    from app.config import config
    app_config = config['default']

    print(f"\n\n=== TEST CONFIGURATION ===")
    print(f"Test Database Path: /mnt/data4/master.duckdb")
    print(f"Test RepSeq BASELINE Path: {repseq_fasta}")
    print(f"Test RepSeq v2 Path: {repseq_v2_fasta}")

    print(f"\n=== APP CONFIGURATION ===")
    print(f"App Database Path: {app_config.DUCKDB_PATH}")
    print(f"App RepSeq Path: {app_config.REPSEQ_FASTA}")

    print(f"\n=== QUERY SEQUENCE ===")
    seq = AUI41147_SEQUENCE.replace('\n', '').strip()
    print(f"Accession: AUI41147")
    print(f"Length: {len(seq)} amino acids")
    print(f"Sequence: {seq}")


def test_mmseqs2_search_baseline(repseq_fasta):
    """Run MMSeqs2 search against BASELINE representative sequences (current production)."""
    sequence = AUI41147_SEQUENCE.replace('\n', '').strip()

    with tempfile.TemporaryDirectory() as temp_dir:
        input_file = os.path.join(temp_dir, "input.fasta")
        output_file = os.path.join(temp_dir, "output.tsv")
        tmp_dir = os.path.join(temp_dir, "tmp")

        # Write query sequence
        with open(input_file, 'w') as f:
            f.write(f">query\n{sequence}\n")

        # Run MMSeqs2 with default parameters from the app
        mmseqs_command = [
            "mmseqs", "easy-search", input_file, repseq_fasta, output_file, tmp_dir,
            "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tseq",
            "-s", "4.0",
            "-e", "0.001",
            "-c", "0.0",
            "--max-seqs", "300",
            "--alignment-mode", "3",
            "--mask", "1",
            "--min-seq-id", "0.0"
        ]

        process = subprocess.run(mmseqs_command, capture_output=True, text=True)

        assert process.returncode == 0, f"MMSeqs2 failed: {process.stderr}"
        assert os.path.exists(output_file), "MMSeqs2 did not produce output file"

        # Parse results
        df = pd.read_csv(output_file, sep='\t', names=[
            'query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen',
            'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'tseq'
        ])

        print(f"\n\nMMSeqs2 search returned {len(df)} hits")
        print(f"\nTop 10 hits:")
        print(df[['target', 'pident', 'evalue', 'bits']].head(10).to_string())

        # Store for next test
        return df


def test_mmseqs2_baseline_hits_in_database(repseq_fasta, db_connection):
    """Check what sequences MMSeqs2 found with BASELINE and examine them in the database."""
    sequence = AUI41147_SEQUENCE.replace('\n', '').strip()

    # Run MMSeqs2 search
    with tempfile.TemporaryDirectory() as temp_dir:
        input_file = os.path.join(temp_dir, "input.fasta")
        output_file = os.path.join(temp_dir, "output.tsv")
        tmp_dir = os.path.join(temp_dir, "tmp")

        with open(input_file, 'w') as f:
            f.write(f">query\n{sequence}\n")

        mmseqs_command = [
            "mmseqs", "easy-search", input_file, repseq_fasta, output_file, tmp_dir,
            "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tseq",
            "-s", "4.0",
            "-e", "0.001",
            "-c", "0.0",
            "--max-seqs", "300",
            "--alignment-mode", "3",
            "--mask", "1",
            "--min-seq-id", "0.0"
        ]

        process = subprocess.run(mmseqs_command, capture_output=True, text=True)
        assert process.returncode == 0, f"MMSeqs2 failed: {process.stderr}"

        df = pd.read_csv(output_file, sep='\t', names=[
            'query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen',
            'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'tseq'
        ])

    # Get top hits
    top_hits = df.nsmallest(10, 'evalue')
    seqhash_ids = top_hits['target'].tolist()

    print(f"\n\nTop 10 MMSeqs2 hits by E-value:")
    print(top_hits[['target', 'pident', 'evalue', 'bits']].to_string())

    # Query database for these sequences
    query = f"""
    SELECT
        s.seqhash_id,
        s.sample_id,
        s.length,
        s.is_representative,
        s.repseq_id,
        m.organism,
        a.description,
        a.preferred_name
    FROM sequences s
    LEFT JOIN sra_metadata m ON s.sample_id = m.sample_id
    LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
    WHERE s.seqhash_id IN ({','.join([f"'{x}'" for x in seqhash_ids])})
    ORDER BY s.length DESC
    """

    result = db_connection.execute(query).fetchdf()
    print(f"\n\nDatabase info for top MMSeqs2 hits:")
    print(result.to_string())

    # Check cluster information
    cluster_query = f"""
    SELECT
        cm.cluster_id,
        cm.seqhash_id,
        c.size as cluster_size,
        c.representative_seqhash_id
    FROM cluster_members cm
    JOIN clusters c ON cm.cluster_id = c.cluster_id
    WHERE cm.seqhash_id IN ({','.join([f"'{x}'" for x in seqhash_ids])})
    ORDER BY c.size DESC
    """

    cluster_result = db_connection.execute(cluster_query).fetchdf()
    print(f"\n\nCluster information for top hits:")
    print(cluster_result.to_string())


def test_mmseqs2_search_v2(repseq_v2_fasta, db_connection):
    """Run MMSeqs2 search against v2 representative sequences (after cluster fix)."""
    sequence = AUI41147_SEQUENCE.replace('\n', '').strip()

    # Run MMSeqs2 search
    with tempfile.TemporaryDirectory() as temp_dir:
        input_file = os.path.join(temp_dir, "input.fasta")
        output_file = os.path.join(temp_dir, "output.tsv")
        tmp_dir = os.path.join(temp_dir, "tmp")

        with open(input_file, 'w') as f:
            f.write(f">query\n{sequence}\n")

        mmseqs_command = [
            "mmseqs", "easy-search", input_file, repseq_v2_fasta, output_file, tmp_dir,
            "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tseq",
            "-s", "4.0",
            "-e", "0.001",
            "-c", "0.0",
            "--max-seqs", "300",
            "--alignment-mode", "3",
            "--mask", "1",
            "--min-seq-id", "0.0"
        ]

        process = subprocess.run(mmseqs_command, capture_output=True, text=True)
        assert process.returncode == 0, f"MMSeqs2 failed: {process.stderr}"

        df = pd.read_csv(output_file, sep='\t', names=[
            'query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen',
            'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'tseq'
        ])

    # Get top hits
    top_hits = df.nsmallest(10, 'evalue')
    seqhash_ids = top_hits['target'].tolist()

    print(f"\n\n=== V2 DATABASE RESULTS ===")
    print(f"Total hits: {len(df)}")
    print(f"\nTop 10 MMSeqs2 hits by E-value:")
    print(top_hits[['target', 'pident', 'evalue', 'bits']].to_string())

    # Query database for these sequences
    query = f"""
    SELECT
        s.seqhash_id,
        s.sample_id,
        s.length,
        s.is_representative,
        s.repseq_id,
        m.organism,
        a.description,
        a.preferred_name
    FROM sequences s
    LEFT JOIN sra_metadata m ON s.sample_id = m.sample_id
    LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
    WHERE s.seqhash_id IN ({','.join([f"'{x}'" for x in seqhash_ids])})
    ORDER BY s.length DESC
    """

    result = db_connection.execute(query).fetchdf()
    print(f"\n\nDatabase info for top MMSeqs2 hits (v2):")
    print(result.to_string())

    # Check cluster information
    cluster_query = f"""
    SELECT
        cm.cluster_id,
        cm.seqhash_id,
        c.size as cluster_size,
        c.representative_seqhash_id
    FROM cluster_members cm
    JOIN clusters c ON cm.cluster_id = c.cluster_id
    WHERE cm.seqhash_id IN ({','.join([f"'{x}'" for x in seqhash_ids])})
    ORDER BY c.size DESC
    """

    cluster_result = db_connection.execute(cluster_query).fetchdf()
    print(f"\n\nCluster information for top hits (v2):")
    print(cluster_result.to_string())

    # Check specifically for Rhodiola hits
    rhodiola_hits = result[result['organism'].str.contains('Rhodiola', case=False, na=False)]
    print(f"\n\n=== RHODIOLA HITS IN V2 ===")
    print(f"Found {len(rhodiola_hits)} Rhodiola hits in top 10")
    if not rhodiola_hits.empty:
        print(rhodiola_hits[['seqhash_id', 'length', 'organism', 'preferred_name']].to_string())


def test_compare_baseline_vs_v2(repseq_fasta, repseq_v2_fasta, db_connection):
    """Compare BASELINE vs V2 clustering results side-by-side."""
    sequence = AUI41147_SEQUENCE.replace('\n', '').strip()

    print(f"\n\n{'='*80}")
    print(f"COMPARISON: BASELINE vs V2 CLUSTERING")
    print(f"Query: AUI41147 (498 aa Rhodiola rosea UGT)")
    print(f"{'='*80}")

    results = {}

    for label, fasta_path in [("BASELINE", repseq_fasta), ("V2", repseq_v2_fasta)]:
        # Run MMSeqs2 search
        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = os.path.join(temp_dir, "input.fasta")
            output_file = os.path.join(temp_dir, "output.tsv")
            tmp_dir = os.path.join(temp_dir, "tmp")

            with open(input_file, 'w') as f:
                f.write(f">query\n{sequence}\n")

            mmseqs_command = [
                "mmseqs", "easy-search", input_file, fasta_path, output_file, tmp_dir,
                "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tseq",
                "-s", "4.0", "-e", "0.001", "-c", "0.0", "--max-seqs", "300",
                "--alignment-mode", "3", "--mask", "1", "--min-seq-id", "0.0"
            ]

            process = subprocess.run(mmseqs_command, capture_output=True, text=True)
            assert process.returncode == 0, f"MMSeqs2 failed for {label}: {process.stderr}"

            df = pd.read_csv(output_file, sep='\t', names=[
                'query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen',
                'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'tseq'
            ])

        # Get top 10 hits
        top_hits = df.nsmallest(10, 'evalue')
        seqhash_ids = top_hits['target'].tolist()

        # Query database
        query = f"""
        SELECT
            s.seqhash_id,
            s.length,
            m.organism,
            a.preferred_name
        FROM sequences s
        LEFT JOIN sra_metadata m ON s.sample_id = m.sample_id
        LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
        WHERE s.seqhash_id IN ({','.join([f"'{x}'" for x in seqhash_ids])})
        """

        db_result = db_connection.execute(query).fetchdf()

        # Merge with MMSeqs2 results
        merged = top_hits.merge(db_result, left_on='target', right_on='seqhash_id', how='left')
        merged = merged[['target', 'pident', 'evalue', 'bits', 'length', 'organism', 'preferred_name']]

        results[label] = merged

    # Print side-by-side comparison
    print(f"\n{'='*80}")
    print("BASELINE (current production)")
    print(f"{'='*80}")
    print(results['BASELINE'].to_string(index=True))

    print(f"\n{'='*80}")
    print("V2 (new clustering with full-length preference)")
    print(f"{'='*80}")
    print(results['V2'].to_string(index=True))

    # Analyze Rhodiola hits
    print(f"\n{'='*80}")
    print("RHODIOLA ROSEA COMPARISON")
    print(f"{'='*80}")

    baseline_rhodiola = results['BASELINE'][
        results['BASELINE']['organism'].str.contains('Rhodiola', case=False, na=False)
    ]
    v2_rhodiola = results['V2'][
        results['V2']['organism'].str.contains('Rhodiola', case=False, na=False)
    ]

    print(f"\nBASELINE Rhodiola hits in top 10:")
    if not baseline_rhodiola.empty:
        for idx, row in baseline_rhodiola.iterrows():
            print(f"  Rank #{idx+1}: {row['length']} aa, E-value: {row['evalue']:.2e}, {row['preferred_name']}")
    else:
        print("  None")

    print(f"\nV2 Rhodiola hits in top 10:")
    if not v2_rhodiola.empty:
        for idx, row in v2_rhodiola.iterrows():
            print(f"  Rank #{idx+1}: {row['length']} aa, E-value: {row['evalue']:.2e}, {row['preferred_name']}")
    else:
        print("  None")

    # Quality metrics
    print(f"\n{'='*80}")
    print("QUALITY METRICS")
    print(f"{'='*80}")

    for label, df in results.items():
        avg_length = df['length'].mean()
        median_length = df['length'].median()
        min_length = df['length'].min()
        max_length = df['length'].max()

        print(f"\n{label}:")
        print(f"  Average sequence length: {avg_length:.1f} aa")
        print(f"  Median sequence length: {median_length:.1f} aa")
        print(f"  Length range: {min_length}-{max_length} aa")
        print(f"  Sequences < 250 aa: {len(df[df['length'] < 250])}/10")
        print(f"  Sequences >= 400 aa: {len(df[df['length'] >= 400])}/10")

    # Verdict
    print(f"\n{'='*80}")
    print("ASSESSMENT")
    print(f"{'='*80}")

    baseline_avg = results['BASELINE']['length'].mean()
    v2_avg = results['V2']['length'].mean()

    if v2_avg > baseline_avg:
        improvement = ((v2_avg - baseline_avg) / baseline_avg) * 100
        print(f"✓ V2 shows IMPROVEMENT: {improvement:.1f}% longer sequences on average")
    else:
        decline = ((baseline_avg - v2_avg) / baseline_avg) * 100
        print(f"✗ V2 shows DECLINE: {decline:.1f}% shorter sequences on average")

    baseline_short = len(results['BASELINE'][results['BASELINE']['length'] < 250])
    v2_short = len(results['V2'][results['V2']['length'] < 250])

    if v2_short < baseline_short:
        print(f"✓ V2 reduces truncated sequences: {baseline_short} → {v2_short} (sequences < 250 aa)")
    elif v2_short > baseline_short:
        print(f"✗ V2 increases truncated sequences: {baseline_short} → {v2_short} (sequences < 250 aa)")
    else:
        print(f"= Same number of truncated sequences in both: {v2_short}")

    # Rhodiola specific
    if not baseline_rhodiola.empty and not v2_rhodiola.empty:
        baseline_rhodiola_rank = baseline_rhodiola.index[0] + 1
        v2_rhodiola_rank = v2_rhodiola.index[0] + 1
        baseline_rhodiola_length = baseline_rhodiola.iloc[0]['length']
        v2_rhodiola_length = v2_rhodiola.iloc[0]['length']

        print(f"\nRhodiola rosea hit:")
        print(f"  BASELINE: Rank #{baseline_rhodiola_rank}, {baseline_rhodiola_length} aa")
        print(f"  V2: Rank #{v2_rhodiola_rank}, {v2_rhodiola_length} aa")

        if v2_rhodiola_rank > baseline_rhodiola_rank:
            print(f"  ✓ V2 pushes truncated Rhodiola DOWN in rankings (better)")
        else:
            print(f"  ✗ V2 does not improve Rhodiola ranking")


def test_detailed_cluster_analysis(db_connection):
    """Detailed analysis of specific clusters to understand why V2 has shorter sequences."""

    print(f"\n\n{'='*80}")
    print("DETAILED CLUSTER ANALYSIS")
    print(f"{'='*80}")

    # Let's look at specific sequences from the comparison
    # Pick a few representative cases to understand what's happening

    sequences_to_analyze = [
        # BASELINE top hit (451 aa) - what happened to this in V2?
        'v1_DLS_fdcf2e4a73bbfc0989fe9358272fce47d76318ec32d9581d3799236695368ae1.p1',
        # V2 top hit (213 aa) - why was this chosen?
        'v1_DLS_8049c2d6159cd0017081c4c9ea0044ddcdc20d077ba71190ea604bf019c3b7c9.p1',
        # BASELINE Rhodiola (241 aa)
        'v1_DLS_38e612b34e68d0fc6d2620a80d0b20bbee0e86c1aa057e221db11bb3caafb8cb.p2',
    ]

    for seq_id in sequences_to_analyze:
        print(f"\n{'-'*80}")
        print(f"Analyzing: {seq_id}")
        print(f"{'-'*80}")

        # Get basic info about this sequence
        seq_query = f"""
        SELECT
            s.seqhash_id,
            s.length,
            s.is_representative,
            s.repseq_id,
            m.organism,
            a.preferred_name
        FROM sequences s
        LEFT JOIN sra_metadata m ON s.sample_id = m.sample_id
        LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
        WHERE s.seqhash_id = '{seq_id}'
        """

        seq_info = db_connection.execute(seq_query).fetchdf()
        if seq_info.empty:
            print(f"  NOT FOUND in database")
            continue

        row = seq_info.iloc[0]
        print(f"  Length: {row['length']} aa")
        print(f"  Organism: {row['organism']}")
        print(f"  Is Representative: {row['is_representative']}")
        print(f"  Representative ID: {row['repseq_id']}")
        print(f"  Annotation: {row['preferred_name']}")

        # Get cluster information
        cluster_query = f"""
        SELECT
            cm.cluster_id,
            c.size as cluster_size,
            c.representative_seqhash_id
        FROM cluster_members cm
        JOIN clusters c ON cm.cluster_id = c.cluster_id
        WHERE cm.seqhash_id = '{seq_id}'
        """

        cluster_info = db_connection.execute(cluster_query).fetchdf()
        if not cluster_info.empty:
            cluster_row = cluster_info.iloc[0]
            print(f"\n  Cluster ID: {cluster_row['cluster_id']}")
            print(f"  Cluster Size: {cluster_row['cluster_size']} members")
            print(f"  Cluster Representative: {cluster_row['representative_seqhash_id']}")

            # Get all members of this cluster to see length distribution
            members_query = f"""
            SELECT
                s.seqhash_id,
                s.length,
                s.is_representative,
                m.organism
            FROM cluster_members cm
            JOIN sequences s ON cm.seqhash_id = s.seqhash_id
            LEFT JOIN sra_metadata m ON s.sample_id = m.sample_id
            WHERE cm.cluster_id = '{cluster_row['cluster_id']}'
            ORDER BY s.length DESC
            LIMIT 20
            """

            members = db_connection.execute(members_query).fetchdf()
            print(f"\n  Top 20 cluster members by length:")
            print(f"    Total members: {len(members)}")
            print(f"    Length range: {members['length'].min()}-{members['length'].max()} aa")
            print(f"    Average length: {members['length'].mean():.1f} aa")
            print(f"    Median length: {members['length'].median():.1f} aa")

            # Show top 5 longest
            print(f"\n  Top 5 longest members:")
            for idx, member in members.head(5).iterrows():
                rep_marker = "★ REP" if member['is_representative'] else ""
                print(f"    {member['length']:4d} aa - {member['organism'][:40]:40s} {rep_marker}")

            # Check if there are longer sequences that should have been representative
            longer_than_rep = members[members['length'] > row['length']]
            if not longer_than_rep.empty and row['is_representative']:
                print(f"\n  ⚠️  WARNING: This is representative but {len(longer_than_rep)} longer sequences exist in cluster!")
                print(f"      Longest alternative: {longer_than_rep.iloc[0]['length']} aa")


def test_check_clustering_algorithm(db_connection):
    """Check how representatives were chosen across all clusters."""

    print(f"\n\n{'='*80}")
    print("CLUSTERING ALGORITHM ANALYSIS")
    print(f"{'='*80}")

    # For each cluster, check if the representative is the longest
    query = """
    WITH cluster_stats AS (
        SELECT
            cm.cluster_id,
            c.representative_seqhash_id,
            MAX(s.length) as max_length,
            MIN(s.length) as min_length,
            AVG(s.length) as avg_length,
            COUNT(*) as member_count
        FROM cluster_members cm
        JOIN sequences s ON cm.seqhash_id = s.seqhash_id
        JOIN clusters c ON cm.cluster_id = c.cluster_id
        GROUP BY cm.cluster_id, c.representative_seqhash_id
    ),
    rep_lengths AS (
        SELECT
            cs.cluster_id,
            cs.representative_seqhash_id,
            s.length as rep_length,
            cs.max_length,
            cs.min_length,
            cs.avg_length,
            cs.member_count
        FROM cluster_stats cs
        JOIN sequences s ON cs.representative_seqhash_id = s.seqhash_id
    )
    SELECT
        COUNT(*) as total_clusters,
        SUM(CASE WHEN rep_length = max_length THEN 1 ELSE 0 END) as reps_are_longest,
        SUM(CASE WHEN rep_length < max_length THEN 1 ELSE 0 END) as reps_not_longest,
        SUM(CASE WHEN rep_length = min_length THEN 1 ELSE 0 END) as reps_are_shortest,
        AVG(rep_length) as avg_rep_length,
        AVG(max_length) as avg_max_length,
        AVG(avg_length) as avg_cluster_avg_length
    FROM rep_lengths
    WHERE member_count > 1
    """

    stats = db_connection.execute(query).fetchdf()

    print(f"\nOverall clustering statistics (multi-member clusters):")
    row = stats.iloc[0]
    print(f"  Total clusters: {row['total_clusters']}")
    print(f"  Representatives that ARE the longest: {row['reps_are_longest']} ({row['reps_are_longest']/row['total_clusters']*100:.1f}%)")
    print(f"  Representatives that are NOT longest: {row['reps_not_longest']} ({row['reps_not_longest']/row['total_clusters']*100:.1f}%)")
    print(f"  Representatives that are shortest: {row['reps_are_shortest']} ({row['reps_are_shortest']/row['total_clusters']*100:.1f}%)")
    print(f"\n  Average representative length: {row['avg_rep_length']:.1f} aa")
    print(f"  Average maximum length in clusters: {row['avg_max_length']:.1f} aa")
    print(f"  Average cluster mean length: {row['avg_cluster_avg_length']:.1f} aa")

    # Find worst offenders - clusters where representative is much shorter than longest
    print(f"\n{'-'*80}")
    print("Worst offenders - clusters with short reps when long sequences available:")
    print(f"{'-'*80}")

    worst_query = """
    WITH cluster_stats AS (
        SELECT
            cm.cluster_id,
            c.representative_seqhash_id,
            MAX(s.length) as max_length,
            COUNT(*) as member_count
        FROM cluster_members cm
        JOIN sequences s ON cm.seqhash_id = s.seqhash_id
        JOIN clusters c ON cm.cluster_id = c.cluster_id
        GROUP BY cm.cluster_id, c.representative_seqhash_id
    ),
    rep_lengths AS (
        SELECT
            cs.cluster_id,
            cs.representative_seqhash_id,
            s.length as rep_length,
            cs.max_length,
            cs.member_count,
            (cs.max_length - s.length) as length_deficit
        FROM cluster_stats cs
        JOIN sequences s ON cs.representative_seqhash_id = s.seqhash_id
        WHERE cs.member_count > 1
    )
    SELECT
        cluster_id,
        representative_seqhash_id,
        rep_length,
        max_length,
        length_deficit,
        member_count
    FROM rep_lengths
    WHERE length_deficit > 200
    ORDER BY length_deficit DESC
    LIMIT 10
    """

    worst = db_connection.execute(worst_query).fetchdf()

    if not worst.empty:
        print(f"\nTop 10 clusters with biggest length deficits:")
        print(worst.to_string(index=False))

        print(f"\n⚠️  PROBLEM IDENTIFIED:")
        print(f"    {len(worst)} clusters have representatives >200 aa shorter than longest member")
        print(f"    This suggests the clustering algorithm is NOT preferring full-length sequences")
    else:
        print(f"\n✓ No major length deficits found")


def test_search_all_sequences(all_sequences_fasta, db_connection):
    """Search against ALL sequences (not just representatives) to see if we find better hits."""
    sequence = AUI41147_SEQUENCE.replace('\n', '').strip()

    print(f"\n\n{'='*80}")
    print(f"SEARCH ALL SEQUENCES (NOT JUST REPRESENTATIVES)")
    print(f"Query: AUI41147 (498 aa Rhodiola rosea UGT)")
    print(f"{'='*80}")

    # Run MMSeqs2 search
    with tempfile.TemporaryDirectory() as temp_dir:
        input_file = os.path.join(temp_dir, "input.fasta")
        output_file = os.path.join(temp_dir, "output.tsv")
        tmp_dir = os.path.join(temp_dir, "tmp")

        with open(input_file, 'w') as f:
            f.write(f">query\n{sequence}\n")

        mmseqs_command = [
            "mmseqs", "easy-search", input_file, all_sequences_fasta, output_file, tmp_dir,
            "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tseq",
            "-s", "4.0", "-e", "0.001", "-c", "0.0", "--max-seqs", "300",
            "--alignment-mode", "3", "--mask", "1", "--min-seq-id", "0.0"
        ]

        process = subprocess.run(mmseqs_command, capture_output=True, text=True)
        assert process.returncode == 0, f"MMSeqs2 failed: {process.stderr}"

        df = pd.read_csv(output_file, sep='\t', names=[
            'query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen',
            'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'tseq'
        ])

    # Get top 10 hits
    top_hits = df.nsmallest(10, 'evalue')
    seqhash_ids = top_hits['target'].tolist()

    print(f"\nTotal hits: {len(df)}")
    print(f"\nTop 10 hits by E-value:")
    print(top_hits[['target', 'pident', 'evalue', 'bits']].to_string())

    # Query database
    query = f"""
    SELECT
        s.seqhash_id,
        s.length,
        s.is_representative,
        m.organism,
        a.preferred_name
    FROM sequences s
    LEFT JOIN sra_metadata m ON s.sample_id = m.sample_id
    LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
    WHERE s.seqhash_id IN ({','.join([f"'{x}'" for x in seqhash_ids])})
    """

    db_result = db_connection.execute(query).fetchdf()

    # Merge with MMSeqs2 results
    merged = top_hits.merge(db_result, left_on='target', right_on='seqhash_id', how='left')
    merged = merged[['target', 'pident', 'evalue', 'bits', 'length', 'is_representative', 'organism', 'preferred_name']]

    print(f"\nDatabase info for top hits:")
    print(merged.to_string(index=True))

    # Check for Rhodiola
    rhodiola_hits = merged[merged['organism'].str.contains('Rhodiola', case=False, na=False)]
    print(f"\n{'='*80}")
    print("RHODIOLA ROSEA HITS")
    print(f"{'='*80}")

    if not rhodiola_hits.empty:
        print(f"\nFound {len(rhodiola_hits)} Rhodiola hits in top 10:")
        for idx, row in rhodiola_hits.iterrows():
            rep_status = "REP" if row['is_representative'] else "NON-REP"
            print(f"  Rank #{idx+1}: {row['length']} aa [{rep_status}], E-value: {row['evalue']:.2e}, {row['preferred_name']}")
    else:
        print(f"\nNo Rhodiola hits in top 10")

    # Quality metrics
    print(f"\n{'='*80}")
    print("QUALITY METRICS (searching all sequences)")
    print(f"{'='*80}")

    avg_length = merged['length'].mean()
    median_length = merged['length'].median()
    representatives = len(merged[merged['is_representative'] == True])

    print(f"  Average sequence length: {avg_length:.1f} aa")
    print(f"  Median sequence length: {median_length:.1f} aa")
    print(f"  Sequences < 250 aa: {len(merged[merged['length'] < 250])}/10")
    print(f"  Sequences >= 400 aa: {len(merged[merged['length'] >= 400])}/10")
    print(f"  How many are representatives: {representatives}/10")


def test_compare_all_versions(repseq_fasta, repseq_v2_fasta, repseq_v3_fasta, db_connection):
    """Compare BASELINE vs V2 vs V3 clustering results."""
    sequence = AUI41147_SEQUENCE.replace('\n', '').strip()

    print(f"\n\n{'='*80}")
    print(f"COMPREHENSIVE COMPARISON: BASELINE vs V2 vs V3")
    print(f"Query: AUI41147 (498 aa Rhodiola rosea UGT)")
    print(f"{'='*80}")

    results = {}

    for label, fasta_path in [
        ("BASELINE (mode 0)", repseq_fasta),
        ("V2 (mode 2 only)", repseq_v2_fasta),
        ("V3 (mode 2 + cov 1)", repseq_v3_fasta)
    ]:
        # Run MMSeqs2 search
        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = os.path.join(temp_dir, "input.fasta")
            output_file = os.path.join(temp_dir, "output.tsv")
            tmp_dir = os.path.join(temp_dir, "tmp")

            with open(input_file, 'w') as f:
                f.write(f">query\n{sequence}\n")

            mmseqs_command = [
                "mmseqs", "easy-search", input_file, fasta_path, output_file, tmp_dir,
                "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tseq",
                "-s", "4.0", "-e", "0.001", "-c", "0.0", "--max-seqs", "300",
                "--alignment-mode", "3", "--mask", "1", "--min-seq-id", "0.0"
            ]

            process = subprocess.run(mmseqs_command, capture_output=True, text=True)
            if process.returncode != 0:
                print(f"\n{label} FAILED: {process.stderr}")
                continue

            df = pd.read_csv(output_file, sep='\t', names=[
                'query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen',
                'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'tseq'
            ])

        # Get top 10 hits
        top_hits = df.nsmallest(10, 'evalue')
        seqhash_ids = top_hits['target'].tolist()

        # Query database
        query = f"""
        SELECT
            s.seqhash_id,
            s.length,
            m.organism,
            a.preferred_name
        FROM sequences s
        LEFT JOIN sra_metadata m ON s.sample_id = m.sample_id
        LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
        WHERE s.seqhash_id IN ({','.join([f"'{x}'" for x in seqhash_ids])})
        """

        db_result = db_connection.execute(query).fetchdf()

        # Merge with MMSeqs2 results
        merged = top_hits.merge(db_result, left_on='target', right_on='seqhash_id', how='left')
        merged = merged[['target', 'pident', 'evalue', 'bits', 'length', 'organism', 'preferred_name']]

        results[label] = merged

    # Print comparison
    for label in ["BASELINE (mode 0)", "V2 (mode 2 only)", "V3 (mode 2 + cov 1)"]:
        if label not in results:
            continue

        print(f"\n{'='*80}")
        print(f"{label}")
        print(f"{'='*80}")
        print(results[label].to_string(index=True))

        # Rhodiola hits
        rhodiola = results[label][results[label]['organism'].str.contains('Rhodiola', case=False, na=False)]
        if not rhodiola.empty:
            print(f"\nRhodiola hits:")
            for idx, row in rhodiola.iterrows():
                print(f"  Rank #{idx+1}: {row['length']} aa, E-value: {row['evalue']:.2e}, {row['preferred_name']}")

    # Summary comparison
    print(f"\n{'='*80}")
    print("SUMMARY COMPARISON")
    print(f"{'='*80}")

    summary_data = []
    for label, df in results.items():
        avg_length = df['length'].mean()
        median_length = df['length'].median()
        short_count = len(df[df['length'] < 250])
        long_count = len(df[df['length'] >= 400])

        rhodiola = df[df['organism'].str.contains('Rhodiola', case=False, na=False)]
        rhodiola_rank = rhodiola.index[0] + 1 if not rhodiola.empty else None
        rhodiola_length = rhodiola.iloc[0]['length'] if not rhodiola.empty else None

        summary_data.append({
            'Version': label,
            'Avg Length': f"{avg_length:.1f}",
            'Median Length': f"{median_length:.1f}",
            '< 250 aa': short_count,
            '>= 400 aa': long_count,
            'Rhodiola Rank': rhodiola_rank,
            'Rhodiola Length': rhodiola_length
        })

    summary_df = pd.DataFrame(summary_data)
    print(f"\n{summary_df.to_string(index=False)}")

    # Final verdict
    print(f"\n{'='*80}")
    print("VERDICT")
    print(f"{'='*80}")

    if "V3 (mode 2 + cov 1)" in results:
        v3 = results["V3 (mode 2 + cov 1)"]
        baseline = results["BASELINE (mode 0)"]

        v3_avg = v3['length'].mean()
        baseline_avg = baseline['length'].mean()

        if v3_avg > baseline_avg:
            improvement = ((v3_avg - baseline_avg) / baseline_avg) * 100
            print(f"✅ V3 shows IMPROVEMENT: {improvement:.1f}% longer sequences on average")
        else:
            print(f"❌ V3 does not improve average length")

        v3_short = len(v3[v3['length'] < 250])
        baseline_short = len(baseline[baseline['length'] < 250])

        if v3_short < baseline_short:
            print(f"✅ V3 reduces truncated sequences: {baseline_short} → {v3_short}")
        else:
            print(f"❌ V3 does not reduce truncated sequences")

        v3_long = len(v3[v3['length'] >= 400])
        baseline_long = len(baseline[baseline['length'] >= 400])

        if v3_long > baseline_long:
            print(f"✅ V3 increases full-length sequences: {baseline_long} → {v3_long}")
        else:
            print(f"⚠️  V3 does not increase full-length sequences")
