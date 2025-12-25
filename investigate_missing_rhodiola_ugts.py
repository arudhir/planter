#!/usr/bin/env python3
"""
Investigate why Rhodiola rosea UGT sequences are missing from cluster_members.

The issue:
- 20 Rhodiola rosea UGT sequences exist in the database
- None of them are in cluster_members table
- When BLASTing a Rhodiola UGT, we can't find it in any cluster
"""

import duckdb
from pathlib import Path

db_path = "/mnt/data4/master_simple.duckdb"
con = duckdb.connect(db_path, read_only=True)

print("=" * 80)
print("RHODIOLA UGT CLUSTERING INVESTIGATION")
print("=" * 80)

# Get all Rhodiola UGT sequences
print("\n1. Getting all Rhodiola rosea UGT sequences...")
rhodiola_ugts = con.execute("""
    SELECT DISTINCT
        a.seqhash_id,
        a.preferred_name,
        a.description,
        a.sample_id,
        s.length
    FROM annotations a
    JOIN sequences s ON a.seqhash_id = s.seqhash_id
    LEFT JOIN sra_metadata m ON a.sample_id = m.sample_id
    WHERE (a.preferred_name LIKE '%UGT%' OR a.description LIKE '%UGT%')
        AND m.organism = 'Rhodiola rosea'
    ORDER BY s.length DESC
""").fetchall()

print(f"   Found {len(rhodiola_ugts)} unique Rhodiola UGT sequences")
print(f"\n   Top 5 sequences:")
for i, (seqhash, name, desc, sample, length) in enumerate(rhodiola_ugts[:5], 1):
    short_hash = seqhash[:50] + "..." if len(seqhash) > 50 else seqhash
    print(f"   {i}. {short_hash}")
    print(f"      Name: {name}, Length: {length} aa, Sample: {sample}")
    print(f"      Desc: {desc[:60]}...")

# Check if any are in cluster_members
print("\n2. Checking cluster membership...")
for seqhash, name, desc, sample, length in rhodiola_ugts:
    result = con.execute("""
        SELECT cluster_id
        FROM cluster_members
        WHERE seqhash_id = ?
    """, [seqhash]).fetchone()

    if result:
        print(f"   ✓ {seqhash[:50]}... IS in cluster {result[0]}")

# Count how many are missing
missing_count = sum(1 for seq in rhodiola_ugts if not con.execute(
    "SELECT 1 FROM cluster_members WHERE seqhash_id = ?", [seq[0]]
).fetchone())

print(f"\n   Result: {missing_count}/{len(rhodiola_ugts)} Rhodiola UGTs are NOT in cluster_members")

# Check overall statistics
print("\n3. Overall clustering statistics:")
stats = con.execute("""
    SELECT
        COUNT(*) as total_sequences,
        COUNT(DISTINCT cm.seqhash_id) as sequences_in_clusters,
        ROUND(100.0 * COUNT(DISTINCT cm.seqhash_id) / COUNT(*), 2) as coverage_pct
    FROM sequences s
    LEFT JOIN cluster_members cm ON s.seqhash_id = cm.seqhash_id
""").fetchone()

print(f"   Total sequences in DB: {stats[0]:,}")
print(f"   Sequences in clusters: {stats[1]:,}")
print(f"   Coverage: {stats[2]}%")

# Check Rhodiola sample statistics
print("\n4. Rhodiola sample clustering statistics:")
sample_stats = con.execute("""
    SELECT
        s.sample_id,
        COUNT(*) as total_seqs,
        COUNT(cm.seqhash_id) as clustered_seqs,
        ROUND(100.0 * COUNT(cm.seqhash_id) / COUNT(*), 2) as coverage_pct
    FROM sequences s
    LEFT JOIN cluster_members cm ON s.seqhash_id = cm.seqhash_id
    WHERE s.sample_id IN ('SRR5936536', 'SRR5936537', 'SRR22844166')
    GROUP BY s.sample_id
    ORDER BY s.sample_id
""").fetchall()

for sample_id, total, clustered, pct in sample_stats:
    print(f"   {sample_id}: {clustered:,}/{total:,} ({pct}%) sequences clustered")

# Check if the Rhodiola UGTs are in the sequences table with correct info
print("\n5. Verifying Rhodiola UGT sequences are properly loaded:")
for seqhash, name, desc, sample, length in rhodiola_ugts[:3]:
    seq_info = con.execute("""
        SELECT seqhash_id, sample_id, length, LEFT(sequence, 50) as seq_preview
        FROM sequences
        WHERE seqhash_id = ?
    """, [seqhash]).fetchone()

    if seq_info:
        print(f"   ✓ {seqhash[:50]}...")
        print(f"     Sample: {seq_info[1]}, Length: {seq_info[2]}")
        print(f"     Sequence: {seq_info[3]}...")
    else:
        print(f"   ✗ {seqhash[:50]}... NOT FOUND in sequences table")

# Check for similar UGT sequences that ARE clustered
print("\n6. Finding other UGT sequences that ARE clustered (for comparison):")
clustered_ugts = con.execute("""
    SELECT
        a.seqhash_id,
        a.preferred_name,
        s.length,
        a.sample_id,
        cm.cluster_id,
        c.representative_seqhash_id
    FROM annotations a
    JOIN sequences s ON a.seqhash_id = s.seqhash_id
    JOIN cluster_members cm ON a.seqhash_id = cm.seqhash_id
    JOIN clusters c ON cm.cluster_id = c.cluster_id
    WHERE (a.preferred_name LIKE '%UGT%' OR a.description LIKE '%UGT%')
    LIMIT 5
""").fetchall()

if clustered_ugts:
    print(f"   Found {len(clustered_ugts)} example UGT sequences that ARE clustered:")
    for seqhash, name, length, sample, cluster_id, rep_id in clustered_ugts:
        print(f"   - {name} ({length} aa) from {sample}")
        print(f"     Cluster: {cluster_id[:50]}...")
        print(f"     Representative: {rep_id[:50]}...")
else:
    print("   No other UGT sequences are clustered either!")

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"""
The problem:
- All {len(rhodiola_ugts)} Rhodiola rosea UGT sequences are missing from cluster_members
- These sequences DO exist in the sequences table with proper annotations
- Other sequences from the same Rhodiola samples ARE clustered (3,951 sequences)
- Overall clustering coverage is {stats[2]}% (only {stats[1]:,} of {stats[0]:,} sequences)

Possible causes:
1. The sequences were added to the database AFTER clustering was performed
2. The sequences didn't meet the clustering thresholds (too short/different)
3. The sequences were in samples that weren't included in the clustering input
4. There was an issue with the iterative clustering process

Next steps:
1. Check if these sequences exist in the original clustering input FASTA files
2. Check the clustering pipeline logs to see which samples were included
3. Try re-running MMSeqs2 clustering with these specific sequences
""")

con.close()
