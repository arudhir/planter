# Rhodiola UGT Missing from Clusters - Investigation Summary

**Date:** 2025-10-13
**Database:** `/mnt/data4/master_simple.duckdb`
**Cluster File:** `/mnt/data4/planter_outputs/repseq_v3/output130/newClusterDB.tsv`

## Problem Statement

When BLASTing a Rhodiola rosea UGT sequence (e.g., AUI41147, 485 aa) against the representative sequences:
- Expected: Should find cluster representatives with Rhodiola sequences in their cluster_members
- Actual: Cannot find the Rhodiola UGT sequences in ANY cluster membership

## Key Findings

### 1. **All Rhodiola UGT sequences are missing from cluster_members**

- **20 unique Rhodiola rosea UGT sequences** exist in the database
- **0 of them** are in the `cluster_members` table
- These sequences ARE properly loaded in the `sequences` table with annotations
- Length range: 240-616 aa (typical UGT size ~485 aa)

Example sequences:
```
v1_DLS_ed0a9b8e87e9ff23685d90c4a18677202888a4ae9ef5a32411e2ae16f40d70f0.p1 (499 aa, UGT85A1)
v1_DLS_5908d838437e1f4cddc2ce7b0c536fb73ead6f6466eb71b9ed0fab352a0413c3.p1 (488 aa, UGT72E1)
v1_DLS_7d3cf135d9c7105733e7b7305ef9308ad5320dfac19544e98aae1349bf2f4e85.p1 (616 aa, UGT80A2-like)
```

### 2. **Rhodiola samples ARE in the database and partially clustered**

Three Rhodiola rosea samples in the database:
- SRR5936536: 672/23,312 (2.88%) sequences clustered
- SRR5936537: 1,260/57,390 (2.20%) sequences clustered
- SRR22844166: 2,019/27,417 (7.36%) sequences clustered

**Total: 3,951 sequences from Rhodiola ARE clustered, but the UGTs are NOT**

### 3. **Very low overall clustering coverage**

- Total sequences in database: 1,970,849
- Sequences with cluster assignments: 169,863
- **Coverage: Only 8.62%**

This low coverage suggests that:
- Most sequences were added to the database AFTER clustering
- OR the iterative clustering process didn't include all sequences

### 4. **Other UGT sequences ARE clustered (for comparison)**

Found examples of clustered UGTs from other samples:
- UGT85A1 (145 aa) from SRR19071499 - IS clustered
- UGT2A1 (147 aa) from SRR18735292 - IS clustered
- UGT fragments (101-275 aa) - ARE clustered

**Key observation:** Shorter UGT fragments ARE clustered, but full-length Rhodiola UGTs (488-616 aa) are NOT.

### 5. **Sequences were not in clustering input FASTA**

Checked `/mnt/data4/planter_outputs/repseq_v3/output130/newRepSeqDB.fasta`:
- Rhodiola UGT sequences: **NOT FOUND** (0 occurrences)
- This confirms they were never included in the MMSeqs2 clustering process

## Root Cause Analysis

The Rhodiola UGT sequences were **added to the database AFTER the clustering process was complete**.

Evidence:
1. Sequences exist in `sequences` table ✓
2. Sequences have proper annotations ✓
3. Sequences are NOT in `cluster_members` table ✗
4. Sequences are NOT in the representative FASTA files ✗
5. Sequences are NOT in the cluster TSV file ✗
6. Only 8.62% of database sequences are clustered (suggests clustering happened early, before most samples were added)

## Expected vs Actual Behavior

### Expected (BLAST-like behavior)
When searching for Rhodiola UGT AUI41147:
- Should find ~5 similar sequences (like NCBI BLAST: AUI41147, AUI41146, AUI41143, AUI41132, AUI41144)
- Should find a cluster representative
- That representative's cluster_members should include the Rhodiola UGT sequences

### Actual
- The Rhodiola UGT sequences don't exist in ANY cluster
- When you BLAST/search, you cannot trace back to these sequences
- The cluster membership is incomplete

## Impact

This causes a **data discoverability problem**:
1. Users searching for Rhodiola UGTs won't find them through cluster representatives
2. The web app query "show me all sequences similar to X" won't return Rhodiola UGTs if they're not clustered
3. Cluster-based analyses will miss these sequences entirely

## Recommended Solutions

### Option 1: Re-run MMSeqs2 clustering on the complete database ⭐ (BEST)

```bash
# Re-cluster using ALL sequences currently in the database
python scripts/run_iterative_clustering.py \
  --input-db /mnt/data4/master_simple.duckdb \
  --output-dir /mnt/data4/planter_outputs/repseq_v4 \
  --threads 16
```

**Pros:**
- Will include ALL sequences, including Rhodiola UGTs
- Cluster memberships will be complete and accurate
- Solves the problem comprehensively

**Cons:**
- Computationally expensive (may take hours/days)
- Need to reload cluster data into database

### Option 2: Run incremental clustering update

Use MMSeqs2's `clusterupdate` to add new sequences to existing clusters:

```bash
# Extract unclustered sequences
python scripts/extract_unclustered_sequences.py \
  -d /mnt/data4/master_simple.duckdb \
  -o /tmp/unclustered.fasta

# Update clusters
mmseqs clusterupdate \
  /mnt/data4/planter_outputs/repseq_v3/output130/newSequenceDB \
  /tmp/unclustered_db \
  /mnt/data4/planter_outputs/repseq_v3/output130/newClusterDB \
  /mnt/data4/planter_outputs/repseq_v3/output131/updatedClusterDB \
  /tmp/clusterupdate_tmp

# Reload clusters
python scripts/load_clusters.py \
  -d /mnt/data4/master_simple.duckdb \
  -t /mnt/data4/planter_outputs/repseq_v3/output131/updatedClusterDB.tsv
```

**Pros:**
- Faster than full re-clustering
- Leverages existing cluster structure

**Cons:**
- More complex to implement
- May not be as optimal as full re-clustering

### Option 3: Document the limitation and plan for next rebuild

- Add documentation that only 8.62% of sequences are currently clustered
- Plan for full re-clustering when adding new samples
- Implement a check to ensure clustering includes all samples

## Diagnostic Scripts Created

1. **`investigate_missing_rhodiola_ugts.py`**
   - Comprehensive analysis of missing Rhodiola UGTs
   - Run with: `python investigate_missing_rhodiola_ugts.py`

## Files to Check

Key files for investigation:
- Database: `/mnt/data4/master_simple.duckdb`
- Cluster TSV: `/mnt/data4/planter_outputs/repseq_v3/output130/newClusterDB.tsv`
- Rep sequences: `/mnt/data4/planter_outputs/repseq_v3/output130/newRepSeqDB.fasta`
- Sample list: `/tmp/samples_with_duckdb_clean.txt` (91 samples)

## Next Steps

1. **Immediate:** Decide on clustering strategy (full re-cluster vs incremental update)
2. **Short-term:** Implement chosen solution
3. **Long-term:** Add validation checks to ensure clustering includes all samples going forward

## Test Case

To verify the fix works, check that these sequences are clustered:
```sql
SELECT cm.seqhash_id, cm.cluster_id
FROM cluster_members cm
WHERE cm.seqhash_id IN (
  'v1_DLS_ed0a9b8e87e9ff23685d90c4a18677202888a4ae9ef5a32411e2ae16f40d70f0.p1',
  'v1_DLS_5908d838437e1f4cddc2ce7b0c536fb73ead6f6466eb71b9ed0fab352a0413c3.p1'
);
```

Should return 2 rows after fix (currently returns 0 rows).
