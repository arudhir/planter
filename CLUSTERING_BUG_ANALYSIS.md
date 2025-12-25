# Clustering Bug Analysis & Fix

**Date:** 2025-10-16
**Issue:** Sequences added to database after initial clustering remain unclustered indefinitely
**Impact:** ~91% of sequences in production database have no cluster assignments
**Status:** Root cause identified, fix implemented

---

## Executive Summary

A critical bug in the incremental clustering workflow caused newly added samples to be excluded from future clustering runs. This resulted in the vast majority of sequences (1.8M out of 2.0M) having no cluster assignments, making them invisible to BLAST searches against representative sequences.

**Root Cause:** The `extract_representative_sequences()` function only extracted sequences where `repseq_id = seqhash_id`, excluding all unclustered sequences (`repseq_id IS NULL`).

**Fix:** Modified the extraction query to include both existing representatives AND unclustered sequences.

---

## Problem Description

### What Users Observed

When searching for Rhodiola rosea UGT sequences:
- Sequences existed in the database ✓
- Sequences had proper annotations ✓
- But BLAST searches didn't find them ✗
- No cluster assignments (`cluster_id = NULL`) ✗
- Not marked as representatives (`is_representative = FALSE`) ✗

### Database State

**Production database (`master.duckdb`):**
```
Total sequences:        1,970,849
Clustered sequences:      169,863
Clustering coverage:         8.62%
Unclustered sequences:  1,800,986  (91.38%)
```

**Rhodiola UGT sequences specifically:**
```
Total Rhodiola UGT sequences:     102
Full-length UGTs (>400aa):         44
UGTs with cluster assignments:      0  ⬅️ THE PROBLEM
```

---

## Root Cause Analysis

### The Incremental Clustering Workflow

The Snakemake workflow (`planter/workflow/rules/finalize.smk`) implements incremental clustering:

```python
rule mmseqs_clustering:
    # 1. Fetch master.duckdb from S3
    # 2. Extract representative sequences
    # 3. Run mmseqs clusterupdate with new samples
    # 4. Load updated clusters back into database
```

### The Bug

In `planter/database/utils/duckdb_utils.py`, the `extract_representative_sequences()` function (line 910):

```python
def extract_representative_sequences(
    db_path: Union[str, Path],
    output_path: Union[str, Path]
) -> None:
    """Extracts representative sequences from a DuckDB database."""

    query = """
        SELECT seqhash_id, sequence
        FROM sequences
        WHERE repseq_id = seqhash_id;  ⬅️ BUG: Only existing representatives
    """
```

### Why This Breaks Incremental Clustering

**Timeline of events:**

1. **Initial clustering** (early 2025):
   - Clustered ~76 samples
   - Created 169,863 cluster assignments
   - All sequences got `repseq_id` values

2. **Rhodiola samples added** (March-April 2025):
   - Individual `.duckdb` files created ✓
   - Merged into `master.duckdb` ✓
   - Sequences added with `repseq_id = NULL` (no cluster info yet)

3. **Next clustering run**:
   - `extract_representative_sequences()` called
   - Query: `WHERE repseq_id = seqhash_id`
   - **Rhodiola sequences excluded** (because `repseq_id IS NULL`)
   - MMSeqs2 only compares new samples against OLD representatives
   - Rhodiola never enters the comparison pool

4. **Result**:
   - Rhodiola sequences stay with `repseq_id = NULL` forever
   - They're never extracted for future clustering runs
   - They accumulate as "orphaned" sequences

### Verification of Root Cause

Checked extraction on current database:

```sql
-- Current (buggy) query
SELECT COUNT(*) FROM sequences WHERE repseq_id = seqhash_id;
-- Result: 157,506 (only old representatives)

-- Including unclustered sequences
SELECT COUNT(*) FROM sequences WHERE repseq_id = seqhash_id OR repseq_id IS NULL;
-- Result: 1,958,492 (representatives + unclustered)
```

**The bug excluded 1.8M sequences from clustering!**

---

## Impact Assessment

### Affected Systems

1. **BLAST Search**
   - Representative sequence FASTA built from `extract_representative_sequences()`
   - Missing 91% of sequences
   - Users can't find Rhodiola UGTs or any other recently added sequences

2. **Cluster-based Queries**
   - Web app queries that filter by cluster membership
   - Sequences with `cluster_id = NULL` excluded from results
   - Data discoverability severely impaired

3. **Data Completeness**
   - Only 8.62% of sequences have cluster information
   - Breaks assumption that all sequences are clustered

### Samples Affected

Any sample added AFTER the initial clustering run (output130):
- Rhodiola samples (SRR5936536, SRR5936537, SRR22844166)
- Approximately 55 other samples
- ~1.8M individual sequences

### Why It Wasn't Caught Earlier

1. **Silent Failure**: Incremental clustering appeared to work (no errors)
2. **Partial Success**: Some sequences from new samples DID cluster (if similar to existing representatives)
3. **Search Still Worked**: Annotation-based search worked fine (uses `annotations` table)
4. **Low Visibility**: Users primarily searched by annotation, not by sequence similarity

---

## The Fix

### Code Change

**File:** `planter/database/utils/duckdb_utils.py`
**Function:** `extract_representative_sequences()` (line 930)

```python
# BEFORE (buggy)
query = """
    SELECT seqhash_id, sequence
    FROM sequences
    WHERE repseq_id = seqhash_id;
"""

# AFTER (fixed)
query = """
    SELECT seqhash_id, sequence
    FROM sequences
    WHERE repseq_id = seqhash_id           -- Existing representatives
       OR repseq_id IS NULL;                -- Unclustered sequences (new samples)
"""
```

### Why This Works

**Next clustering run after fix:**

1. Extract representatives:
   - Gets existing representatives (repseq_id = seqhash_id) ✓
   - Gets unclustered sequences (repseq_id IS NULL) ✓

2. MMSeqs2 clustering:
   - Compares new sample against ALL sequences (old + unclustered)
   - Rhodiola sequences compared against everything
   - If similar → cluster together
   - If dissimilar → become singleton representatives

3. Result:
   - All sequences get cluster assignments
   - No orphaned sequences
   - True hands-off incremental clustering

### Immediate Workaround

The current full re-clustering (131 samples from scratch) bypasses the bug entirely:
- Processes all `.pep` files directly
- Doesn't rely on `extract_representative_sequences()`
- Every sequence gets clustered

---

## Testing

### Test Case: 10-Sample Subset

Ran clustering on 10 samples including all 3 Rhodiola samples:

```bash
python planter/scripts/iterative_cluster.py \
  -g "/tmp/test_clustering_peps/*.pep" \
  -o /tmp/test_clustering_output
```

**Results:**
- Total sequences: 283,341
- Clustered: 83,081 (29.3%)
- Singletons: 60,811 (73.2% of clustered sequences)
- Time: 13:41 minutes

**Key Finding:** MMSeqs2 DOES create singleton entries in the TSV as `seq_id seq_id`, so the workflow design was correct. The bug was purely in the extraction function.

### Verification Tests

Created `/tmp/verify_clustering_test.py` to validate:
1. ✅ Representatives are longest sequences in clusters (91.6% success rate)
2. ✅ Singletons are included as self-representatives
3. ✅ All clustered sequences appear in TSV

---

## Migration Plan

### Phase 1: Full Re-clustering (In Progress)

Running complete re-clustering on all 131 samples:

```bash
python planter/scripts/iterative_cluster.py \
  -g "/mnt/data4/recombia.planter/*/transdecoder/*.pep" \
  -o /mnt/data4/planter_outputs/repseq_v6
```

**Estimated time:** 3-4 hours
**Output:** `/mnt/data4/planter_outputs/repseq_v6/output130/newClusterDB.tsv`

### Phase 2: Build New Master Database

```bash
# Option A: Build from all 131 samples
python build_master_131.py

# Option B: Use existing master_simple.duckdb (76 samples)
# This already has all Rhodiola samples

# Load new cluster data
python scripts/load_clusters.py \
  -d /mnt/data4/master_v6.duckdb \
  -t /mnt/data4/planter_outputs/repseq_v6/output130/newClusterDB.tsv
```

### Phase 3: Deploy Fix

Apply the code fix to `extract_representative_sequences()` and deploy to production.

### Phase 4: Validation

After deployment, verify:

```sql
-- Check clustering coverage
SELECT
    COUNT(*) as total_sequences,
    COUNT(CASE WHEN repseq_id IS NOT NULL THEN 1 END) as clustered,
    ROUND(COUNT(CASE WHEN repseq_id IS NOT NULL THEN 1 END) * 100.0 / COUNT(*), 1) as coverage_pct
FROM sequences;

-- Expected: >95% coverage

-- Check Rhodiola UGTs specifically
SELECT
    a.sample_id,
    COUNT(*) as total_ugts,
    COUNT(CASE WHEN cm.cluster_id IS NOT NULL THEN 1 END) as clustered_ugts
FROM annotations a
JOIN sequences s ON a.seqhash_id = s.seqhash_id
LEFT JOIN cluster_members cm ON a.seqhash_id = cm.seqhash_id
WHERE a.sample_id IN ('SRR5936536', 'SRR5936537', 'SRR22844166')
  AND (a.preferred_name LIKE '%UGT%' OR a.description LIKE '%glucosyltransferase%')
GROUP BY a.sample_id;

-- Expected: All UGTs should have cluster assignments
```

---

## Future Recommendations

### 1. Add Monitoring

Add a database health check to detect unclustered sequences:

```python
def check_clustering_health(db_path):
    """Alert if >10% of sequences are unclustered."""
    con = duckdb.connect(db_path)
    result = con.execute("""
        SELECT
            COUNT(*) as total,
            COUNT(CASE WHEN repseq_id IS NULL THEN 1 END) as unclustered,
            COUNT(CASE WHEN repseq_id IS NULL THEN 1 END) * 100.0 / COUNT(*) as pct
        FROM sequences
    """).fetchone()

    if result[2] > 10:
        logger.warning(f"ALERT: {result[2]:.1f}% of sequences are unclustered!")

    return result
```

### 2. Add Unit Tests

```python
def test_extract_representative_sequences_includes_unclustered():
    """Ensure unclustered sequences are included in extraction."""
    # Create test database with unclustered sequences
    # Extract representatives
    # Assert unclustered sequences are in output
```

### 3. Periodic Full Re-clustering

Consider running full re-clustering:
- Every 50-100 new samples
- Or quarterly
- To prevent accumulation of edge cases

### 4. Add Validation to Snakemake

Add a checkpoint after clustering:

```python
rule validate_clustering:
    input:
        cluster_tsv = rules.mmseqs_clustering.output.cluster_file,
        master_db = rules.update_database.output.done
    run:
        # Verify all sequences in database are in cluster TSV
        con = duckdb.connect(input.master_db)

        # Check for unclustered sequences
        unclustered = con.execute("""
            SELECT COUNT(*) FROM sequences WHERE repseq_id IS NULL
        """).fetchone()[0]

        if unclustered > 0:
            logger.warning(f"Found {unclustered} unclustered sequences after clustering!")
```

---

## Related Documentation

- **Investigation Details:** `RHODIOLA_UGT_INVESTIGATION_SUMMARY.md`
- **Lab Notebook:** `UGT_Investigation.md`
- **Test Scripts:** `/tmp/verify_clustering_test.py`, `/tmp/generate_meeting_stats.sh`
- **Snakemake Rules:** `planter/workflow/rules/finalize.smk`
- **Database Utils:** `planter/database/utils/duckdb_utils.py`

---

## Lessons Learned

1. **Silent failures are dangerous**: The clustering appeared to work but was silently excluding data
2. **Test edge cases**: Initial implementation worked for first batch but failed on incremental updates
3. **Monitor data completeness**: Should have caught the 8.62% coverage metric earlier
4. **NULL handling matters**: The difference between `repseq_id = seqhash_id` and `repseq_id IS NULL` was critical
5. **Full integration tests needed**: Unit tests alone wouldn't have caught this workflow bug

---

## Acknowledgments

This bug was discovered during investigation of why Rhodiola rosea UGT sequences were not appearing in BLAST search results, despite being present in the database with proper annotations.

The investigation revealed a systemic issue affecting all samples added after the initial clustering run, not just Rhodiola samples.
