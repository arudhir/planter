# Clustering Workflow: Issues and Solutions

## Current Problems

### Problem 1: Wrong Representatives Selected (FIXED)
**Issue:** MMSeqs2 was selecting SHORTEST sequences as cluster representatives instead of longest.

**Root Cause:**
- `mmseqs cluster` had correct parameters: `--cluster-mode 2 --cov-mode 1` ✅
- `mmseqs clusterupdate` was missing these parameters ❌

**Fix Applied:** Added `--cluster-mode 2 --cov-mode 1` to `clusterupdate` in `planter/scripts/mmseqs_cluster_update.py`

**Status:** ✅ FIXED (commit 65d4180)

### Problem 2: Cluster Membership Data Never Imported (UNFIXED)
**Issue:** Only 232,051 sequences have cluster assignments, but 2,902,248 sequences exist in database (92% missing cluster data!)

**Root Cause:** The workflow never calls `load_clusters_from_tsv()` after clustering completes.

**Impact:**
- Cluster relationships are invisible in database queries
- Can't see which sequences belong to the same cluster
- Can't navigate from representative to members

**Status:** ❌ NEEDS FIX

### Problem 3: Incremental Clustering Limitations (ARCHITECTURAL)
**Issue:** `mmseqs clusterupdate` preserves existing representatives by design, even when longer sequences are added.

**Example:**
- Iteration 6: 296aa sequence becomes representative
- Iteration 91: 535aa sequence added to same cluster
- Result: 296aa stays as representative (wrong!)

**Root Cause:** MMSeqs2's `clusterupdate` is optimized for efficiency by preserving cluster structure.

**Status:** ❌ FUNDAMENTAL LIMITATION - needs architectural decision

---

## Current Workflow

```
┌─────────────────────────────────────────────────────────────────┐
│ 1. Build Database (per sample)                                  │
│    - Load sequences from .pep files                             │
│    - Load annotations from .emapper files                        │
│    - Load expression from quant.sf files                         │
│    - Mark all sequences as is_representative=TRUE initially      │
├─────────────────────────────────────────────────────────────────┤
│ 2. Iterative Clustering (separate process)                      │
│    - Run mmseqs_cluster_update.py for each sample iteratively   │
│    - Creates output130/newClusterDB.tsv                          │
│    - Creates output130/newRepSeqDB.fasta                         │
│                                                                  │
│    ⚠️  Issue: clusterupdate preserves old representatives        │
├─────────────────────────────────────────────────────────────────┤
│ 3. ??? (MISSING STEP)                                           │
│    - Should load cluster data into database                      │
│    - Should update is_representative flags                       │
│                                                                  │
│    ❌ This step never happens!                                   │
├─────────────────────────────────────────────────────────────────┤
│ 4. Web Search                                                    │
│    - Searches representatives only (by default)                  │
│    - Returns incomplete results due to wrong representatives     │
│    - Can't show cluster relationships (data missing)             │
└─────────────────────────────────────────────────────────────────┘
```

---

## Proposed Solutions

### Option A: Fix Current Workflow (Incremental Updates)
**Keep iterative clustering but fix the gaps**

#### Changes Required:

1. **After clustering completes, automatically load cluster data:**

```python
# In planter/scripts/iterative_cluster.py (add at end of main())
def main():
    # ... existing code ...

    # After clustering completes
    if args.load_to_database:
        logger.info("Loading cluster data into database")
        from planter.database.builder import SequenceDBBuilder

        cluster_tsv = os.path.join(base_dir, f"output{len(files)}", "newClusterDB.tsv")

        with SequenceDBBuilder(args.database_path, output_dir=args.output_dir) as builder:
            builder.load_clusters_from_tsv(cluster_tsv)

        logger.info("✓ Cluster data loaded into database")
```

2. **Add database path argument to iterative_cluster.py:**

```python
parser.add_argument(
    "-db", "--database-path",
    help="Path to DuckDB database (if provided, will load clusters after completion)"
)
```

3. **Accept that representatives might not be optimal:**
   - Due to incremental updates, some representatives may be suboptimal
   - Trade-off: Speed vs Accuracy

**Pros:**
- Fast incremental updates
- Minimal changes to existing code
- Works for adding new samples

**Cons:**
- Representatives may not be longest sequences (architectural limitation)
- Example: 296aa stays as rep even when 535aa is added

---

### Option B: Periodic Full Re-clustering
**Use incremental updates day-to-day, but periodically re-cluster from scratch**

#### Changes Required:

1. **Use Option A for daily updates** (fast, but suboptimal)

2. **Weekly/Monthly: Run full re-clustering:**

```bash
# Use the new full_recluster.py script
./planter/scripts/full_recluster.py \
  -g "/mnt/data4/recombia.planter/*/transdecoder/*.pep" \
  -o /mnt/data4/planter_outputs/repseq_full
```

3. **Load the optimal clusters:**

```bash
python scripts/load_clusters.py \
  -d /mnt/data4/master.duckdb \
  -t /mnt/data4/planter_outputs/repseq_full/clusters.tsv
```

**Pros:**
- Best of both worlds: fast updates + optimal representatives
- Guaranteed longest sequences as representatives after full re-clustering

**Cons:**
- Need to run full re-clustering periodically (takes hours)
- Two clustering approaches to maintain

---

### Option C: Always Full Re-clustering
**Abandon incremental updates entirely**

#### Changes Required:

1. **Replace iterative_cluster.py usage with full_recluster.py:**

```bash
# Instead of:
./planter/scripts/iterative_cluster.py -g "*.pep" -o /output

# Use:
./planter/scripts/full_recluster.py -g "*.pep" -o /output
```

2. **After clustering, load into database:**

```bash
python scripts/load_clusters.py \
  -d /mnt/data4/master.duckdb \
  -t /output/clusters.tsv
```

**Pros:**
- Guaranteed optimal representatives (longest sequences)
- Simpler workflow (one clustering method)
- No architectural limitations

**Cons:**
- Takes 3-4 hours to cluster 131 samples from scratch
- Must re-cluster all sequences when adding even 1 new sample
- Not scalable for frequent updates

---

## Recommended Solution: **Option B** (Hybrid Approach)

### Rationale:
- Daily/weekly: Use incremental updates for speed (Option A)
- Monthly: Run full re-clustering for optimal representatives
- Best balance of speed and accuracy

### Implementation Steps:

#### 1. Fix iterative clustering to load cluster data

```bash
# Modify planter/scripts/iterative_cluster.py
# Add --database-path argument
# Add automatic load_clusters_from_tsv() at end
```

#### 2. Add full re-clustering script (already created)

```bash
# Already exists: planter/scripts/full_recluster.py
# Use for periodic optimization
```

#### 3. Add cluster loading script (already created)

```bash
# Already exists: scripts/load_clusters.py
# Use after any clustering operation
```

#### 4. Update documentation and workflow

```bash
# Daily workflow (fast):
./planter/scripts/iterative_cluster.py \
  -g "/mnt/data4/recombia.planter/*/transdecoder/*.pep" \
  -o /mnt/data4/planter_outputs/repseq_daily \
  --database-path /mnt/data4/master.duckdb

# Monthly workflow (optimal):
./planter/scripts/full_recluster.py \
  -g "/mnt/data4/recombia.planter/*/transdecoder/*.pep" \
  -o /mnt/data4/planter_outputs/repseq_monthly

python scripts/load_clusters.py \
  -d /mnt/data4/master.duckdb \
  -t /mnt/data4/planter_outputs/repseq_monthly/clusters.tsv
```

---

## Immediate Actions Needed

### To Fix Current Database:

```bash
# Load the existing cluster data (even if suboptimal)
python scripts/load_clusters.py \
  -d /mnt/data4/master.duckdb \
  -t /mnt/data4/planter_outputs/repseq_v3/output130/newClusterDB.tsv
```

This will:
- Add cluster membership data for all 2.9M sequences
- Update is_representative flags correctly
- Enable cluster relationship queries

### To Get Optimal Representatives:

```bash
# Run full re-clustering from scratch
tmux new -s full_recluster

./planter/scripts/full_recluster.py \
  -g "/mnt/data4/recombia.planter/*/transdecoder/*.pep" \
  -o /mnt/data4/planter_outputs/repseq_v4

# After completion, load into database
python scripts/load_clusters.py \
  -d /mnt/data4/master.duckdb \
  -t /mnt/data4/planter_outputs/repseq_v4/clusters.tsv
```

This will ensure:
- Longest sequences are representatives (535aa for Rhodiola cluster)
- All cluster relationships are in database
- Optimal search results

---

## Testing the Fix

After loading cluster data, verify with:

```bash
# Check cluster membership is loaded
duckdb /mnt/data4/master.duckdb -c "
SELECT
    (SELECT COUNT(*) FROM sequences) as total_sequences,
    (SELECT COUNT(*) FROM cluster_members) as clustered_sequences,
    (SELECT COUNT(*) FROM clusters) as total_clusters
"

# Check Rhodiola cluster specifically
duckdb /mnt/data4/master.duckdb -c "
SELECT s.seqhash_id, s.length, s.is_representative, cm.cluster_id
FROM sequences s
LEFT JOIN cluster_members cm ON s.seqhash_id = cm.seqhash_id
WHERE s.seqhash_id LIKE '%810d761d368e1bd880fb7abf%'  -- 535aa
   OR s.seqhash_id LIKE '%12a38ce2a433f0e8a119290564%'  -- 308aa
   OR s.seqhash_id LIKE '%2d9f22e9b819e2baba8337c974%'  -- 296aa
ORDER BY s.length DESC
"
```

Expected result with full re-clustering:
- 535aa sequence: is_representative=TRUE
- 308aa sequence: is_representative=FALSE, same cluster_id as 535aa
- 296aa sequence: is_representative=FALSE, same cluster_id as 535aa

---

## Questions to Answer:

1. **How often do you add new samples?**
   - Daily → Option A (incremental with periodic full re-clustering)
   - Weekly/Monthly → Option C (always full re-clustering)

2. **How important is representative optimality?**
   - Critical → Option C or B
   - Nice to have → Option A

3. **How long can you wait for clustering?**
   - Must be fast → Option A
   - Can wait hours → Option C
   - Hybrid → Option B (recommended)

4. **Do you need cluster relationships in database queries?**
   - Yes → Must load cluster data after every clustering run
   - No → Can skip cluster loading (but limits functionality)
