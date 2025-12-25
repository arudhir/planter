# Rhodiola UGT Search Discrepancy - Root Cause Analysis

## Problem Statement

When searching for Rhodiola rosea UGT enzyme AUI41147 (498 aa):
- **NCBI BLAST** returns 5 full-length hits: AUI41147, AUI41146, AUI41143, AUI41132, AUI41144
- **MMSeqs2/Planter** returns 1 truncated hit: AUI41142 (240 aa)

## Root Cause

**The clustering algorithm is selecting truncated sequence fragments as cluster representatives instead of full-length sequences.**

## Investigation Findings

### 1. Query Sequence Not in Database
- AUI41147 (498 aa) does **not exist** in the database
- None of the expected NCBI hits (AUI41147, AUI41146, AUI41143, AUI41132, AUI41144) are present
- This explains why we can't find exact matches

### 2. Database Contains Many Rhodiola Sequences
- **Total Rhodiola rosea sequences:** 81,779
- **Truncated sequences (230-250 aa):** 3,405
  - 232 marked as representatives (6.8%)
- **Full-length UGT-sized (475-510 aa):** 3,010
  - 112 marked as representatives (3.7%)

### 3. MMSeqs2 Search Results
Running the actual search with AUI41147:
- **73 total hits** found
- **Only 2 Rhodiola hits**, both TRUNCATED:
  - `v1_DLS_38e612b34e68d0fc6d2620a80d0b20bbee0e86c1aa057e221db11bb3caafb8cb.p2` (241 aa, 55.6% identity)
  - `v1_DLS_60e8c75f5cf1e6d418a596f084f21654458e3f6ea937c476cd0f17c55c0436a8.p1` (110 aa, 45.1% identity)
- **ZERO full-length Rhodiola UGTs found**

### 4. Cluster Analysis Reveals the Problem

**Example: 241 aa representative cluster**
```
Cluster members (4 total):
- 539 aa (full-length!)
- 241 aa (representative) ← Truncated sequence chosen as rep
- 240 aa
- 240 aa
```

**The clustering selected the 241 aa truncated sequence as representative even though a 539 aa full-length sequence exists in the same cluster!**

### 5. Full-Length UGTs with Truncated Representatives

Analyzing 20 full-length Rhodiola UGTs (475-510 aa):
- **6 out of 20 (30%)** have truncated representatives (< 400 aa)
- Examples:
  - 510 aa sequence → 244 aa representative
  - 510 aa sequence → 201 aa representative
  - 510 aa sequence → 241 aa representative
  - 510 aa sequence → 214 aa representative
  - 510 aa sequence → 370 aa representative

**These full-length sequences are invisible to MMSeqs2 search because only representatives are searched.**

## Why This Differs from NCBI BLAST

| Aspect | NCBI BLAST | Planter/MMSeqs2 |
|--------|-----------|-----------------|
| Database | All sequences | Only cluster representatives |
| Coverage | Comprehensive | Reduced (by clustering) |
| Full-length bias | None | No preference for full-length |
| Sensitivity | High | Depends on representative selection |

## Recommendations

### Immediate Fixes

1. **Fix Representative Selection in Clustering**
   - Prefer longer sequences as representatives
   - Add length-based weighting to representative selection
   - Options:
     - Re-cluster with `--cluster-mode 2` (longest sequence as rep)
     - Use `--cov-mode 0` (coverage of representative)
     - Add post-processing to swap truncated reps with full-length members

2. **Remove `.dropna()` in Merge Logic** (app/main.py:825)
   ```python
   # BEFORE (drops results missing cluster info):
   merged_df = merged_df.merge(cluster_df, on='target', how='left').dropna()

   # AFTER:
   merged_df = merged_df.merge(cluster_df, on='target', how='left')
   # Re-enable the fillna logic (lines 828-836)
   ```

3. **Search Against Full Database Option**
   - Provide option to search full sequences, not just representatives
   - Trade-off: slower but more comprehensive

### Long-Term Solutions

1. **Add Quality Control for Representatives**
   - Flag clusters where representative is significantly shorter than members
   - Report: "X% of clusters have truncated representatives"
   - Automated re-selection of representatives based on length

2. **Multi-tier Search Strategy**
   - First: Search representatives (fast)
   - If few results: Search full database (comprehensive)
   - Or: Always search cluster members of top hits

3. **Documentation**
   - Add warning in UI that results may miss sequences due to clustering
   - Display cluster size and suggest checking cluster members
   - Link to view all cluster members for each hit

## Test Files Created

1. `tests/database/test_rhodiola_ugt_debug.py` - Database queries to analyze the issue
2. `test_search_aui41147.py` - Run actual MMSeqs2 search and analyze results

## Commands to Reproduce

```bash
# Run database analysis
uv run pytest tests/database/test_rhodiola_ugt_debug.py -v -s

# Run MMSeqs2 search test
uv run python test_search_aui41147.py
```

## Conclusion

The discrepancy between NCBI BLAST and Planter is caused by:

1. **Clustering reducing search space** - Only representatives are searched
2. **Poor representative selection** - Truncated sequences chosen over full-length
3. **Full-length sequences hidden** - 30% of full-length UGTs have truncated representatives

The fix requires either:
- Improving clustering to prefer full-length representatives, OR
- Searching full database instead of just representatives, OR
- Expanding search to cluster members after initial representative search

**Priority:** HIGH - This affects search quality for any protein family where truncated fragments exist in the database.
