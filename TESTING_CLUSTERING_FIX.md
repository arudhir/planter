# Testing the Clustering Bug Fix

**Date:** 2025-10-16
**Related:** CLUSTERING_BUG_ANALYSIS.md
**Status:** All tests passing ✓

---

## Overview

This document describes the comprehensive test suite created to verify the clustering bug fix. The fix addresses the issue where unclustered sequences (repseq_id IS NULL) were excluded from incremental clustering runs.

**Bug:** `extract_representative_sequences()` only extracted sequences WHERE repseq_id = seqhash_id
**Fix:** Modified to extract WHERE repseq_id = seqhash_id OR repseq_id IS NULL
**Impact:** Ensures 1.8M unclustered sequences (91% of database) get included in future clustering

---

## Test Strategy

The test suite uses a **two-layer approach**:

1. **Unit Tests** - Verify the SQL query fix in isolation
2. **Integration Tests** - Verify the complete incremental clustering workflow

### Why This Approach?

- **Unit tests** catch regressions in the specific function
- **Integration tests** ensure the fix works in the real workflow context
- **Both** reproduce the exact bug scenario that caused Rhodiola UGTs to remain unclustered

---

## Unit Tests

**Location:** `tests/database/test_extract_representatives.py`
**Tests:** 8 tests, all passing ✓

### Test Coverage

#### 1. `test_extract_includes_unclustered_sequences` ⭐ PRIMARY TEST
**Purpose:** Verify the fix includes unclustered sequences (repseq_id IS NULL)

**Test Data:**
- 3 existing representatives (repseq_id = seqhash_id)
- 2 clustered members (repseq_id points to a rep)
- 5 unclustered sequences (repseq_id IS NULL) ← **THE FIX**

**Assertion:**
```python
# Expected: 3 reps + 5 unclustered = 8 sequences
assert len(extracted_ids) == 8
assert 'new1' in extracted_ids  # Unclustered sequence MUST be included
```

**What it tests:** The core fix - that unclustered sequences are now extracted

---

#### 2. `test_extract_excludes_clustered_members`
**Purpose:** Verify clustered members (non-representatives) are still excluded

**Assertion:**
```python
# These should NOT be extracted (they have representatives already)
assert 'member1' not in extracted_ids
assert 'member2' not in extracted_ids
```

**What it tests:** The fix doesn't break existing behavior

---

#### 3. `test_old_behavior_would_fail`
**Purpose:** Demonstrate the OLD buggy behavior

**Test:**
```sql
-- Old buggy query
SELECT seqhash_id FROM sequences WHERE repseq_id = seqhash_id;
-- Returns: Only 3 reps (excludes 5 unclustered) ← THE BUG
```

**What it tests:** Documents the exact bug that was fixed

---

#### 4. `test_new_behavior_includes_unclustered`
**Purpose:** Demonstrate the NEW fixed behavior

**Test:**
```sql
-- New fixed query
SELECT seqhash_id FROM sequences
WHERE repseq_id = seqhash_id OR repseq_id IS NULL;
-- Returns: 3 reps + 5 unclustered = 8 ✓
```

**What it tests:** The fix correctly includes unclustered sequences

---

#### 5. `test_extract_sequences_are_valid`
**Purpose:** Verify extracted sequences have valid sequence data

**What it tests:** Data integrity after extraction

---

#### 6. `test_extract_empty_database`
**Purpose:** Edge case - extracting from empty database doesn't crash

**What it tests:** Robustness

---

#### 7. `test_extract_all_unclustered`
**Purpose:** Edge case - all sequences unclustered (first clustering run)

**Test Data:** 10 sequences, all with repseq_id = NULL

**Assertion:**
```python
# Should extract all 10 sequences
assert len(extracted_ids) == 10
```

**What it tests:** First-time clustering scenario works

---

#### 8. `test_clustering_coverage_metric`
**Purpose:** Verify the health check metric from CLUSTERING_BUG_ANALYSIS.md

**Test:**
```sql
SELECT
    COUNT(*) as total,
    COUNT(CASE WHEN repseq_id IS NULL THEN 1 END) as unclustered,
    ROUND(COUNT(CASE WHEN repseq_id IS NULL THEN 1 END) * 100.0 / COUNT(*), 1) as coverage_pct
FROM sequences
```

**What it tests:** Monitoring capability to detect this bug in the future

---

## Integration Tests

**Location:** `tests/clustering/test_incremental_clustering.py`
**Tests:** 5 tests, all passing ✓

### Test Coverage

#### 1. `test_incremental_workflow_bug_reproduction` ⭐ PRIMARY TEST
**Purpose:** Reproduce the EXACT bug scenario

**Workflow:**
1. **Phase 1:** Initial state - 6 sequences clustered (output130)
2. **Phase 2:** Add Rhodiola samples - 5 new sequences with repseq_id = NULL
3. **Phase 3:** Extract representatives
4. **Phase 4:** Verify extraction includes BOTH reps AND unclustered

**Critical Assertions:**
```python
# Rhodiola UGTs MUST be in extraction (the fix)
assert 'rhodiola_ugt1' in extracted_ids
assert 'rhodiola_ugt2' in extracted_ids
assert 'rhodiola_ugt3' in extracted_ids
```

**What it tests:** The exact scenario that caused the production bug

---

#### 2. `test_old_bug_behavior`
**Purpose:** Demonstrate OLD bug excluded Rhodiola sequences

**Test:**
```sql
-- Old buggy query
SELECT seqhash_id FROM sequences WHERE repseq_id = seqhash_id
```

**Assertion:**
```python
# THE BUG: Rhodiola sequences were EXCLUDED
assert 'rhodiola_ugt1' not in old_extracted_ids
assert 'rhodiola_ugt2' not in old_extracted_ids
```

**What it tests:** Documents the production bug behavior

---

#### 3. `test_full_incremental_clustering_cycle` ⭐ COMPREHENSIVE TEST
**Purpose:** Test complete end-to-end clustering cycle

**Workflow:**
1. Start with clustered database (6 sequences)
2. Add new Rhodiola samples (2 sequences, repseq_id = NULL)
3. Extract representatives (includes unclustered after fix)
4. Simulate MMSeqs2 clustering output
5. Load clusters into database
6. **Verify:** ALL sequences now have cluster assignments

**Critical Assertions:**
```python
# ALL sequences should be clustered (no orphans)
assert clustered == 8  # All 8 sequences
assert unclustered == 0  # Zero unclustered

# Rhodiola UGTs have cluster assignments
assert ugt1[1] == 'rhodiola_ugt1'  # ugt1 is representative
assert ugt2[1] == 'rhodiola_ugt1'  # ugt2 clusters with ugt1

# cluster_members table includes Rhodiola
assert len(rhodiola_members) == 2
```

**What it tests:**
- Complete workflow works end-to-end
- No sequences left orphaned
- Database state is correct after clustering

---

#### 4. `test_clustering_health_check`
**Purpose:** Test the monitoring metric for clustering health

**Test:**
1. Start with 100% clustered (healthy)
2. Add 50 unclustered sequences (simulating bug)
3. Check health metric

**Assertion:**
```python
# Should detect 89.3% unclustered (unhealthy!)
assert unclustered_pct == 89.3
if unclustered_pct > 10.0:
    # Alert threshold exceeded
    assert True  # Health check works
```

**What it tests:**
- Monitoring can detect this bug in production
- Alert threshold (10%) is appropriate

---

#### 5. `test_multiple_incremental_rounds`
**Purpose:** Test multiple sequential clustering rounds

**Workflow:**
- **Round 1:** Add batch A (2 sequences) → cluster → update
- **Round 2:** Add batch B (1 sequence) → cluster → update
- **Verify:** All sequences eventually get clustered

**What it tests:**
- Incremental clustering works over multiple iterations
- No cumulative degradation
- True "hands-off" incremental clustering

---

## Running the Tests

### Run All Tests
```bash
# Both unit and integration tests
pytest tests/database/test_extract_representatives.py tests/clustering/test_incremental_clustering.py -v

# Total: 13 tests, all passing ✓
```

### Run Unit Tests Only
```bash
pytest tests/database/test_extract_representatives.py -v

# 8 tests, ~2 seconds
```

### Run Integration Tests Only
```bash
pytest tests/clustering/test_incremental_clustering.py -v

# 5 tests, ~2 seconds
```

### Run Specific Test
```bash
# Run the primary bug reproduction test
pytest tests/clustering/test_incremental_clustering.py::test_incremental_workflow_bug_reproduction -v -s
```

---

## Test Results

**Date Run:** 2025-10-16
**Environment:** Linux 6.8.0-1023-aws, Python 3.12.0, DuckDB 1.1.3

```
tests/database/test_extract_representatives.py::test_extract_includes_unclustered_sequences PASSED
tests/database/test_extract_representatives.py::test_extract_excludes_clustered_members PASSED
tests/database/test_extract_representatives.py::test_extract_sequences_are_valid PASSED
tests/database/test_extract_representatives.py::test_old_behavior_would_fail PASSED
tests/database/test_extract_representatives.py::test_new_behavior_includes_unclustered PASSED
tests/database/test_extract_representatives.py::test_extract_empty_database PASSED
tests/database/test_extract_representatives.py::test_extract_all_unclustered PASSED
tests/database/test_extract_representatives.py::test_clustering_coverage_metric PASSED

tests/clustering/test_incremental_clustering.py::test_incremental_workflow_bug_reproduction PASSED
tests/clustering/test_incremental_clustering.py::test_old_bug_behavior PASSED
tests/clustering/test_incremental_clustering.py::test_full_incremental_clustering_cycle PASSED
tests/clustering/test_incremental_clustering.py::test_clustering_health_check PASSED
tests/clustering/test_incremental_clustering.py::test_multiple_incremental_rounds PASSED

======================== 13 passed in 4.06s ========================
```

**Status:** ✅ All tests passing

---

## What the Tests Verify

### ✅ The Fix Works
- Unclustered sequences (repseq_id IS NULL) are now extracted
- Rhodiola UGTs will be included in future clustering runs
- No sequences get permanently orphaned

### ✅ No Regressions
- Existing representatives still extracted correctly
- Clustered members still excluded (as intended)
- Edge cases handled (empty database, all unclustered)

### ✅ Complete Workflow
- End-to-end incremental clustering works
- Database state correct after clustering
- Multiple rounds of incremental clustering successful

### ✅ Monitoring
- Health check metric works
- Can detect high percentage of unclustered sequences
- Alert threshold (10%) is testable

---

## Future Test Additions

### Recommended Tests to Add

1. **Performance Test**
   - Test with large database (1M+ sequences)
   - Verify extraction doesn't timeout
   - Measure memory usage

2. **Concurrency Test**
   - Test concurrent extractions
   - Verify thread safety

3. **Schema Migration Test**
   - Test with old schema (no is_representative column)
   - Verify backward compatibility

4. **Real Data Test**
   - Use actual Rhodiola samples
   - Verify against production database state

---

## Test Maintenance

### When to Update Tests

- **Schema changes:** Update test database creation
- **Query changes:** Update SQL assertions
- **New edge cases:** Add new test cases
- **Performance issues:** Add performance tests

### How to Debug Failing Tests

1. **Unit test failure:**
   - Check SQL query in `extract_representative_sequences()`
   - Verify test database setup
   - Compare expected vs actual IDs

2. **Integration test failure:**
   - Check each phase separately
   - Verify MMSeqs2 simulation
   - Check database state after each step

3. **Flaky test:**
   - Check for file cleanup issues
   - Verify no shared state between tests
   - Add more assertions

---

## Key Takeaways

1. **Test the bug scenario exactly:** The integration tests reproduce the EXACT workflow that caused the production bug

2. **Test both old and new behavior:** Documenting the old buggy behavior helps prevent regressions

3. **Test the complete workflow:** Unit tests verify the fix, integration tests verify it works in context

4. **Test monitoring capabilities:** The health check tests ensure we can detect this bug class in the future

5. **Edge cases matter:** Empty database, all unclustered, multiple rounds - all tested

---

## Related Documentation

- **Bug Analysis:** CLUSTERING_BUG_ANALYSIS.md
- **Code Fix:** `planter/database/utils/duckdb_utils.py:930`
- **Snakemake Workflow:** `planter/workflow/rules/finalize.smk`
- **Investigation:** RHODIOLA_UGT_INVESTIGATION_SUMMARY.md
