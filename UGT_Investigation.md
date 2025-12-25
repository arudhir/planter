# UGT Investigation Lab Notebook

**Date:** 2025-10-13
**Investigator:** Claude
**Goal:** Find why Rhodiola rosea UGT sequences are missing from cluster representatives

---

## Investigation Plan

1. Survey all samples in `/mnt/data4/recombia.planter` (using `*RR*` glob for SRA IDs)
2. Verify each sample has a DuckDB file
3. Use `./planter/scripts/iterative_clustering.py` for clustering
4. Use MMSeqs2 to search for Rhodiola UGT hits and analyze their sizes
5. Audit `/mnt/data4/master.duckdb`

---

## Step 1: Sample Survey & Database Audit

### Finding all SRA samples

**Command:** `ls -d /mnt/data4/recombia.planter/*RR* | wc -l`
**Result:** 155 samples found on disk ✓

### Auditing master.duckdb

**Command:** `duckdb /mnt/data4/master.duckdb "SELECT COUNT(DISTINCT sample_id) FROM sequences"`
**Result:** 115 samples in database

**Discrepancy:** 155 samples on disk vs 115 in master.duckdb

### Missing Samples Analysis

**Total missing from master.duckdb:** 40 samples

**Missing samples with .duckdb files ready:** 31 samples
```
ERR12954282, ERR12954283, ERR12954285, SRR11002815, SRR11002816,
SRR14338344, SRR14338345, SRR18735291, SRR25582083, SRR25582084,
SRR25582086, SRR307776, SRR307777, SRR307779, SRR314553,
SRR5007079, SRR5557829, SRR5557833, SRR5557843, SRR5557844,
SRR5557845, SRR5557848, SRR5557852, SRR5557859, SRR5557860,
SRR5557863, SRR5557864, SRR5557869, SRR5557873, SRR5557876,
SRR8053131
```

**Missing samples WITHOUT .duckdb files:** 9 samples
```
SRR128115, SRR128116, SRR128117, SRR128118, SRR128119,
SRR128121, SRR13329606, SRR1660509, SRR1660513
```

### Rhodiola Samples Status

**In master.duckdb:**
- ✓ SRR5936536 (Rhodiola rosea)
- ✓ SRR5936537 (Rhodiola rosea)
- ✗ SRR22844166 (Rhodiola rosea) - **MISSING!**

### File Requirements (from Snakefile)

For database creation, each sample needs:
1. `.pep` file in `transdecoder/` directory
2. `.emapper.annotations` file in `eggnog/` directory
3. `.duckdb` file (created by pipeline)

**Verified for ERR12954282:**
- ✓ Has .pep file (15 MB)
- ✓ Has .emapper.annotations (12 MB)
- ✓ Has .duckdb file (100 MB)

### Key Finding #1

**31 samples have complete data and .duckdb files but are NOT in master.duckdb**

This includes SRR22844166 (one of the Rhodiola samples), which explains why we're missing some Rhodiola data!

---

## Step 2: Comparing master.duckdb vs master_simple.duckdb

### master.duckdb (production - with problems)
- Samples: 115
- Rhodiola samples: 2/3 (missing SRR22844166)

### master_simple.duckdb (your rebuild last night)
- Samples: 76
- Rhodiola samples: 3/3 ✓ (ALL present!)
- Missing: 79 samples

### Summary Table

| Database | Total Samples | Rhodiola Complete? | Missing Samples |
|----------|--------------|-------------------|-----------------|
| Disk | 155 | N/A | N/A |
| master.duckdb | 115 | NO (2/3) | 40 |
| master_simple.duckdb | 76 | YES (3/3) | 79 |

### Key Finding #2

**master_simple.duckdb has ALL 3 Rhodiola samples, but master.duckdb is missing SRR22844166!**

This confirms that the rebuild approach works for including Rhodiola data, but we need a complete rebuild with all 155 samples.

---

## Step 3: Checking Rhodiola UGT Clustering in master_simple.duckdb

### Query: Are Rhodiola UGTs in cluster_members?

**Result:** ❌ **NO** - All Rhodiola UGT sequences have NULL cluster_id

Example sequences checked:
- `v1_DLS_7d3cf135d9c7105733e7b7305ef9308ad5320dfac19544e98aae1349bf2f4e85.p1` (616 aa, UGT80A2-like) - NOT clustered
- `v1_DLS_ed0a9b8e87e9ff23685d90c4a18677202888a4ae9ef5a32411e2ae16f40d70f0.p1` (499 aa, UGT85A1) - NOT clustered

### Clustering Coverage in master_simple.duckdb

- Total sequences: 1,970,849
- Sequences with clusters: 169,863
- **Coverage: 8.62%** (same as master.duckdb!)

### Key Finding #3: THE REAL PROBLEM

**Both master.duckdb AND master_simple.duckdb loaded the SAME old cluster data!**

The problem is NOT just missing samples in the database. The problem is:

1. ✓ Samples were added to database (sequences + annotations loaded)
2. ✓ master_simple rebuild included all 3 Rhodiola samples
3. ❌ **Clustering was NEVER re-run** - old cluster data from 169K sequences was reused
4. ❌ The Rhodiola UGT sequences (and 1.8M other sequences) were never clustered

**This means:**
- The cluster_members table is incomplete and outdated
- It only covers ~8.62% of sequences in the database
- The Rhodiola UGTs exist in the database but are "invisible" to cluster-based searches
- We need to **re-run MMSeqs2 clustering** on ALL 1.97M sequences, not just reload old clusters

---

## Step 4: Understanding the Correct Clustering Workflow

### After reviewing `iterative_cluster.py` and `mmseqs_cluster_update.py`:

The workflow operates on **.pep files** (protein sequences), NOT DuckDB files!

**Critical Understanding:**
- MMSeqs2 works with FASTA files (`.pep` files contain protein sequences)
- The DuckDB files are for STORAGE and QUERYING, not for clustering input
- Clustering happens FIRST, then results are loaded INTO DuckDB

### Correct Workflow (3 Steps)

#### Step 1: Prepare Input (Extract .pep files from all 155 samples)

```bash
# Create list of all .pep files
find /mnt/data4/recombia.planter/*RR*/transdecoder -name "*.pep" | \
  grep -v "transdecoder_dir" | \
  sort > /tmp/all_pep_files.txt

# Verify count
wc -l /tmp/all_pep_files.txt
# Should be 155 (or close - some samples might be missing .pep files)
```

#### Step 2: Run Iterative Clustering with MMSeqs2

```bash
# Run iterative clustering on ALL samples
python planter/scripts/iterative_cluster.py \
  -g "/mnt/data4/recombia.planter/*/transdecoder/*.pep" \
  -o /mnt/data4/planter_outputs/repseq_v4 \
  --database /mnt/data4/master_complete.duckdb

# This will:
# 1. Take first .pep file as initial representatives
# 2. Iteratively add each subsequent .pep file
# 3. Use mmseqs clusterupdate to maintain clusters
# 4. Generate final newClusterDB.tsv with ALL sequences
# 5. Load clusters into database automatically (if --database provided)
```

**What this does internally:**
- For each `.pep` file (iteration):
  - `mmseqs createdb` - convert FASTA to MMSeqs format
  - `mmseqs clusterupdate` - add new sequences to existing clusters
  - Prefers **longest sequences** as representatives (`--cluster-mode 2 --cov-mode 1`)
  - Updates cluster membership
- Final output: `newClusterDB.tsv` with cluster assignments for ALL sequences

#### Step 3: Create Complete Database

**Option A: Start fresh (RECOMMENDED)**
```bash
# 1. Create complete master database from all 155 samples
python smart_merge.py  # (modify to use all 155 samples)

# 2. Load clusters from v4 output
python scripts/load_clusters.py \
  -d /mnt/data4/master_complete.duckdb \
  -t /mnt/data4/planter_outputs/repseq_v4/output154/newClusterDB.tsv
```

**Option B: Use existing master_simple.duckdb**
```bash
# Just reload with NEW cluster data
python scripts/load_clusters.py \
  -d /mnt/data4/master_simple.duckdb \
  -t /mnt/data4/planter_outputs/repseq_v4/output154/newClusterDB.tsv
```

### Key Points About the Workflow

1. **Input is .pep files**, NOT .duckdb files
2. **Clustering is separate from database**, happens on protein sequences
3. **iterative_cluster.py handles everything**:
   - Calls mmseqs_cluster_update.py for each iteration
   - Uses `mmseqs clusterupdate` (NOT `mmseqs cluster` from scratch)
   - Can optionally load results into DuckDB at the end
4. **Longest sequences become representatives** (`--cov-mode 1`)
5. **Output is newClusterDB.tsv** - a 2-column file: `representative_id\tmember_id`

### Why Your Current Approach Is Wrong

❌ **Wrong:** "Concatenate duckdb; cluster_members empty → MMSeqs2 → populate cluster_members"

✓ **Correct:** "Run MMSeqs2 on .pep files → Get TSV → Load TSV into DuckDB cluster_members"

The DuckDB database is the END RESULT, not the input for clustering!

---

## Step 5: Ready-to-Run Commands

### Verification: Check .pep file availability
```bash
find /mnt/data4/recombia.planter/*RR*/transdecoder -name "*.pep" 2>/dev/null | \
  grep -v "transdecoder_dir" | wc -l
```
**Result:** 131 .pep files available (out of 155 samples)

### Command 1: Run Complete Clustering (LONG RUNNING - Hours to days)

```bash
cd /home/ubuntu/planter

# Run iterative clustering on all 131 samples with .pep files
python planter/scripts/iterative_cluster.py \
  -g "/mnt/data4/recombia.planter/*/transdecoder/*.pep" \
  -o /mnt/data4/planter_outputs/repseq_v4

# Do NOT provide --database flag yet
# We'll load clusters AFTER creating a complete merged database
```

**Expected output:**
- Directory: `/mnt/data4/planter_outputs/repseq_v4/output130/` (or output131, depending on iterations)
- Key file: `newClusterDB.tsv` - contains ALL cluster assignments
- Representative sequences: `newRepSeqDB.fasta`

**Time estimate:** Several hours to complete (depends on sequence count)

### Command 2: Create Complete Master Database

**Option A: Merge all 155 samples (requires modifying smart_merge.py)**
```bash
# Create list of all samples with .duckdb files
for sample in $(ls -d /mnt/data4/recombia.planter/*RR*/); do
  sample_id=$(basename $sample)
  if [ -f "$sample/$sample_id.duckdb" ]; then
    echo $sample_id
  fi
done | sort > /tmp/all_samples_with_duckdb.txt

wc -l /tmp/all_samples_with_duckdb.txt
# Then modify smart_merge.py to use this list
```

**Option B: Use existing master_simple.duckdb (76 samples, has all Rhodiola)**
```bash
# This already has all 3 Rhodiola samples
# Just needs new cluster data loaded
cp /mnt/data4/master_simple.duckdb /mnt/data4/master_v2.duckdb
```

### Command 3: Load New Cluster Data

```bash
# After clustering completes, load the TSV
python scripts/load_clusters.py \
  -d /mnt/data4/master_v2.duckdb \
  -t /mnt/data4/planter_outputs/repseq_v4/output130/newClusterDB.tsv

# Expected output:
#   Coverage: ~100% (or close - all sequences that were in .pep files will be clustered)
```

### Command 4: Verify Rhodiola UGTs Are Now Clustered

```bash
duckdb /mnt/data4/master_v2.duckdb "
  SELECT
    a.seqhash_id,
    a.preferred_name,
    s.length,
    cm.cluster_id
  FROM annotations a
  JOIN sequences s ON a.seqhash_id = s.seqhash_id
  LEFT JOIN sra_metadata m ON a.sample_id = m.sample_id
  LEFT JOIN cluster_members cm ON a.seqhash_id = cm.seqhash_id
  WHERE (a.preferred_name LIKE '%UGT%' OR a.description LIKE '%UGT%')
    AND m.organism = 'Rhodiola rosea'
  ORDER BY s.length DESC
  LIMIT 10
"

# cluster_id should NOT be NULL anymore!
```

---

## Summary

### The Root Causes Identified

1. **master.duckdb missing samples:** 40 samples on disk not included (including Rhodiola SRR22844166)
2. **Outdated cluster data:** Both databases loaded old cluster TSV covering only 169K sequences
3. **Never re-clustered:** When samples were added, clustering wasn't re-run on the complete dataset

### The Solution

1. ✓ Run `iterative_cluster.py` on ALL .pep files (131 samples)
2. ✓ Create/use complete master database with all samples
3. ✓ Load NEW cluster TSV into database
4. ✓ Verify Rhodiola UGTs are now discoverable via cluster membership

### Next Steps After Fix

Test that the Rhodiola UGT query (using AUI41147 sequence) now returns:
- Cluster representative IDs
- Other Rhodiola sequences in the same clusters
- Similar functionality to NCBI BLAST results

