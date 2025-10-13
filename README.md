# planter: work in progress

# vision
![full](images/planter.png "planter")

we want to assemble and annotate transcripts

# quickstart

## requirements

* ~12gb ram
* docker & docker-compose


These two commands will do a transcriptome assembly and annotation for Mesoplasma florum:
```console
$ make image  # Build the docker image
$ docker-compose run --rm planter \
    snakemake \
        --cores 16 \
        --config samples="SRR12068547" \
        outdir="outputs" \
        s3_bucket=$S3_BUCKET  # Run the pipeline
```

Example output structure:
```console
├── SRR12068547_stats.json
├── eggnog
│   ├── SRR12068547.emapper.annotations
│   ├── SRR12068547.emapper.annotations.xlsx
│   └── plots
│       ├── cog_category_counts.csv
│       └── cog_distribution_plot.png
├── illumina
│   ├── processed                                       # Trimmed, rRNA filtered, normalized; used for assembly
│   │   ├── SRR12068547.1.fq.gz                 
│   │   ├── SRR12068547.2.fq.gz
│   │   ├── SRR12068547_fastp.html
│   │   └── SRR12068547_fastp.json
│   └── rrna_filtered                                   # Trimmed, rRNA filtered; used for quant
│       ├── SRR12068547.1.fq.gz
│       ├── SRR12068547.2.fq.gz
│       ├── SRR12068547_rRNA_filter.log
│       └── SRR12068547_rRNA_filter.stats
├── quants
│   ├── SRR12068547.quant.tsv                           # NumReads (for DESeq2) and TPM counts
=├── rnaspades
│   ├── SRR12068547_transcripts_renamed.fasta           # Headers renamed to seqhashes
│   ├── transcripts.fasta                            
└── transdecoder
    ├── SRR12068547.bed
    ├── SRR12068547.cds
    ├── SRR12068547.gff3
    ├── SRR12068547.pep                                 # Predicted proteins; used for eggNOG
```
# Workflow

## Read Processing

`workflow/rules/reads.smk`

### 1. Download the reads using `fastq-dump`

**n.b** The SRA experiment needs to be paired end. There is a command-line utility in `scripts/get_srr_metadata.py` that will get the metadata for SRA IDs:

```console
$ ./scripts/get_srr_metadata.py --help
usage: get_srr_metadata.py [-h] --srr SRR [SRR ...]

Retrieve comprehensive SRA information for one or multiple SRR
IDs

options:
  -h, --help           show this help message and exit
  --srr SRR [SRR ...]  One or more SRR IDs to look up
```

### 2. Compress the reads with `pigz`

Fast compression to `*.fastq.gz`

### 3. Preprocess the reads with `fastp`

[`fastp`](https://github.com/OpenGene/fastp) is an all-in-one read preprocessing tool. It does adapter trimming, quality filtering, and gets quality statistics.

### 4. Filter rRNA reads with `bbduk`

The rRNA sequences are obtained from the [SILVA database](https://www.arb-silva.de/download/arb-files/). We remove reads that come from rRNA.

### 5. Normalize read coverage with `bbnorm`

We even out the read coverage to make assembly more efficient.

**n.b.** This normalization is done for assembly, not for read quantification, where we would of course want to retain uneven coverage.

### 6. Get final read statistics with `fastp`

Just a sanity check.

## Assembling Genes

`workflow/rules/assembly.smk`

### 1. Assemble transcripts with `rnaSPAdes`

We opt for [`rnaSPAdes`](https://gensoft.pasteur.fr/docs/SPAdes/3.14.0/rnaspades_manual.html) over [`Trinity`](https://github.com/trinityrnaseq/trinityrnaseq) for speed and consistency. Modern reviews of RNAseq assemblers show it having comprable performance.

**n.b.** This is why we explicitly normalize the reads before `rnaSPAdes`. `Trinity` would normalize as part of the pipeline.

**n.b.** We also rename the output file from the generic `transcripts.fasta` to `{SRR_ID}_transcripts.fasta` for clarity.

### 2. Rename headers with `seqhash`

The FASTA headers from `rnaSPAdes` contain metadata information about the assembled transcript itself. For example:

`>NODE_1_length_11323_cov_82.795022_g0_i0`

Although this has useful information, we want something that:

1. Uniquely identifies the sequence
2. Is a stable identifier

To this end, we use [`seqhash`](https://keonigandall.com/posts/seqhash.html). The headers will then look like:

`>v1_DLS_f4afc4cb2e95c54f170cc1fe8bd33bb73b12b378ddccde59a53559c84a6d05bf NODE_1_length_11323_cov_82.795022_g0_i0`

i.e., `{seqhash} {original_header}`. Downstream tools utilize the seqhash to identify the sequence.

### 3. Predict ORFs with `TransDecoder`

[`TransDecoder`](https://github.com/TransDecoder/TransDecoder/wiki) will predict ORFs in the transcripts. It generates four files: a .bed file, a .gff3 file, a .cds (fasta) file, and a .pep (fasta) file. We use the .pep file for annotations.

**n.b.** According to the GitHub, `TransDecoder` is no longer maintained as of March, 2024. This may pose an issue in the future.

## Annotation

`workflow/rules/annotate.smk`

### 1. Annotate with `eggnog-mapper`

Predicted proteins (i.e. the `*.pep` file from `TransDecoder`) are functionally annotated with [`eggnog-mapper`](https://github.com/eggnogdb/eggnog-mapper) using precomputed orthologs from the [eggNOG database](http://eggnog5.embl.de/#/app/home).

`eggnog-mapper` requires pre-downloaded eggNOG databases. We have downloaded these already an EBS attached to the cloud instance located in `/mnt/data`. When running the containerized pipeline, our `docker-compose.yml` mounts this volume.

`eggnog-mapper` provides a script to download the databases. We downloaded MMSeqs2, HMMER, PFAM, and DIAMOND (for novel families) databases with:
```console
/opt/eggnog-mapper/download_eggnog_data.py -M -H -F -P --data_dir /mnt/data
```

The primary output of `eggnog-mapper` is a `{SRR_ID}.emapper.annotations` tab-delimited file and also `{SRR_ID}.emapper.annotations.xlsx` Excel file with the same contents.

There are many fields in the annotation file. A full description can be found on the [eggNOG wiki](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v1#project_nameemapperannotations-file). I think there are a few columns that are of particular interest to us:

- `query_name`: The seqhash for the transcript
- `seed_ortholog`: Best predicted protein match in eggNOG
- `Description`: The functional description
- `Preferred_name`: The "standard" gene name
- `COG_category`: The COG (clusters of orthologous groups)category of the ortholog. You can find descriptions of these [here](https://www.ncbi.nlm.nih.gov/research/cog#).
- `GOs`: The gene ontology terms
- `EC`: The EC number
- `KEGG_Pathway, KEGG_Module, KEGG_Reaction, KEGG_rclass`: KEGG identifiers
- `PFAMs`: Pfam domains

### 2. Analyze the annotations with `parse_eggnog.py`

A quick Python script to count the number of hits for each COG category and visualize the results.

## Quantification

`workflow/rules/quant.smk`

### 1. Quantify expression with `salmon`

We primarily use this for QC-purposes. We will use the RNAspades transcripts and the un-normalized reads.

## Finalizing and uploading the results

`workflow/rules/finalize.smk`

### 1. QC Stats

We gather some key statistics from our various analyses and upload the results to S3.

For example, for our example Mesoplasma sample, we get:

```console
$ jq < SRR12068547_stats.json 
{
  "Number of Orthologs": 251,
  "Number of Secondary Metabolite Genes": 0,
  "Percent Mapped": 96.8,
  "Number of Reads Mapped": 6771889.0,
  "Total Reads (Before)": 16089680,
  "Total Reads (After)": 14926556,
  "Sequencing Type": "paired end (50 cycles + 50 cycles)",
  "Number of Transcripts": 3894,
  "Total Assembly Length": 1570556,
  "N50": 386,
  "ExN50": 290,
  "Average Expression (TPM)": 256.81
}
```

If we use a more interesting sample, such as SRR14292007, which is RNAseq of a lichen strain *Cladonia macilenta* that produces biruloquinone, we get:
```
$ jq < SRR14292007/SRR14292007_stats.json 
{
  "Number of Orthologs": 16763,
  "Number of Secondary Metabolite Genes": 771,
  "Percent Mapped": 97.91,
  "Number of Reads Mapped": 50725894.0,
  "Total Reads (Before)": 104107964,
  "Total Reads (After)": 103881978,
  "Sequencing Type": "paired end (100 cycles + 100 cycles)",
  "Number of Transcripts": 15757,
  "Total Assembly Length": 44849063,
  "N50": 4125,
  "ExN50": 7708,
  "Average Expression (TPM)": 63.46
}
```

We see that eggNOG identifies 16763 orthologs with 771 of them being known secondary metabolite genes! Maybe worth looking into.

### 3. Upload to S3

We then zip the output directory and upload to an S3 bucket.

# snakemake dag
![dag](images/dag.png "dag")

## Creating and Updating the Transcriptome Database

We would like to iteratively collect transcripts we assemble into a non-redundant database we can use for downstream homology searches. For this, a utility script is provided in `scripts/ mmseqs_cluster_update.py` that takes in two inputs: the old FASTA file of transcripts and the new file we want to update it with.

```console
$ python scripts/mmseqs_cluster_update.py --help
usage: mmseqs_cluster_update.py [-h] --old OLD --new NEW [-o OUTPUT_DIR] [-t TMP_DIR]

Update MMSeqs2 clusters with new sequences.

options:
  -h, --help            show this help message and exit
  --old OLD             Path to the old representative sequences FASTA file.
  --new NEW             Path to the new sequences FASTA file to add.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory to store outputs.
  -t TMP_DIR, --tmp_dir TMP_DIR
                        Directory for temporary files.
```

This script wraps `mmseqs clusterupdate` and associated commands and statistics to allow easy updating of the transcriptome database.

**n.b.** We are using `mmseqs clusterupdate` to ensure the stability of cluster assignments and avoid redundant clustering.

### Critical Discovery: MMSeqs2 Representative Selection Bug

**Problem Identified (October 2025)**

Through database analysis and testing, we discovered that **65% of cluster representatives were the SHORTEST sequences in their clusters**, not the longest. This caused major issues:

- Full-length sequences were hidden from searches (only representatives are searched)
- Truncated fragments were returned instead of complete sequences
- NCBI BLAST returned different results because it searches all sequences

**Root Cause Analysis**

Investigation revealed the issue was in how MMSeqs2 selects cluster representatives:

1. **Default behavior (`--cluster-mode 0`)**: Uses greedy set cover algorithm that picks the sequence with the MOST alignments as representative, regardless of length
   - Our unit tests proved this actually prefers SHORTER sequences
   - Database analysis showed 65% of representatives were the shortest in their cluster

2. **Initial fix attempt (`--cluster-mode 2` only)**: Added greedy clustering by sequence length
   - This STILL didn't work correctly
   - Created multiple small clusters instead of grouping similar sequences

3. **Complete fix (`--cluster-mode 2 + --cov-mode 1`)**:
   - `--cluster-mode 2`: Sorts sequences by decreasing length, clusters longest first
   - `--cov-mode 1`: Calculates coverage on TARGET sequence only, allowing shorter sequences to join longer sequences' clusters
   - **This combination ensures the longest sequence is ALWAYS selected as representative**

**How We Discovered This**

1. **User reported issue**: Searching for 498 aa Rhodiola UGT returned 241 aa truncated sequence
2. **Database analysis**: Wrote `tests/database/test_rhodiola_ugt_debug.py` that showed:
   - Only 38.5% of cluster representatives were the longest sequence
   - 64.9% were actually the SHORTEST sequence (opposite of desired!)
   - Some clusters had 261 aa rep when 4,250 aa sequence was available
3. **Unit tests**: Created `tests/clustering/test_mmseqs_representative_selection.py` with synthetic data proving:
   - Cluster-mode 0 (default) selects SHORTEST sequences
   - Cluster-mode 2 alone creates too many clusters
   - Cluster-mode 2 + cov-mode 1 correctly selects longest

**The Fix**

Updated `planter/scripts/mmseqs_cluster_update.py` line 132:
```python
# OLD (broken)
"--cluster-mode", "2"

# NEW (correct)
"--cluster-mode", "2",  # Greedy clustering by sequence length
"--cov-mode", "1",      # Select longest sequence as representative
```

**Verification**

Run the unit tests to verify correct behavior:
```bash
python tests/clustering/test_mmseqs_representative_selection.py
```

Or with pytest:
```bash
pytest tests/clustering/test_mmseqs_representative_selection.py -v -s
```

The tests use synthetic sequences of known lengths to prove:
- ✅ Full-length sequences selected over truncated versions
- ✅ Longest sequence chosen when multiple lengths present
- ✅ Rhodiola-like scenario (5 full-length + 1 truncated) correctly handled
- 🔴 Cluster-mode 0 incorrectly selects shortest sequences

**Re-clustering Required**

Existing databases must be re-clustered with the corrected parameters:
```bash
python ./planter/scripts/iterative_cluster.py \
  -g "/path/to/samples/*/transdecoder/*.pep" \
  -o /path/to/output/repseq_fixed
```

## Homology Search

Flask app in `app/`, `python app/main.py`

![mmseqs-search](images/mmseqs-search.png "mmseqs-search")
![mmseqs-search-results](images/mmseqs-search-results.png "mmseqs-search-results")

# Database

![simplified-vision](images/simplified-vision.png "simplified-vision")
![schema](images/schema.png "schema")

# Tests

The project includes comprehensive test suites for various components:

## Database Tests

### Core Database Builder (`tests/database/test_database_builder.py`)

Tests for the `SequenceDBBuilder` class which handles database construction from raw sequence data:

- **Database Initialization Tests**:
  - Verifies proper creation of all required tables (sequences, annotations, go_terms, etc.)
  - Tests correct loading of sequence data, annotation data, GO terms, EC numbers
  - Validates expression data loading and linking
  - Verifies gene-protein mapping relationships

- **Database Summary Tests**:
  - Tests the database summary functionality that reports statistics
  - Validates counts of total sequences, samples, annotations, and various annotation types

### Expression Queries (`tests/database/test_expression_queries.py`) 

Tests for the `QueryManager` class, specifically focusing on expression data queries:

- **Expression Summary Tests**:
  - Tests summarizing expression data across samples
  - Validates expression level categorization (low, medium, high, very high)
  - Tests retrieval of top-expressed sequences with proper ordering

- **Expression Data Retrieval Tests**:
  - Tests fetching expression data for specific gene/protein sequences
  - Validates linking between annotation and expression data
  - Tests expression level statistics and distribution analysis

- **Sequence Search Tests**:
  - Tests searching sequences by various criteria (description, length, sample)
  - Validates handling of non-existent data in search functions

### Schema Migrations (`tests/database/test_schema_migrations.py`)

Tests for the `SchemaManager` class which handles database schema migrations:

- **Migration Application Tests**:
  - Verifies migrations are applied in the correct sequence
  - Tests that all required tables are created by migrations
  - Validates the structure of tables created by migrations

- **Schema Validation Tests**:
  - Tests schema of expression table (columns, types)
  - Verifies foreign key constraints between tables
  - Validates data integrity across the database

### Sequence Search (`tests/database/test_sequence_search.py`)

Tests for the web interface's sequence search functionality:

- **API Endpoint Tests**:
  - Tests the `load_example` endpoint for loading example sequences
  - Verifies proper error handling when files are not found

### Database Utility Tests (`tests/database/utils/`)

#### DuckDB Utilities (`test_duckdb_utils.py`)

Tests for database operations in the utilities module:

- **Database Merging Tests**:
  - Tests merging of multiple DuckDB databases
  - Verifies all tables and data are correctly combined
  - Validates data integrity after merging

- **Cluster Update Tests**:
  - Tests updating sequence cluster information from MMSeqs2 TSV output
  - Verifies correct representative sequence assignment
  - Validates cluster member relationships

- **Integration Tests**:
  - Tests with real sample data from fixtures
  - End-to-end MMSeqs2 clustering integration
  - Tests full workflow from clustering to database update

#### S3 Utilities (`test_s3_utils.py`)

Tests for S3 interaction functionality:

- **Archive Creation Tests**:
  - Tests creating zip archives of output directories

- **S3 Upload Tests**:
  - Tests successful upload to S3
  - Verifies handling of files that already exist in S3
  - Tests error handling during uploads

## Pipeline Tests

### DuckDB Creation Pipeline (`tests/pipeline/test_duckdb_creation.py`)

Tests for the database creation pipeline:

- **Database Creation Tests**:
  - Tests the end-to-end process of creating a DuckDB database
  - Verifies all required tables and data are loaded correctly
  - Tests integration with fixture data from real samples

- **Expression Data Pipeline Tests**:
  - Tests loading expression data from JSON files into the database
  - Validates queries for expression statistics
  - Tests distribution analysis and top expressed sequence retrieval

## Clustering Tests

### MMSeqs2 Representative Selection (`tests/clustering/test_mmseqs_representative_selection.py`)

Unit tests that verify MMSeqs2 clustering correctly selects the longest sequence as representative:

- **Simple Truncation Test**: Full-length vs truncated sequence
- **Multiple Length Test**: 5 sequences of different lengths (500, 400, 300, 200, 100 aa)
- **Rhodiola Scenario Test**: Realistic case with 5 full-length + 1 truncated sequence
- **Cluster Mode Comparison**: Demonstrates difference between modes 0, 2, and 3

Run these tests to verify clustering behavior:
```bash
python tests/clustering/test_mmseqs_representative_selection.py
```

### Rhodiola UGT Debug Tests (`tests/database/test_rhodiola_ugt_debug.py`)

Tests that analyze the actual database to identify clustering issues:

- **Database Analysis**: Checks cluster statistics across all clusters
- **Baseline vs V2 Comparison**: Compares search results between clustering versions
- **All Sequences Search**: Tests searching against all sequences vs just representatives

These tests revealed the 65% shortest-sequence representative problem.

## Running Tests

Run all tests with:
```bash
python -m pytest tests/
```

Run specific test files:
```bash
python -m pytest tests/database/test_database_builder.py
```

Run a specific test:
```bash
python -m pytest tests/database/test_database_builder.py::TestSequenceDBBuilder::test_database_initialization
```

Run tests with increased verbosity:
```bash
python -m pytest tests/ -v
```

# TODO

- [ x ] Write script that finalizes the output metadata
    - [ x ] Number of reads
    - [ x ] Number of rRNA reads
    - [ x ] Number of contigs
    - [ ] Number of TransDecoder predicted ORFs
    - [ ] Number of TransDecoder predicted proteins
    - [ ] Number of EggNOG Orthologs
    - [ ] Distribution of COG annotations
    - [ x ] Percentage of reads mapped to the SPAdes transcripts
- [ x ] Create database of transcripts
    - [ x ] Generate clustered representative transcripts
    - [ x ] Come up with S3 scheme
- [ x ] Create Flask app to serve the homology search

## Misc Notes

- SRR5936537: https://www.ncbi.nlm.nih.gov/sra/?term=SRR5936537
- SRR8053131: https://www.ncbi.nlm.nih.gov/sra/?term=SRR8053131
- SRR29142729: https://www.ncbi.nlm.nih.gov/sra/?term=SRR29142729
- SRR22420515: https://www.ncbi.nlm.nih.gov/sra/?term=SRR22420515

### Downloading Viridiplantae, green plants
```bash
datasets download genome taxon 33090 --include gtf,cds,protein,rna --reference
```

SourMash
```bash
find signatures -name '*.sig.gz' | xargs sourmash gather SRR2103848_trinity.Trinity.fasta.sig -o gather_results.csv
```

### Random Installation Notes

- You need to install SPAdes from source to resolve the error code -11 thing: https://github.com/ablab/spades/issues/1297


### Test Enzymes

- SRR5936537|: Genbank AUI41117.1|A0A2I6B3N5|Rhodiola rosea
- SRR8053131|Piper methysticum: Genbank QCX36371.1|A0A384E132|Piper methysticum
- SRR29142729|Gerbera hybrid cultivar: Genbank QCX36376.1|A0A4Y5QR90|Piper methysticum
- SRR22420515: Genbank BBI55602.1|BBI55602.1|

- SRR12068547 is Mesoplasma. JC got it!
