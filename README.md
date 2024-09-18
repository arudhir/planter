# planter: work in progress

# vision
![full](images/planter.png "planter")

we want to assemble and annotate transcripts

# quickstart

These two commands will do a transcriptome assembly and annotation for Mesoplasma florum:
```console
$ make image  # Build the docker image
$ docker-compose run --rm planter \
    snakemake \
        --cores 16 \
        --config samples="SRR12068547" \
        outdir="outputs" \
        s3_bucket="recombia.planter"  # Run the pipeline
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

`fastp`![https://github.com/OpenGene/fastp] is an all-in-one read preprocessing tool. It does adapter trimming, quality filtering, and gets quality statistics.

### 4. Filter rRNA reads with `bbduk`

The rRNA sequences are obtained from the SILVA database [https://www.arb-silva.de/download/arb-files/]. We remove reads that come from rRNA.

### 5. Normalize read coverage with `bbnorm`

We even out the read coverage to make assembly more efficient.

**n.b.** This normalization is done for assembly, not for read quantification, where we would of course want to retain uneven coverage.

### 6. Get final read statistics with `fastp`

Just a sanity check.

## Assembly



## Annotation

## Expression

# snakemake dag
![dag](images/dag.png "dag")


You can get a quick summary the pipeline's outputs:
```console
$ jq < output/SRR8053131/SRR8053131_stats.json
{
  "Number of Orthologs": 42707,
  "Number of Secondary Metabolite Genes": 945,
  "Percent Mapped": 88.21,
  "Number of Reads Mapped": 79206071.0,
  "Total Reads (Before)": 182858938,
  "Total Reads (After)": 180926144,
  "Sequencing Type": "paired end (50 cycles + 50 cycles)",
  "Number of Transcripts": 126768,
  "Total Assembly Length": 110773877,
  "N50": 1372,
  "ExN50": 980,
  "Average Expression (TPM)": 7.89
}
```

We see SRR8053131 has 42,707 orthologs and 945 secondary metabolite genes.


# TODO

- [ ] Write script that finalizes the output metadata
    - [ ] Number of reads
    - [ ] Number of rRNA reads
    - [ ] Number of contigs
    - [ ] Number of TransDecoder predicted ORFs
    - [ ] Number of TransDecoder predicted proteins
    - [ ] Number of EggNOG Orthologs
    - [ ] Distribution of COG annotations
    - [ ] Percentage of reads mapped to the SPAdes transcripts
- [ ] Create database of transcripts
    - [ ] Generate clustered representative transcripts
- [ ] Create Flask app to serve the homology search

# Test Samples

SRR5936537: https://www.ncbi.nlm.nih.gov/sra/?term=SRR5936537
SRR8053131: https://www.ncbi.nlm.nih.gov/sra/?term=SRR8053131
SRR29142729: https://www.ncbi.nlm.nih.gov/sra/?term=SRR29142729
SRR22420515: https://www.ncbi.nlm.nih.gov/sra/?term=SRR22420515


Notes:

- SRR2103848 implies that there are some contaminants. i notice insect orthologs mixed in with the flower. we should screen with sourmash to figure out what we ought to use to filter. or perhaps another way to filter reads. that should make assemblies faster
- rnaspades bc its slow. but we'd need to probably normalize
- sourmash --> download reference genome --> filter?

Downloading Viridiplantae, green plants
```bash
datasets download genome taxon 33090 --include gtf,cds,protein,rna --reference
```

SourMash
```bash
find signatures -name '*.sig.gz' | xargs sourmash gather SRR2103848_trinity.Trinity.fasta.sig -o gather_results.csv
```

1. Download reads
2. Trimmomatic -- I am curious whether default trimming parameters are the best
3. Mash to plant reference genomes
4. Download suspected reference genome, filter reads that match the reference
5. Normalize?
6. 

# Random Installation Notes

- You need to install SPAdes from source to resolve the error code -11 thing: https://github.com/ablab/spades/issues/1297


# Test Enzymes

SRR5936537|: Genbank AUI41117.1|A0A2I6B3N5|Rhodiola rosea
SRR8053131|Piper methysticum: Genbank QCX36371.1|A0A384E132|Piper methysticum
SRR29142729|Gerbera hybrid cultivar: Genbank QCX36376.1|A0A4Y5QR90|Piper methysticum
SRR22420515: Genbank BBI55602.1|BBI55602.1|




# NOTE

SRR12068547 is Mesoplasma. JC got it!

# TODO

1. Work out how to get the other members of a cluster
2. Get the eggNOG annotations as part of the search output
3. Attach the metadata for the SRA ID for the seqhash hit
