






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
