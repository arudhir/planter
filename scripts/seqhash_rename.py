#!/usr/bin/env python
import sys
from pathlib import Path
from Bio import SeqIO
import seqhash

fasta_file = Path(sys.argv[1])
seqrecords = list(SeqIO.parse(fasta_file, format='fasta'))

renamed_seqrecords = []
for sr in seqrecords:
    seq_hash = seqhash.seqhash(str(sr.seq))
    sr.id = seq_hash
    renamed_seqrecords.append(sr)

out_filename = f'{fasta_file.stem}_renamed.fa'
SeqIO.write(renamed_seqrecords, out_filename, format='fasta')
