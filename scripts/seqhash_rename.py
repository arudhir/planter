#!/usr/bin/env python
import argparse
from pathlib import Path
from Bio import SeqIO
import seqhash

def rename_sequences(input_file: Path, output_file: Path):
    seqrecords = list(SeqIO.parse(input_file, format='fasta'))

    renamed_seqrecords = []
    for sr in seqrecords:
        seq_hash = seqhash.seqhash(str(sr.seq))
        sr.id = seq_hash
        renamed_seqrecords.append(sr)

    SeqIO.write(renamed_seqrecords, output_file, format='fasta')

def main():
    parser = argparse.ArgumentParser(description="Rename FASTA sequences using seqhash.")
    parser.add_argument('--input', type=Path, required=True, help="Input FASTA file")
    parser.add_argument('--output', type=Path, required=True, help="Output FASTA file")

    args = parser.parse_args()

    rename_sequences(args.input, args.output)
    print(f"Renamed sequences written to {args.output}")

if __name__ == "__main__":
    main()
