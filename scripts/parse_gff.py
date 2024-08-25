#!/usr/bin/env python3

# Read in the GFF file using BioPython
import Bio
from BCBio import GFF
from collections import defaultdict

def parse_gff(file_path):
    with open(file_path) as handle:
        for rec in GFF.parse(handle):
            yield rec

# Example usage
if __name__ == "__main__":
    gff_file = "../tests/test-output/transdecoder/SRR29142729_transdecoder/longest_orfs.gff3"
    
    # Collect the records into a defaultdict(list)
    # Aggregate on the gene level by splitting on ~~
    records = defaultdict(list)
    for rec in parse_gff(gff_file):
        for feature in rec.features:
            if feature.type == "gene":
                gene_id = feature.id.split("~~")[0]
                records[gene_id].append(feature)
    