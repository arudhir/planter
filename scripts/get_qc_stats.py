#!/usr/bin/env python

from pathlib import Path
import pandas as pd
import json
import argparse
from Bio import SeqIO
from pprint import pprint

def parse_eggnog(eggnog_file):
    stats = {}
    df = pd.read_csv(eggnog_file, sep='\t', comment='#', header=None)
    df.columns = ['query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 'max_annot_lvl', 'COG_category', 'Description', 
                  'Preferred_name', 'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass', 
                  'BRITE', 'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs']
    cog_category_count = df.groupby('COG_category')['COG_category'].agg('count').to_dict()

    stats['Number of Orthologs'] = df.shape[0]
    stats['Number of Secondary Metabolite Genes'] = cog_category_count.get('Q', 0)
    return stats

def parse_salmon(salmon_qc):
    salmon_output = json.load(open(salmon_qc, 'r'))
    stats = {}
    stats['Percent Mapped'] = round(float(salmon_output['percent_mapped']), 2)
    stats['Number of Reads Mapped'] = float(salmon_output['num_mapped'])
    return stats

def parse_fastp(fastp_file):
    fastp_output = json.load(open(fastp_file, 'r'))
    stats = {}
    stats['Total Reads (Before)'] = fastp_output['summary']['before_filtering']['total_reads']
    stats['Total Reads (After)'] = fastp_output['summary']['after_filtering']['total_reads']
    stats['Sequencing Type'] = fastp_output['summary']['sequencing']
    return stats

def extract_transcript_id(description):
    """Extracts the v1_DLS contig ID from a SequenceRecord description."""
    return description.split(' ')[0]

def calculate_n50(lengths):
    """Calculates the N50 metric."""
    lengths = sorted(lengths, reverse=True)
    total_length = sum(lengths)
    running_sum = 0
    for length in lengths:
        running_sum += length
        if running_sum >= total_length / 2:
            return length

def calculate_exn50(lengths, expressions):
    """Calculates the ExN50, weighted by expression."""
    total_expression = sum(expressions)
    sorted_indices = sorted(range(len(lengths)), key=lambda i: expressions[i] * lengths[i], reverse=True)
    
    running_sum = 0
    for i in sorted_indices:
        running_sum += expressions[i] * lengths[i]
        if running_sum >= total_expression / 2:
            return lengths[i]

def get_total_assembly_length(lengths):
    """Returns the total length of all transcripts."""
    return sum(lengths)

def get_number_of_transcripts(transcripts):
    """Returns the total number of transcripts."""
    return len(transcripts)

def get_average_expression(expressions):
    """Returns the average expression (TPM) of transcripts."""
    return sum(expressions) / len(expressions) if expressions else 0

def get_transcript_stats(rnaspades_file, salmon_quantsf):
    """Gathers core transcriptome QC stats."""
    # Load SPAdes contigs and extract lengths
    transcripts = list(SeqIO.parse(rnaspades_file, 'fasta'))
    transcript_ids = [extract_transcript_id(transcript.description) for transcript in transcripts]
    lengths = [len(transcript.seq) for transcript in transcripts]

    # Load Salmon quant.sf file
    quantsf = pd.read_csv(salmon_quantsf, sep='\t')

    # Match contigs between SPAdes and Salmon outputs
    matched_lengths = []
    matched_expressions = []

    for idx, row in quantsf.iterrows():
        name = row['Name']
        if name in transcript_ids:
            index = transcript_ids.index(name)
            matched_lengths.append(lengths[index])
            matched_expressions.append(row['TPM'])

    # Collect stats
    stats = {}
    stats['Number of Transcripts'] = get_number_of_transcripts(transcripts)
    stats['Total Assembly Length'] = get_total_assembly_length(lengths)
    stats['N50'] = calculate_n50(lengths)
    
    if matched_lengths and matched_expressions:
        stats['ExN50'] = calculate_exn50(matched_lengths, matched_expressions)
        stats['Average Expression (TPM)'] = round(get_average_expression(matched_expressions), 2)
    else:
        stats['ExN50'] = 'N/A'
        stats['Average Expression (TPM)'] = 'N/A'

    return stats

def parse_arguments():
    parser = argparse.ArgumentParser(description="Get QC stats")
    parser.add_argument("--eggnog", help="Path to eggNOG annotation file", required=True)
    parser.add_argument("--salmon_metadata", help="Path to salmon quant metadata file", required=True)
    parser.add_argument("--fastp", help="Path to fastp JSON file", required=True)
    parser.add_argument("--transcripts", help="Path to rnaspades assembly file", required=True)
    parser.add_argument("--quantsf", help="Path to salmon quant.sf file", required=True)
    parser.add_argument("--output_file", help="Path to output file", required=True)
    return parser.parse_args()

def get_stats(args):
    annotation_stats = parse_eggnog(args.eggnog)
    mapping_stats = parse_salmon(args.salmon_metadata)
    read_stats = parse_fastp(args.fastp)
    transcript_stats = get_transcript_stats(args.transcripts, args.quantsf)
    stats = {**annotation_stats, **mapping_stats, **read_stats, **transcript_stats}
    return stats

def main():
    args = parse_arguments()
    stats = get_stats(args)
    # Write to output file
    with open(args.output_file, 'w') as outfile:
        json.dump(stats, outfile, indent=4)

if __name__ == "__main__":
    main()
