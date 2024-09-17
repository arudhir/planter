#!/usr/bin/env python

import argparse
import subprocess
import os
import sys
import shutil
import logging

def run_command(cmd, log_file):
    """Run a shell command and log its output."""
    logging.info(f"Running command: {' '.join(cmd)}")
    with open(log_file, 'a') as lf:
        lf.write(f"Running command: {' '.join(cmd)}\n")
        process = subprocess.Popen(cmd, stdout=lf, stderr=lf)
        process.communicate()
        if process.returncode != 0:
            sys.exit(f"Command failed: {' '.join(cmd)}")

def main():
    parser = argparse.ArgumentParser(description="Update MMSeqs2 clusters with new sequences.")
    parser.add_argument('--old', required=True, help="Path to the old representative sequences FASTA file.")
    parser.add_argument('--new', required=True, help="Path to the new sequences FASTA file to add.")
    parser.add_argument('-o', '--output_dir', default='mmseqs_output', help="Directory to store outputs.")
    parser.add_argument('-t', '--tmp_dir', default='mmseqs_tmp', help="Directory for temporary files.")
    args = parser.parse_args()

    # Set up directories
    output_dir = args.output_dir
    tmp_dir = args.tmp_dir
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)

    # Set up logging to both file and stdout
    log_file = os.path.join(output_dir, 'pipeline.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

    logging.info("Starting MMSeqs2 clustering update pipeline.")

    # Paths for databases and files
    sequenceDB = os.path.join(output_dir, 'sequenceDB')
    clusterDB = os.path.join(output_dir, 'clusterDB')
    addedSequenceDB = os.path.join(output_dir, 'addedSequenceDB')
    allSequenceDB = os.path.join(output_dir, 'allSequenceDB')
    newSequenceDB = os.path.join(output_dir, 'newSequenceDB')
    newClusterDB = os.path.join(output_dir, 'newClusterDB')
    repSeqDB = os.path.join(output_dir, 'repSeqDB')
    newRepSeqDB = os.path.join(output_dir, 'newRepSeqDB')
    alignDB = os.path.join(output_dir, 'alignDB')
    newAlignDB = os.path.join(output_dir, 'newAlignDB')

    # Initialize numerical statistics variables
    initial_clusters = 0
    updated_clusters = 0
    new_reps_added = 0
    reps_removed = 0

    # Initial Clustering and Alignment
    logging.info("Step 1: Initial clustering and alignment.")

    logging.info("Creating MMSeqs2 database from old representative sequences.")
    run_command(['mmseqs', 'createdb', args.old, sequenceDB], log_file)

    logging.info("Clustering the old sequences.")
    run_command(['mmseqs', 'cluster', sequenceDB, clusterDB, tmp_dir], log_file)

    logging.info("Creating TSV file of initial clusters.")
    run_command(['mmseqs', 'createtsv', sequenceDB, sequenceDB, clusterDB, os.path.join(output_dir, 'clusterDB.tsv')], log_file)

    logging.info("Aligning the old sequences.")
    run_command(['mmseqs', 'align', sequenceDB, sequenceDB, clusterDB, alignDB, '-a'], log_file)

    logging.info("Converting initial alignments to m8 format.")
    run_command(['mmseqs', 'convertalis', sequenceDB, sequenceDB, alignDB, os.path.join(output_dir, 'align.m8')], log_file)
    logging.info("Step 1 completed.")

    # Concatenate FASTA
    logging.info("Step 2: Concatenating FASTA files.")

    logging.info("Creating MMSeqs2 database from new sequences.")
    run_command(['mmseqs', 'createdb', args.new, addedSequenceDB], log_file)

    logging.info("Concatenating old and new sequence databases.")
    run_command(['mmseqs', 'concatdbs', sequenceDB, addedSequenceDB, allSequenceDB], log_file)

    logging.info("Concatenating header databases.")
    run_command(['mmseqs', 'concatdbs', f"{sequenceDB}_h", f"{addedSequenceDB}_h", f"{allSequenceDB}_h"], log_file)

    logging.info("Converting concatenated database to FASTA.")
    run_command(['mmseqs', 'convert2fasta', allSequenceDB, os.path.join(output_dir, 'allSequenceDB.fasta')], log_file)
    logging.info("Step 2 completed.")

    # Update the Cluster
    logging.info("Step 3: Updating clusters with new sequences.")

    logging.info("Running clusterupdate to incorporate new sequences.")
    run_command(['mmseqs', 'clusterupdate', sequenceDB, allSequenceDB, clusterDB, newSequenceDB, newClusterDB, tmp_dir], log_file)
    logging.info("Step 3 completed.")

    # Count the Number of Updated Clusters
    logging.info("Step 4: Counting clusters.")

    logging.info("Creating TSV file of updated clusters.")
    run_command(['mmseqs', 'createtsv', newSequenceDB, newSequenceDB, newClusterDB, os.path.join(output_dir, 'newClusterDB.tsv')], log_file)

    # Count clusters
    logging.info("Counting initial clusters.")
    initial_clusters = int(subprocess.check_output(f"cut -f1 {os.path.join(output_dir, 'clusterDB.tsv')} | sort | uniq | wc -l", shell=True).decode().strip())

    logging.info("Counting updated clusters.")
    updated_clusters = int(subprocess.check_output(f"cut -f1 {os.path.join(output_dir, 'newClusterDB.tsv')} | sort | uniq | wc -l", shell=True).decode().strip())

    logging.info("Step 4 completed.")

    # Get the Representative Sequences
    logging.info("Step 5: Extracting representative sequences.")

    logging.info("Extracting initial representative sequences.")
    run_command(['mmseqs', 'result2repseq', sequenceDB, clusterDB, repSeqDB], log_file)

    logging.info("Converting initial representative sequences to FASTA.")
    run_command(['mmseqs', 'convert2fasta', repSeqDB, os.path.join(output_dir, 'repSeqDB.fasta')], log_file)

    logging.info("Extracting updated representative sequences.")
    run_command(['mmseqs', 'result2repseq', newSequenceDB, newClusterDB, newRepSeqDB], log_file)

    logging.info("Converting updated representative sequences to FASTA.")
    run_command(['mmseqs', 'convert2fasta', newRepSeqDB, os.path.join(output_dir, 'newRepSeqDB.fasta')], log_file)
    logging.info("Step 5 completed.")

    # See Additions
    logging.info("Step 6: Comparing representative sequences.")

    logging.info("Extracting headers from initial representative sequences.")
    original_reps = os.path.join(output_dir, 'original_reps.txt')
    with open(os.path.join(output_dir, 'repSeqDB.fasta')) as f_in, open(original_reps, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                f_out.write(line.strip() + '\n')

    logging.info("Extracting headers from updated representative sequences.")
    updated_reps = os.path.join(output_dir, 'updated_reps.txt')
    with open(os.path.join(output_dir, 'newRepSeqDB.fasta')) as f_in, open(updated_reps, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                f_out.write(line.strip() + '\n')

    logging.info("Sorting header files for comparison.")
    sorted_original_reps = os.path.join(output_dir, 'original_reps_sorted.txt')
    sorted_updated_reps = os.path.join(output_dir, 'updated_reps_sorted.txt')
    with open(original_reps, 'r') as f_in, open(sorted_original_reps, 'w') as f_out:
        for line in sorted(f_in):
            f_out.write(line)
    with open(updated_reps, 'r') as f_in, open(sorted_updated_reps, 'w') as f_out:
        for line in sorted(f_in):
            f_out.write(line)

    logging.info("Comparing headers to find new and removed representative sequences.")
    # Get counts of unique headers
    new_reps_added = int(subprocess.check_output(f"comm -13 {sorted_original_reps} {sorted_updated_reps} | wc -l", shell=True).decode().strip())
    reps_removed = int(subprocess.check_output(f"comm -23 {sorted_original_reps} {sorted_updated_reps} | wc -l", shell=True).decode().strip())

    logging.info("Step 6 completed.")

    # Get Updated Alignments
    logging.info("Step 7: Generating updated alignments.")

    logging.info("Aligning updated sequences.")
    run_command(['mmseqs', 'align', newSequenceDB, newSequenceDB, newClusterDB, newAlignDB, '-a'], log_file)

    logging.info("Converting updated alignments to m8 format.")
    run_command(['mmseqs', 'convertalis', newSequenceDB, newSequenceDB, newAlignDB, os.path.join(output_dir, 'newAlign.m8')], log_file)
    logging.info("Step 7 completed.")

    # Summary of Numerical Statistics
    logging.info("\nSummary of Numerical Statistics:")
    logging.info(f"- Initial clusters: {initial_clusters}")
    logging.info(f"- Updated clusters: {updated_clusters}")
    logging.info(f"- New representative sequences added: {new_reps_added}")
    logging.info(f"- Representative sequences removed: {reps_removed}")

    # Clean up temporary files
    logging.info("Cleaning up temporary files.")
    shutil.rmtree(tmp_dir)

    logging.info("MMSeqs2 clustering update pipeline completed successfully.")

if __name__ == '__main__':
    main()

