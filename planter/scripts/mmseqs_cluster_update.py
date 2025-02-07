#!/usr/bin/env python3

import argparse
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple

# Set up logging at the module level before any classes
def setup_logging():
    """Configure root logger to output to both stdout and file."""
    # Clear any existing handlers
    root_logger = logging.getLogger()
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Set the logging level
    root_logger.setLevel(logging.INFO)
    
    # Create console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_formatter)
    root_logger.addHandler(console_handler)

# Initialize logging
setup_logging()

class MMseqsClusterUpdater:
    """Class to handle MMseqs2 cluster updates for a single pair of files."""
    
    def __init__(self, output_dir: str, tmp_dir: Optional[str] = None):
        self.output_dir = output_dir
        self.tmp_dir = tmp_dir or os.path.join(output_dir, 'tmp')
        self._setup_directories()
        self._setup_logging()
        self._initialize_paths()
        
    def _setup_directories(self) -> None:
        """Create output and temporary directories."""
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.tmp_dir, exist_ok=True)
        
    def _setup_logging(self) -> None:
        """Configure logging to both file and stdout."""
        log_file = os.path.join(self.output_dir, 'pipeline.log')
        logging.basicConfig(
            level=logging.INFO,
            format='%(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        
    def _initialize_paths(self) -> None:
        """Initialize paths for databases and files."""
        self.paths = {
            'sequenceDB': os.path.join(self.output_dir, 'sequenceDB'),
            'clusterDB': os.path.join(self.output_dir, 'clusterDB'),
            'addedSequenceDB': os.path.join(self.output_dir, 'addedSequenceDB'),
            'allSequenceDB': os.path.join(self.output_dir, 'allSequenceDB'),
            'newSequenceDB': os.path.join(self.output_dir, 'newSequenceDB'),
            'newClusterDB': os.path.join(self.output_dir, 'newClusterDB'),
            'repSeqDB': os.path.join(self.output_dir, 'repSeqDB'),
            'newRepSeqDB': os.path.join(self.output_dir, 'newRepSeqDB'),
            'alignDB': os.path.join(self.output_dir, 'alignDB'),
            'newAlignDB': os.path.join(self.output_dir, 'newAlignDB')
        }
        
    def run_command(self, cmd: List[str]) -> None:
        """Run a shell command and log its output."""
        log_file = os.path.join(self.output_dir, 'pipeline.log')
        logging.info(f"Running command: {' '.join(cmd)}")
        with open(log_file, 'a') as lf:
            lf.write(f"Running command: {' '.join(cmd)}\n")
            process = subprocess.Popen(cmd, stdout=lf, stderr=lf, shell=isinstance(cmd, str))
            process.communicate()
            if process.returncode != 0:
                sys.exit(f"Command failed: {' '.join(cmd)}")

    def update_clusters(self, old_seqs: str, new_seqs: str, output_new: str, output_removed: str) -> Tuple[int, int, int, int]:
        """
        Run the complete MMseqs2 clustering update pipeline.
        Returns tuple of (initial_clusters, updated_clusters, new_reps_added, reps_removed)
        """
        logging.info("Starting MMSeqs2 clustering update pipeline.")
        
        # Step 1: Initial clustering and alignment
        self._run_initial_clustering(old_seqs)
        
        # Step 2: Add new sequences
        self._add_new_sequences(new_seqs)
        
        # Step 3: Update clusters
        self._update_clusters()
        
        # Step 4: Count clusters
        initial_clusters, updated_clusters = self._count_clusters()
        
        # Step 5: Get representative sequences
        self._extract_representative_sequences()
        
        # Step 6: Compare representatives
        new_reps_added, reps_removed = self._compare_representatives(output_new, output_removed)
        
        # Step 7: Generate updated alignments
        self._generate_alignments()
        
        # Clean up
        self._cleanup()
        
        return initial_clusters, updated_clusters, new_reps_added, reps_removed

    def _run_initial_clustering(self, old_seqs: str) -> None:
        """Run initial clustering steps."""
        logging.info("Step 1: Initial clustering and alignment.")
        self.run_command(['mmseqs', 'createdb', old_seqs, self.paths['sequenceDB']])
        self.run_command(['mmseqs', 'cluster', self.paths['sequenceDB'], self.paths['clusterDB'], self.tmp_dir])
        self.run_command(['mmseqs', 'createtsv', self.paths['sequenceDB'], self.paths['sequenceDB'], 
                         self.paths['clusterDB'], os.path.join(self.output_dir, 'clusterDB.tsv')])
        self.run_command(['mmseqs', 'align', self.paths['sequenceDB'], self.paths['sequenceDB'], 
                         self.paths['clusterDB'], self.paths['alignDB'], '-a'])
        self.run_command(['mmseqs', 'convertalis', self.paths['sequenceDB'], self.paths['sequenceDB'], 
                         self.paths['alignDB'], os.path.join(self.output_dir, 'align.m8')])

    def _add_new_sequences(self, new_seqs: str) -> None:
        """Add new sequences to the database."""
        logging.info("Step 2: Adding new sequences.")
        self.run_command(['mmseqs', 'createdb', new_seqs, self.paths['addedSequenceDB']])
        self.run_command(['mmseqs', 'concatdbs', self.paths['sequenceDB'], self.paths['addedSequenceDB'], 
                         self.paths['allSequenceDB']])
        self.run_command(['mmseqs', 'concatdbs', f"{self.paths['sequenceDB']}_h", 
                         f"{self.paths['addedSequenceDB']}_h", f"{self.paths['allSequenceDB']}_h"])

    def _update_clusters(self) -> None:
        """Update clusters with new sequences."""
        logging.info("Step 3: Updating clusters.")
        self.run_command(['mmseqs', 'clusterupdate', self.paths['sequenceDB'], self.paths['allSequenceDB'],
                         self.paths['clusterDB'], self.paths['newSequenceDB'], self.paths['newClusterDB'], 
                         self.tmp_dir])

    def _count_clusters(self) -> Tuple[int, int]:
        """Count initial and updated clusters."""
        logging.info("Step 4: Counting clusters.")
        self.run_command(['mmseqs', 'createtsv', self.paths['newSequenceDB'], self.paths['newSequenceDB'],
                         self.paths['newClusterDB'], os.path.join(self.output_dir, 'newClusterDB.tsv')])
        
        initial = int(subprocess.check_output(
            f"cut -f1 {os.path.join(self.output_dir, 'clusterDB.tsv')} | sort | uniq | wc -l", 
            shell=True).decode().strip())
        updated = int(subprocess.check_output(
            f"cut -f1 {os.path.join(self.output_dir, 'newClusterDB.tsv')} | sort | uniq | wc -l", 
            shell=True).decode().strip())
        return initial, updated

    def _extract_representative_sequences(self) -> None:
        """Extract representative sequences from clusters."""
        logging.info("Step 5: Extracting representative sequences.")
        
        # Extract initial representatives
        logging.info("Extracting initial representative sequences.")
        self.run_command([
            'mmseqs', 
            'result2repseq',
            self.paths['sequenceDB'],
            self.paths['clusterDB'],
            self.paths['repSeqDB']
        ])
        self.run_command([
            'mmseqs',
            'convert2fasta',
            self.paths['repSeqDB'],
            os.path.join(self.output_dir, 'repSeqDB.fasta')
        ])
        
        # Extract updated representatives
        logging.info("Extracting updated representative sequences.")
        self.run_command([
            'mmseqs',
            'result2repseq',
            self.paths['newSequenceDB'],
            self.paths['newClusterDB'],
            self.paths['newRepSeqDB']
        ])
        self.run_command([
            'mmseqs',
            'convert2fasta',
            self.paths['newRepSeqDB'],
            os.path.join(self.output_dir, 'newRepSeqDB.fasta')
        ])

    def _compare_representatives(self, output_new: str, output_removed: str) -> Tuple[int, int]:
        """Compare initial and updated representative sequences using seqkit and save outputs."""
        logging.info("Step 6: Comparing representatives using seqkit.")

        original_reps = os.path.join(self.output_dir, 'original_reps.txt')
        updated_reps = os.path.join(self.output_dir, 'updated_reps.txt')

        # Ensure the expected directory exists
        os.makedirs(os.path.dirname(output_new), exist_ok=True)
        os.makedirs(os.path.dirname(output_removed), exist_ok=True)

        # Extract headers using seqkit
        self.run_command(['seqkit', 'seq', '-n', os.path.join(self.output_dir, 'repSeqDB.fasta'), '-o', original_reps])
        self.run_command(['seqkit', 'seq', '-n', os.path.join(self.output_dir, 'newRepSeqDB.fasta'), '-o', updated_reps])

        try:
            # Save new sequences
            self.run_command(f"bash -c 'comm -13 <(sort {original_reps}) <(sort {updated_reps}) > {output_new}'")
            # Save removed sequences
            self.run_command(f"bash -c 'comm -23 <(sort {original_reps}) <(sort {updated_reps}) > {output_removed}'")

            new_added = sum(1 for _ in open(output_new))
            removed = sum(1 for _ in open(output_removed))

            logging.info(f"New representative sequences added: {new_added}")
            logging.info(f"Representative sequences removed: {removed}")

            return new_added, removed
        except subprocess.CalledProcessError as e:
            logging.error(f"Error comparing representative sequences: {e.output.decode()}")
            return 0, 0

    def _generate_alignments(self) -> None:
        """Generate alignments for updated clusters."""
        logging.info("Step 7: Generating alignments.")
        self.run_command(['mmseqs', 'align', self.paths['newSequenceDB'], self.paths['newSequenceDB'],
                         self.paths['newClusterDB'], self.paths['newAlignDB'], '-a'])
        self.run_command(['mmseqs', 'convertalis', self.paths['newSequenceDB'], self.paths['newSequenceDB'],
                         self.paths['newAlignDB'], os.path.join(self.output_dir, 'newAlign.m8')])

    def _cleanup(self) -> None:
        """Clean up temporary files."""
        logging.info("Cleaning up temporary files.")
        shutil.rmtree(self.tmp_dir)
        logging.info("MMSeqs2 clustering update pipeline completed successfully.")


class BatchMMseqsUpdater:
    """Class to handle batch processing of multiple sequence files."""
    
    def __init__(self, files: List[str], base_output_dir: str, initial_reps: Optional[str] = None):
        self.files = files
        self.base_output_dir = base_output_dir
        self.initial_reps = initial_reps
        os.makedirs(base_output_dir, exist_ok=True)
    
    def _get_srr_id(self, filepath: str) -> str:
        """Extract SRR ID from filepath."""
        filename = os.path.basename(filepath)
        # Remove .pep extension
        srr_id = filename.replace('.pep', '')
        return srr_id
        
    def process_files(self) -> None:
        """Process all files in the list."""
        if not self.files:
            sys.exit("No input files provided")
            
        files = self._validate_files()
        old_rep_seqs = self.initial_reps
        
        if not old_rep_seqs and files:
            old_rep_seqs = files[0]
            first_srr = self._get_srr_id(files[0])
            logging.info(f"No initial representative sequences provided. Using {first_srr} ({files[0]}) as initial input.")
            files = files[1:]
        
        for i, file in enumerate(files, start=1):
            srr_id = self._get_srr_id(file)
            output_dir = os.path.join(self.base_output_dir, f"update_{srr_id}")
            logging.info(f"\nProcessing file {i}/{len(files)}: {srr_id}")
            
            updater = MMseqsClusterUpdater(output_dir)
            stats = updater.update_clusters(old_rep_seqs, file)
            
            self._log_iteration_stats(i, srr_id, stats)
            old_rep_seqs = os.path.join(output_dir, "newRepSeqDB.fasta")
    
    def _validate_files(self) -> List[str]:
        """Validate that all files exist and follow the pattern."""
        for file in self.files:
            if not os.path.exists(file):
                sys.exit(f"File not found: {file}")
            if not file.endswith('.pep'):
                sys.exit(f"File {file} does not end with .pep extension")
            filename = os.path.basename(file)
            if not filename.replace('.pep', ''):
                sys.exit(f"Could not extract SRR ID from filename: {filename}")
        return self.files
    
    def _log_iteration_stats(self, iteration: int, srr_id: str, stats: Tuple[int, int, int, int]) -> None:
        """Log statistics for each iteration."""
        initial, updated, added, removed = stats
        logging.info(f"\nIteration {iteration} Statistics for {srr_id}:")
        logging.info(f"- Initial clusters: {initial}")
        logging.info(f"- Updated clusters: {updated}")
        logging.info(f"- New representative sequences added: {added}")
        logging.info(f"- Representative sequences removed: {removed}")


def main():
    parser = argparse.ArgumentParser(description="MMseqs2 clustering update tool")
    parser.add_argument('-i', '--input-old', required=True, help='Old representative sequences file')
    parser.add_argument('-p', '--input-new', required=True, help='New protein sequences file')
    parser.add_argument('-o', '--output-dir', required=True, help='Output directory')
    parser.add_argument('--output-new', required=True, help='Output file for newly added sequences')
    parser.add_argument('--output-removed', required=True, help='Output file for removed sequences')

    args = parser.parse_args()
    
    updater = MMseqsClusterUpdater(args.output_dir)
    updater.update_clusters(args.input_old, args.input_new, args.output_new, args.output_removed)

if __name__ == '__main__':
    main()