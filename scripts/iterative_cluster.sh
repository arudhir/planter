#!/bin/bash

# mmseqs_batch_update.sh
# Usage: ./mmseqs_batch_update.sh -g "*.pep" -o base_output_dir -s path_to_mmseqs_cluster_update.py -i initial_old_reps.fasta

# Parse command-line arguments
while getopts g:o:s:i:h option
do
case "${option}"
in
g) GLOB_PATTERN=${OPTARG};; # Glob pattern for input files
o) BASE_OUTPUT_DIR=${OPTARG};; # Base output directory
s) SCRIPT_PATH=${OPTARG};; # Path to mmseqs_cluster_update.py
i) INITIAL_OLD_REPS=${OPTARG};; # Initial old representative sequences file
h) echo "Usage: $0 -g \"*.pep\" -o base_output_dir -s path_to_mmseqs_cluster_update.py -i initial_old_reps.fasta"
   exit;;
*) echo "Invalid option. Usage: $0 -g \"*.pep\" -o base_output_dir -s path_to_mmseqs_cluster_update.py -i initial_old_reps.fasta"
   exit 1;;
esac
done

# Check required arguments
if [ -z "$GLOB_PATTERN" ] || [ -z "$BASE_OUTPUT_DIR" ] || [ -z "$SCRIPT_PATH" ]; then
    echo "Missing arguments. Usage: $0 -g \"*.pep\" -o base_output_dir -s path_to_mmseqs_cluster_update.py -i initial_old_reps.fasta"
    exit 1
fi

# Create the base output directory if it doesn't exist
mkdir -p "$BASE_OUTPUT_DIR"

# Get the list of files matching the glob pattern
FILES=($GLOB_PATTERN)

# Check if files are found
if [ ${#FILES[@]} -eq 0 ]; then
    echo "No files found matching the pattern $GLOB_PATTERN"
    exit 1
fi

# Sort the files (optional)
IFS=$'\n' FILES_SORTED=($(sort <<<"${FILES[*]}"))
unset IFS

# Initialize variables
OLD_REP_SEQS="$INITIAL_OLD_REPS"
COUNTER=1

for FILE in "${FILES_SORTED[@]}"; do
    OUTPUT_DIR="${BASE_OUTPUT_DIR}/output${COUNTER}"
    mkdir -p "$OUTPUT_DIR"

    # If OLD_REP_SEQS is empty, use the current file as initial --old input
    if [ -z "$OLD_REP_SEQS" ]; then
        OLD_REP_SEQS="$FILE"
        echo "No initial representative sequences provided. Using $FILE as initial --old input."
    fi

    echo "======================================="
    echo "Processing file: $FILE"
    echo "Old representative sequences: $OLD_REP_SEQS"
    echo "Output directory: $OUTPUT_DIR"
    echo "======================================="

    # Run the mmseqs_cluster_update.py script
    python3 "$SCRIPT_PATH" --old "$OLD_REP_SEQS" --new "$FILE" --output "$OUTPUT_DIR"

    # Update OLD_REP_SEQS for the next iteration
    OLD_REP_SEQS="${OUTPUT_DIR}/newRepSeqDB.fasta"

    # Increment counter
    ((COUNTER++))
done

