## Download and trimming 

Raw reads were trimmed to a Phred score of 20 using *TrimGalore* (v.0.6.6) (https://github.com/FelixKrueger/TrimGalore), with the following parameters:

```diff

+ # bash #

#!/bin/bash

# Input and output directories
input_dir="/input_dir/rawreads"
output_dir="/input_dir/trim_reads"

# Create the output directory 
mkdir -p "$output_dir"

for file in "$input_dir"/*_R1.fastq.gz; do
  
    file_R2="${file/_R1/_R2}"

    # Define the base name of the file (without extension)
    base_name=$(basename "$file" _R1.fastq.gz)

    # Run Trim Galore with the specified parameters
    trim_galore --quality 20 --illumina --fastqc --gzip --paired "$file" "$file_R2" --output_dir "$output_dir"
done

```

## Variant calling 

Then, reads were mapped against their assemblies using *Breseq* (v0.38.1) ([https://github.com/ablab/spades](https://github.com/barricklab/breseq)):


```diff

# Exit immediately if a command exits with a non-zero status,
# Treat unset variables as an error, and
# Prevent errors in a pipeline from being masked

set -euo pipefail
IFS=$'\n\t'

# Enable nullglob to ensure non-matching globs result in empty arrays
shopt -s nullglob

# Function to check if a command exists
check_command() {
    command -v "$1" >/dev/null 2>&1 || { echo >&2 "Error: $1 is not installed or not in PATH. Aborting."; exit 1; }
}

# Check for required tools
check_command breseq
check_command parallel

# Determine the base directory where the script is located
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Input file
INPUT_FILE="/storage/MA/B_paper_code/1_files/Assembly_runs.txt"

# Directories (Ensure correct case)
ASSEMBLY_DIR="/storage/MA/2_assemblies/Assemblies"
READS_DIR="/storage/MA/1_reads/Batch3_jero"
PAIRED_DIR="/storage/MA/1_reads/Batch4_paula"
SINGLE_DIR="$READS_DIR/Single/cleaned"
OUTPUT_DIR="$BASE_DIR/Results_breseq_trimme_batch3jero"

# Create the output directory if it does not exist
mkdir -p "$OUTPUT_DIR"

process_entry() {
    ASSEMBLY="$1"
    RUN="$2"

    echo "----------------------------------------"
    echo "Processing Assembly: '$ASSEMBLY', Run: '$RUN'"

    # Define a specific output directory for this run
    SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${RUN}"

    # Check if the output directory already exists
    if [[ -d "$SAMPLE_OUTPUT_DIR" ]]; then
        echo "Output directory already exists for Run: $RUN. Skipping..."
        return
    fi

    # Create the output directory
    mkdir -p "$SAMPLE_OUTPUT_DIR"

    # Search for the reference file corresponding to the Assembly
    REF_FILES=("${ASSEMBLY_DIR}/${ASSEMBLY}"_*.fna)

    # Debug: Print the REF_FILES array
    if [[ ${#REF_FILES[@]} -eq 0 ]]; then
        echo "Files found: None"
    else
        echo "Files found: ${REF_FILES[@]}"
    fi

    if [[ ${#REF_FILES[@]} -eq 0 ]]; then
        echo "No reference file found for Assembly: $ASSEMBLY. Skipping..."
        return
    elif [[ ${#REF_FILES[@]} -gt 1 ]]; then
        echo "Multiple reference files found for Assembly: $ASSEMBLY. Using the first one found: ${REF_FILES[0]}"
    fi

    REF_FILE="${REF_FILES[0]}"

    if [[ ! -f "$REF_FILE" ]]; then
        echo "Invalid reference file: $REF_FILE. Skipping..."
        return
    fi

    # Check if the files are Paired or Single
    PAIR1="$PAIRED_DIR/${RUN}_1_val_1.fq.gz"  
    PAIR2="$PAIRED_DIR/${RUN}_2_val_2.fq.gz"
    SINGLE="$SINGLE_DIR/${RUN}_trimmed.fq.gz"

    if [[ -f "$PAIR1" && -f "$PAIR2" ]]; then
        READ_TYPE="paired"
        READ1="$PAIR1"
        READ2="$PAIR2"
        echo "Paired-end reads found: $READ1 and $READ2."
    elif [[ -f "$SINGLE" ]]; then
        READ_TYPE="single"
        READ1="$SINGLE"
        READ2=""
        echo "Single-end reads found: $READ1."
    else
        echo "No read files found for Run: $RUN. Skipping..."
        return
    fi

    # Run breseq
    echo "Running breseq for Run: $RUN..."

    if [[ "$READ_TYPE" == "paired" ]]; then
        breseq -r "$REF_FILE" -o "$SAMPLE_OUTPUT_DIR" -p "$READ1" "$READ2"
    else
        breseq -r "$REF_FILE" -o "$SAMPLE_OUTPUT_DIR" -p "$READ1"
    fi

    if [[ $? -ne 0 ]]; then
        echo "Error running breseq for Run: $RUN. Skipping..."
        return
    fi

    echo "Processing completed for Run: $RUN. Results in: $SAMPLE_OUTPUT_DIR"
}

export -f process_entry
export ASSEMBLY_DIR READS_DIR PAIRED_DIR SINGLE_DIR OUTPUT_DIR

# Determine the number of CPU cores
TOTAL_CORES=$(nproc)
THREADS_PER_JOB=8
MAX_JOBS=$((TOTAL_CORES / THREADS_PER_JOB))

# Ensure at least one job is allowed
if [[ "$MAX_JOBS" -lt 1 ]]; then
    MAX_JOBS=1
fi

echo "Starting pipeline with a maximum of $MAX_JOBS jobs in parallel."

# Read the input file and pass each line to GNU parallel with a limited number of jobs
tail -n +2 "$INPUT_FILE" | parallel --colsep '\t' -j "$MAX_JOBS" process_entry {1} {2}

echo "Pipeline completed."

```
Finally, mutations were analyzing usin R v4.1.2

