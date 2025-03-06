## Download and trimming 

We downloaded the reads and the assemblies used in the analysis from the  NCBI database (n = 4,124) on 05/12/2023:

```diff

+ # bash #

# For Linux users:
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt

#For MacOS users:
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt

grep "Escherichia coli" assembly_summary_refseq.txt | grep "Complete Genome" > data_escherichia.txt

```
Then, we used **Entrez Direct** (v20.6) (https://www.ncbi.nlm.nih.gov/books/NBK179288/) to download the genome assemblies, with the following arguments:

```diff

+ # bash #

# For Linux users:
esearch -db assembly -query "biosample_ID" \
    | esummary \
    | xtract -pattern DocumentSummary -element FtpPath_GenBank \
    | while read -r line ;
    do
        fname=$(echo $line | grep -o 'GCF_.*' | sed 's/$/_genomic.fna.gz/') ;
        wget "$line/$fname" ;
    done

# For MacOS users:
esearch -db assembly -query "biosample_ID" \
    | esummary \
    | xtract -pattern DocumentSummary -element FtpPath_GenBank \
    | while read -r line ;
    do
        fname=$(echo $line | grep -o 'GCF_.*' | sed 's/$/_genomic.fna.gz/') ;
        curl -O "$line/$fname" ;
    done

```
And we analyzed the number of contigs of each assembly file to filter out those that don't contain plasmids:

```diff

+ # bash #

mkdir LengthAssembly
for file in *.fna; do
    awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' "$file" > "LengthAssembly/$(basename -- "$file" .fna)"
done

+ # Python3 #

import os
current_directory = os.getcwd()

# Create a dictionary to store modified lines
modified_lines = {}

# Modify files
for filename in os.listdir(current_directory):
    if not filename.startswith('MOD_'):
        output_filename = f"MOD_{filename}"
        file_name = os.path.splitext(filename)[0]
        
        with open(filename, 'r') as file:
            with open(output_filename, 'w') as output_file:
                output_file.write("Biosample,contig,sample,length\n")
                for line in file:
                    if line.startswith(">"):
                        parts = line.strip().split(' ', 1)
                        if len(parts) == 2:
                            sample_info = f"{file_name},{parts[0]},{parts[1].rstrip()}"
                            number_line = next(file).strip()
                            output_file.write(f"{sample_info},{number_line}\n")
                            modified_lines[sample_info] = number_line

# Merge modified lines into 'alllengths.txt'
with open('alllengths.txt', 'w') as merged_file:
    merged_file.write("Biosample,contig,sample,length\n")
    for sample_info, number_line in modified_lines.items():
        merged_file.write(f"{sample_info},{number_line}\n")
```

Then, we downloaded the files containing run information for selected assemblies of NCBI, also through esearch: 

```diff

+ # bash #

esearch -db sra -query Biosample | efetch -format runinfo > Biosample.numbers

```

With the Biosample.numbers files, we filter those runs paired-end with "genomic" as a library source and from the platform Illumina, and use the run IDs to download the reads with fasterq-dump from **SRA-toolkit** package (v2.11.3)(https://hpc.nih.gov/apps/sratoolkit.html):


```diff

+ # Python3 #

import multiprocessing
import subprocess

with open('DefNumbers.txt', 'r') as file:    numbers = [(line.strip()) for line in file if line.strip()]

def execute_fasterq_dump(number):
    command = f"fasterq-dump --split-3 {number}"
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Command executed for number: {number}")
    except subprocess.CalledProcessError as e:
        print(f"Error executing command for number {number}: {e}")

if __name__ == "__main__":
    pool = multiprocessing.Pool(processes=8)    
    pool.map(execute_fasterq_dump, numbers)
    pool.close()
    pool.join()

```

Raw reads were trimmed to a Phred score of 20 using **TrimGalore** (v.0.6.6) (https://github.com/FelixKrueger/TrimGalore), with the following parameters:

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

Then, reads were mapped against their assemblies using **Breseq** (v0.38.1) ([https://github.com/ablab/spades](https://github.com/barricklab/breseq)):


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

