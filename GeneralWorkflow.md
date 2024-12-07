Thirty-six evolved populations of the hypermutator strain carrying the plasmid pTA44, and the twenty-four evolved populations of the wild-type strain were sequenced using Illumina. 

## Trimming and assembly 

First, low-quality ends of raw reads were trimmed to a Phred score of 20 using *TrimGalore* (v.0.6.6) (https://github.com/FelixKrueger/TrimGalore), with the following parameters:

```diff

+ # bash #

#!/bin/bash

# Input and output directories
input_dir="/run/user/1000/gvfs/smb-share:server=deepmind.local,share=paula/Mutation_accumulation/readsma/fastqsinadapters"
output_dir="/run/user/1000/gvfs/smb-share:server=deepmind.local,share=paula/Mutation_accumulation/readsma/cutadapt"

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

Then, reads were assembled into contigs using *SPAdes genome assembler* (v3.13.1) (https://github.com/ablab/spades):

```diff
+ # bash #

#!/bin/bash

# Input and output directories
input_dir="/home/paula/Sequences/bsg-ftp.well.ox.ac.uk/220124_A00711_0476_AHTMGTDSX2"
spades_output_base="/home/paula/Sequences/bsg-ftp.well.ox.ac.uk/220124_A00711_0476_AHTMGTDSX2"

# Loop to run SPAdes for each pair of fastq.gz files
for file_R1 in "$input_dir"/*_1.fastq.gz; do
    # Define the R2 file based on R1
    file_R2="${file_R1/_1.fastq.gz/_2.fastq.gz}"

    # Define the output directory based on the base name of the file
    base_name=$(basename "$file_R1" _1.fastq.gz)
    spades_output="$spades_output_base/$base_name"

    # Create the output directory if it doesn't exist
    mkdir -p "$spades_output"

    # Run SPAdes
    /usr/lib/spades/bin/spades.py -1 "$file_R1" -2 "$file_R2" -o "$spades_output"
done

```

And the references were annotated using *Prokka* (v1.14.5) (https://github.com/tseemann/prokka):

```diff
+ # bash #




```
