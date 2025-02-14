![image](https://github.com/user-attachments/assets/debd788e-8821-493d-bc1c-5d0c59e35de5)


## Hipermutator strains analysis

We download a reference assembly for Escherichia coli (ASM584v2, 24/01/2025) and extract the genes involved in the mutator ability: 

```diff
+ # bash #

grep -iE -A 1 "mutS|mutL|mutH|mutT|mutY|mutM|dnaQ|uvrA|uvrB|uvrC" mg1655_bakta/mg1655.ffn > mutator_genes.fasta

# we name the contigs only with the name of the gene

sed -i '/--/d' mutator_genes.fasta
sed -i 's/^>\(.*\) \(.*\)$/>\2/' mutator_genes.fasta

```

We annotate the .fasta file with the genes of interest with *Prokka v.1.14.6* (https://github.com/tseemann/prokka):

```diff
+ # bash #
prokka mutator_genes.fasta  --prefix mutator_genes --out mutator_genes --force

```
And perform a variant calling with *snippy v4.6.0* (https://github.com/tseemann/snippy):

```diff
+ # bash #
for file in *.fna; do snippy --outdir ${file%_genomic.fna} --ref /storage/MA/6_mutator_analysis/mutator_genes/mutator_genes.gbk --ctgs ${file} --force; done

# We rename the .tab with the results as the name of the sample, and add a column with the name of the strain

for f in */*.tab ;do fp=$(dirname "$f"); fl=$(basename "$f"); mv "$fp/$fl" "$fp/$fp"_"$fl"; done

for f in *.tab; do awk -v fName="${f%.tab}" '{printf("%s,%s\n", (FNR==1 ? "filename" : fName), $0)}' "$f" > mod"$f"; done

 cat mod* > snippy_results_mutators.txt

```

Extract information from .vcf results:

```diff
+ # bash #

for file in */*_snps.vcf; do grep -E '^#CHROM|^[^#]' ${file} > ${file%_snps.vcf}_limpio.vcf; done



```

Download annotations: 

```diff
+ # bash #

#!/bin/bash

# Archivo con la lista de ensamblajes
assembly_file="assembly_list.txt"

# Generar las URLs desde el archivo de ensamblajes
cat "$assembly_file" | while read -r assembly; do
    dir1=$(echo $assembly | cut -c 5-7)
    dir2=$(echo $assembly | cut -c 8-10)
    dir3=$(echo $assembly | cut -c 11-13)
    echo "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${dir1}/${dir2}/${dir3}/${assembly}/${assembly}_protein.faa.gz"
done > urls_to_download.txt

# Descargar las URLs en paralelo usando `parallel`
parallel -j 8 wget {} < urls_to_download.txt
```




# Lenski data

Runs:
SRR2584408
SRR2584409
SRR2584410
SRR2584411
SRR2584438
SRR2584462
SRR2584463
SRR2584464
SRR2584466
SRR2584467
SRR2584470
SRR2584490
SRR2584518
SRR2584508
SRR2584519
SRR2584520
SRR2584521
SRR2584522
SRR2086043
SRR2086070
SRR2584524
SRR2584534
SRR2584549
SRR2584550
SRR2584583
SRR2584607
SRR2584608
SRR2584609
SRR2584610
SRR2584611
SRR2584613
SRR2584614
SRR2584637
SRR2584654
SRR2584655
SRR2584659
SRR2584660
SRR2584661
SRR2584662
SRR2584663
SRR2086096
SRR2086097
SRR2584665
SRR2588645
SRR2588658
SRR2584668
SRR2584669
SRR2591034
SRR2584671
SRR2591035
SRR2584673
SRR2588711
SRR2584677
SRR2584678
SRR2591036
SRR2591037
SRR2584681
SRR2584682
SRR2584683
SRR2584684
SRR2584685
SRR2588848
SRR1536189
SRR2584694
SRR2584696
SRR2591038
SRR2584697
SRR2584698
SRR2588900
SRR2588923
SRR2588946
SRR2584702
SRR2584703
SRR2584705
SRR2584706
SRR2588991
SRR2584708
SRR2588992
SRR2584710
SRR2588993
SRR2584712
SRR2591039
SRR2584714
SRR2588995
SRR2086098
SRR2086148
SRR2584774
SRR2591040
SRR2584776
SRR2584777
SRR2584778
SRR2584779
SRR2584780
SRR2584781
SRR2584782
SRR2584784
SRR2584785
SRR2584786
SRR2584787
SRR2584811
SRR2584812
SRR2584813
SRR2584814
SRR2584815
SRR2584816
SRR2584820
SRR2086149
SRR2086150
SRR2584821
SRR2584822
SRR2584823
SRR2584824
SRR2584825
SRR2584826
SRR2584827
SRR2584828
SRR2584829
SRR2584831
SRR2584832
SRR2584833
SRR2584834
SRR2584835
SRR2584836
SRR2584837
SRR2584838
SRR2584839
SRR2584840
SRR2584842
SRR1536190
SRR2584843
SRR2584844
SRR2588999
SRR2591031
SRR2588223
SRR2584817
SRR2584854
SRR2591042
SRR2584876
SRR030253
ERR051713
SRR030254
ERR051724
SRR030255
ERR051718
SRR030256
ERR051719
SRR030257
ERR051714
ERR051722
ERR051721
SRR030258
ERR051717
SRR2584887
SRR2591047
SRR2591050
SRR2589073
SRR2584405
SRR2584465
SRR2591033
SRR2584612
SRR2584664
SRR2584676
SRR2584693
SRR2584704
SRR2584715
SRR2584783
SRR2584819
SRR2584830
SRR2584841
SRR2584846
SRR2584847
SRR2584848
SRR1536187
SRR2086152
SRR2584849
SRR2584851
SRR2584852
SRR2584853
SRR2589001
SRR2584856
SRR2584857
SRR2584858
SRR2584859
SRR2591041
SRR2589044
SRR098281
SRR2589045
SRR098282
SRR098283
SRR2584863
SRR098284
SRR098285
SRR098031
SRR098280
SRR098029
SRR098030
SRR2584864
SRR2584866
SRR2584867
SRR2584868
SRR2584869
SRR2584870
SRR2584871
SRR2591043
SRR2584873
SRR2584874
SRR2591044
SRR2589049
SRR2584878
SRR2589050
SRR2591045
SRR2584881
SRR2584882
SRR2584883
SRR2584884
SRR2584885
SRR1536188
SRR2584886
SRR2584888
SRR2584889
SRR2589054
SRR2589055
SRR2591046
SRR2584893
SRR2584894
SRR2589057
SRR2589058
SRR2589059
SRR2584899
SRR2584900
SRR2589061
SRR2589062
SRR2591048
SRR2591049
SRR2584905
SRR2589064
SRR2589065
SRR2584908
SRR2085290
SRR2086039
SRR2584927
SRR2589067
SRR2584929
SRR2589068
SRR2584931
SRR2584932
SRR2591051
SRR2584934
SRR2589071
SRR2584936
SRR2584938
SRR2591052
SRR2591053
SRR2591054
SRR2584942
SRR2584943
SRR2589074
SRR2584978
SRR2584979
SRR2589077
SRR2086041
SRR2086042
SRR2584407
SRR2584406 


We download reads and assemblies from hypermutators and non-hypermutators strains with fasterq-dump:


```diff
+ # python #
import multiprocessing
import subprocess
import time

def execute_fasterq_dump(number):
    command = f"fasterq-dump --split-3 {number}"
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Command executed successfully for: {number}")
    except subprocess.CalledProcessError as e:
        print(f"Error executing command for {number}: {e}")
        # Opción de reintentar, ajusta el número de intentos como sea necesario
        for i in range(3):  # Intenta hasta 3 veces si falla
            print(f"Retrying for {number}, attempt {i + 1}")
            try:
                subprocess.run(command, shell=True, check=True)
                print(f"Command executed successfully on retry for: {number}")
                break
            except subprocess.CalledProcessError as e:
                print(f"Error retrying for {number}: {e}")
                if i == 2:  # Después de 3 intentos fallidos, no lo volverá a intentar
                    print(f"Failed to execute for {number} after 3 attempts")

if __name__ == "__main__":
    with open('runs.txt', 'r') as file:
        numbers = [(line.strip()) for line in file if line.strip()]

    # Usando el 'with' para un manejo adecuado de recursos
    with multiprocessing.Pool(processes=8) as pool:  # Define el número de procesos según tu máquina
        pool.map(execute_fasterq_dump, numbers)

```

And we trim the reads using *Trim-galore v0.6.10* (https://github.com/FelixKrueger/TrimGalore):

```diff
+ # bash #
#!/bin/bash

# Input and output directories
input_dir="/storage/MA/6_mutator_analysis/Lenski/single"
output_dir="/storage/MA/6_mutator_analysis/Lenski/single/trimreads"

# Create the output directory
mkdir -p "$output_dir"

# Define the trim_galore function to be used by parallel
trim_galore_func() {
    local file="$1"
    local base_name=$(basename "$file" .fastq)

    trim_galore --quality 20 --illumina --fastqc --gzip "$file" --output_dir "$output_dir"
}

export -f trim_galore_func

# Get the number of CPU cores available
num_cores=$(nproc)

# Find all .fastq files and run the trim_galore_func in parallel, using the number of cores
find "$input_dir"/*.fastq | parallel -j "$num_cores" trim_galore_func

```


Then, we assemble the reads with *SPAdes v3.13.1* (https://github.com/ablab/spades): 

```diff
+ # bash #

# For paired-end reads:

#!/bin/bash

# Input and output directories
input_dir="/storage/MA/6_mutator_analysis/Lenski/paired"
spades_output_base="/storage/MA/6_mutator_analysis/Lenski/paired/spades"

# Loop to run SPAdes for each pair of fastq.gz files
for file_R1 in "$input_dir"/*_1_val_1.fq.gz; do
    # Define the R2 file based on R1
    file_R2="${file_R1/_1_val_1.fq.gz/_2_val_2.fq.gz}"

    # Check if the R2 file exists
    if [[ -f "$file_R2" ]]; then
        # Define the output directory based on the base name of the file
        base_name=$(basename "$file_R1" _1_val_1.fq.gz)
        spades_output="$spades_output_base/$base_name"

        # Create the output directory if it doesn't exist
        mkdir -p "$spades_output"

        # Run SPAdes
        /usr/lib/spades/bin/spades.py -1 "$file_R1" -2 "$file_R2" -o "$spades_output"
    else
        echo "R2 file for $file_R1 not found. Skipping."
    fi
done

# For single reads:

#!/bin/bash

# Directorios de entrada y salida
input_dir="/storage/MA/6_mutator_analysis/Lenski/single"
spades_output_base="/storage/MA/6_mutator_analysis/Lenski/single/spades"

# Loop para ejecutar SPAdes para cada archivo de read single
for file_R1 in "$input_dir"/*_trimmed.fq.gz; do
    # Definir el nombre base del archivo
    base_name=$(basename "$file_R1" _trimmed.fq.gz)
    spades_output="$spades_output_base/$base_name"

    # Crear el directorio de salida si no existe
    mkdir -p "$spades_output"

    # Ejecutar SPAdes
    /usr/lib/spades/bin/spades.py -s "$file_R1" -o "$spades_output"
done


```

And parsed the outputs with the following code:

```diff

+ # python #

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VCF Parser Script

This script parses all .vcf files located within 'output' folders inside each SRRXXXX directory
within a master directory. It extracts relevant VCF information and aggregates it into a
summary CSV or Excel file, including the sample name derived from the SRRXXXX folder.

Usage:
    python vcf_parser_corrected.py -m /path/to/master_directory -o summary.xlsx

Author: ChatGPT
Date: 2024-04-27
"""

import os
import argparse
import pandas as pd
import glob

def parse_vcf(vcf_file, sample_name):
    """
    Parses a single VCF file and extracts relevant information.

    Parameters:
        vcf_file (str): Path to the VCF file.
        sample_name (str): Name of the sample (SRRXXXX).

    Returns:
        pd.DataFrame: DataFrame containing parsed VCF data with an added 'Sample' column.
    """
    try:
        with open(vcf_file, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error reading {vcf_file}: {e}")
        return pd.DataFrame()  # Return empty DataFrame on error

    data = []
    headers = []
    for line in lines:
        line = line.strip()
        if line.startswith('##'):
            continue  # Skip meta-information lines
        elif line.startswith('#CHROM'):
            headers = line.lstrip('#').split('\t')
        else:
            if not headers:
                print(f"Warning: No header found before data lines in {vcf_file}. Skipping file.")
                return pd.DataFrame()
            parts = line.split('\t')
            if len(parts) < 8:
                print(f"Warning: Incomplete VCF entry in {vcf_file}: {line}. Skipping this line.")
                continue
            record = dict(zip(headers, parts[:8]))
            info_fields = record['INFO'].split(';')
            info_dict = {}
            for field in info_fields:
                if '=' in field:
                    key, value = field.split('=', 1)
                    info_dict[key] = value
                else:
                    info_dict[field] = True  # Flag fields

            # Extract AF, AD, DP if available
            AF = info_dict.get('AF', 'NA')
            AD = info_dict.get('AD', 'NA')
            DP = info_dict.get('DP', 'NA')

            # Append extracted information
            data.append({
                'Sample': sample_name,
                'CHROM': record.get('CHROM', 'NA'),
                'POS': record.get('POS', 'NA'),
                'ID': record.get('ID', 'NA'),
                'REF': record.get('REF', 'NA'),
                'ALT': record.get('ALT', 'NA'),
                'QUAL': record.get('QUAL', 'NA'),
                'FILTER': record.get('FILTER', 'NA'),
                'AF': AF,
                'AD': AD,
                'DP': DP
            })

    if not data:
        print(f"No data parsed from {vcf_file}.")
        return pd.DataFrame()

    df = pd.DataFrame(data)
    return df

def find_vcf_files_corrected(master_dir):
    """
    Finds all .vcf files within 'output' folders inside SRRXXXX directories.

    Parameters:
        master_dir (str): Path to the master directory.

    Returns:
        list of tuples: Each tuple contains (VCF file path, Sample name).
    """
    vcf_files = []
    # List all directories in master_dir that start with 'SRR'
    sr_dirs = [d for d in os.listdir(master_dir) if os.path.isdir(os.path.join(master_dir, d)) and d.startswith('SRR')]
    for sr_dir in sr_dirs:
        sr_path = os.path.join(master_dir, sr_dir)
        output_dir = os.path.join(sr_path, 'output')
        vcf_path = os.path.join(output_dir, 'output.vcf')
        if os.path.isfile(vcf_path):
            vcf_files.append((vcf_path, sr_dir))
        else:
            print(f"Warning: No VCF file found in {output_dir}. Expected at least 'output.vcf'. Skipping {sr_dir}.")
    return vcf_files

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Parse multiple VCF files from a master directory and aggregate the data.'
    )
    parser.add_argument('-m', '--master_dir', required=True, help='Path to the master directory containing SRRXXXX folders.')
    parser.add_argument('-o', '--output_file', required=True, help='Path for the output summary file (e.g., summary.xlsx or summary.csv).')
    args = parser.parse_args()

    master_dir = args.master_dir
    output_file = args.output_file

    if not os.path.isdir(master_dir):
        print(f"Error: The specified master directory does not exist: {master_dir}")
        return

    # Find all VCF files
    vcf_files = find_vcf_files_corrected(master_dir)
    if not vcf_files:
        print("No VCF files found. Exiting.")
        return

    print(f"Found {len(vcf_files)} VCF files. Starting parsing...")

    # Parse all VCF files and collect data
    all_data = []
    for vcf_path, sample in vcf_files:
        print(f"Parsing {vcf_path} for sample {sample}...")
        df = parse_vcf(vcf_path, sample)
        if not df.empty:
            all_data.append(df)

    if not all_data:
        print("No data parsed from any VCF files. Exiting.")
        return

    # Concatenate all DataFrames
    summary_df = pd.concat(all_data, ignore_index=True)

    # Save to Excel or CSV based on the output file extension
    try:
        if output_file.lower().endswith('.xlsx') or output_file.lower().endswith('.xls'):
            summary_df.to_excel(output_file, index=False)
            print(f"Successfully saved summary to {output_file}")
        elif output_file.lower().endswith('.csv'):
            summary_df.to_csv(output_file, index=False)
            print(f"Successfully saved summary to {output_file}")
        else:
            print("Error: Output file must have a .xlsx, .xls, or .csv extension.")
    except Exception as e:
        print(f"Error saving the summary file: {e}")

if __name__ == "__main__":
    main()
    
```

# Genomic information from paper:

https://www.nature.com/articles/nature24287

We downloaded the reads corresponding to samples from 50,000 to 60,000 generations.

```diff

+ # bash #

cat Runs.txt | parallel -j 8 "fasterq-dump {} --split-files --progress"

```

And we trim the reads using *Trim-galore v0.6.10* (https://github.com/FelixKrueger/TrimGalore):


```diff

+ # bash #

# Input and output directories
input_dir="/mnt/9c115607-eaba-44b4-bf47-a7c423299bfb/Mutation_accumulation/Good_2017/reads/paired"
output_dir="/mnt/9c115607-eaba-44b4-bf47-a7c423299bfb/Mutation_accumulation/Good_2017/reads/paired/trimreads"

# Create the output directory
mkdir -p "$output_dir"

# Define the trim_galore function to be used by parallel
trim_galore_func() {
    local file1="$1"
    local file2="${file1/_1.fastq/_2.fastq}"
    local base_name=$(basename "$file1" _1.fastq)

    trim_galore --quality 20 --illumina --fastqc --gzip --paired "$file1" "$file2" --output_dir "$output_dir"
}

export -f trim_galore_func

# Get the number of CPU cores available
# num_cores=$(nproc)

# Find all _1.fastq files and run the trim_galore_func in parallel
find "$input_dir" -name "*_1.fastq" | parallel -j 16 trim_galore_func

```

Then, we assembled the trimmed reads with spades: 

```diff

+ # bash #

cat list_SRR.txt | parallel -j 16 "spades.py -1 ${}1_val_1.fq.gz -2 ${}2_val_2.fq.gz -o ${}"

```

