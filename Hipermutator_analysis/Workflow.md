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
