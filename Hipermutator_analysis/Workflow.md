## Hipermutator strains analysis

We download a reference assembly for Escherichia coli (ASM584v2, 24/01/2025) and extract the genes involved in the mutator ability: 

```diff
+ # bash #

grep -iE -A 1 "mutS|mutL|mutH|dnaQ|uvrA|uvrB|uvrC|dam[^a-z]" mg1655_bakta/mg1655.ffn > mutator_genes.fasta

sed -i '/--/d' mutator_genes.fasta

# we use for the name of the contigs only the name of the gene

sed -i 's/^>\(.*\) \(.*\)$/>\2/' mutator_genes.fasta


```

We search for mutations in each gene and each of our XXX samples with blast: 


```diff
+ # bash #

# Database with E.coli genomes
makeblastdb -in tus_secuencias.fasta -dbtype nucl -out mi_base_datos
blastn -query genes.fasta -db mi_base_datos -out resultados_blast.txt -outfmt 6


# We run blastn for each gene
mkdir resultados_blast

for gene in $(grep ">" mutator_genes.fasta | sed 's/>//'); do
    grep -A 1 $gene mutator_genes.fasta > ${gene}.fasta
    blastn -query ${gene}.fasta -db mi_base_datos -out resultados_blast/${gene}_blast.txt -outfmt 6
done

```
