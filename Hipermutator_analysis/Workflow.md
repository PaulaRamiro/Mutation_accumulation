## Hipermutator strains analysis

We download a reference assembly for Escherichia coli (ASM584v2, 24/01/2025) and extract the genes involved in the mutator ability: 

```diff
+ # bash #

grep -iE -A 1 "mutS|mutL|mutH|mutT|mutY|mutM|dnaQ|uvrA|uvrB|uvrC" mg1655_bakta/mg1655.ffn > mutator_genes.fasta

sed -i '/--/d' mutator_genes.fasta

# we use for the name of the contigs only the name of the gene

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

```


