## Hipermutator strains analysis

We download a reference assembly for Escherichia coli (ASM584v2, 24/01/2025) and extract the genes involved in the mutator ability: 

```diff
+ # bash #

grep -iE -A 1 "mutS|mutL|mutH|dnaQ|uvrA|uvrB|uvrC|dam[^a-z]" mg1655_bakta/mg1655.ffn > mutator_genes.fasta

sed -i '/--/d' mutator_genes.fasta

# we use for the name of the contigs only the name of the gene

sed -i 's/^>\(.*\) \(.*\)$/>\2/' mutator_genes.fasta


```





