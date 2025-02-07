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

And we trim the reads using Trim-galore :






Then, we download the BioSample accession from the BioProject (PRJNA294072) metadata downloaded in SRA Run Selector (https://www.ncbi.nlm.nih.gov/Traces/study/) and download the assemblies by using : 

```diff


```
