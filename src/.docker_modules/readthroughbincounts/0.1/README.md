# ReadthroughBinCounts

Comptage des readthrough dans des fichiers BAM, à partir de coordonnées BED

## Step 1 : GTF select 'gene' lines 

```
grep -P "\tgene\t" Homo_sapiens.GRCh37.75.gtf > Homo_sapiens.GRCh37.75.genes.gtf
```

## Step 2 : bed pour avoir le gene puis le gene + 10kb 

Packages : 
- tidyverse 
- progress 
- GenomicRanges 
- valr 
- rtracklayer 
- optparse 
- foFuture 
- DESeq2 
- viridis 
- gplots 
- openxlsx


```
Rscript --vanilla ./src/refToKallisto.R --help 
Usage: ./src/refToKallisto.R [options]


Options:
        -f CHARACTER, --file=CHARACTER
                GTF File

        -o CHARACTER, --out=CHARACTER
                output repertory [default= ./]

        -s CHARACTER, --suffix=CHARACTER
                A suffix to specify your results [default= readtrhough]

        -r NUMBER, --RT=NUMBER
                Size of Readthrough search zone after last gene exon (#bases) [default= 10000]

        -c CHR, --chr=CHR
                List of chromosomes you want to keep (default: human chr) 
                [default= 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,X,Y,MT]

        -h, --help
                Show this help message and exit


Rscript --vanilla ./src/refToKallisto.R -f "./data/Homo_sapiens.GRCh37.75.genes.gtf" -o "./results/" 

Loading GTF file ( 1 / 6 )
Export classical Bed file ( 2 / 6 )
Create 10000kb extension ( 3 / 6 )
Overlap on genes : coordinate correction ( 4 / 6 )
 [==================================] 100% [Elapsed time: 00:00:18 || Estimated time remaining:  0s]
Overlap on genes : Stranded correction ( 5 / 6 )
 [==================================] 100% [Elapsed time: 00:01:15 || Estimated time remaining:  0s]
Export Extension Bed file ( 6 / 6 )

```
/!\ faire sur le dernier exon et non tout le gène, reprendre data Xavier 
Version - last exon ( All bases ) 


## Step 3 : création d'un fasta file pour kallisto 

```
bedtools getfasta -name -fi ~/[PATH]/genomes/Homo_sapiens.GRCh37.dna.primary_assembly.fa -bed GenesAndRT_10000bases_tete.bed -fo Homo_sapiens.GRCh37.GenesAndRT.fasta
```

## Step 4 : index & quant kallisto (PSMN) v46.2 

```
kallisto index -i GRCh37_GenesAndRT  ../genome/Homo_sapiens.GRCh37.GenesAndRT.fasta

kallisto quant -i kallisto_idx/GRCh37_GenesAndRT -t 32 -o kallisto/5Y_siDDX5-17_B1 fastq/5Y_siDDX5-17_B1_R1_cutadapt_match.fastq.gz fastq/5Y_siDDX5-17_B1_R2_cutadapt_match.fastq.gz

kallisto quant -i kallisto_idx/GRCh37_GenesAndRT -t 32 -o kallisto/5Y_siDDX5-17_B2 fastq/5Y_siDDX5-17_B2_R1_cutadapt_match.fastq.gz fastq/5Y_siDDX5-17_B2_R2_cutadapt_match.fastq.gz

kallisto quant -i kallisto_idx/GRCh37_GenesAndRT -t 32 -o kallisto/5Y_siDDX5-17_B3 fastq/5Y_siDDX5-17_B3_R1_cutadapt_match.fastq.gz fastq/5Y_siDDX5-17_B3_R2_cutadapt_match.fastq.gz


kallisto quant -i kallisto_idx/GRCh37_GenesAndRT -t 32 -o kallisto/5Y_siGL2_B1 fastq/5Y_siGL2_B1_R1_cutadapt_match.fastq.gz fastq/5Y_siGL2_B1_R2_cutadapt_match.fastq.gz

kallisto quant -i kallisto_idx/GRCh37_GenesAndRT -t 32 -o kallisto/5Y_siGL2_B2 fastq/5Y_siGL2_B2_R1_cutadapt_match.fastq.gz fastq/5Y_siGL2_B2_R2_cutadapt_match.fastq.gz

kallisto quant -i kallisto_idx/GRCh37_GenesAndRT -t 32 -o kallisto/5Y_siGL2_B3 fastq/5Y_siGL2_B3_R1_cutadapt_match.fastq.gz fastq/5Y_siGL2_B3_R2_cutadapt_match.fastq.gz

```

## Step 5 : Formattage table de comptage et Detection des Readthrough par DESeq2

```
Rscript --vanilla ./src/DEG_RT.R --help 
Usage: ./src/DEG_RT.R [options]

Options:
        -d CHARACTER, --dir=CHARACTER
                Directory with kallisto folders

        -o CHARACTER, --out=CHARACTER
                output repertory [default= ./]

        -f NUMBER, --lfch=NUMBER
                log2FoldChange Threshold [default= 4]

        -t THREADS, --threads=THREADS
                Number of process parallelization [default= 1]

        -h, --help
                Show this help message and exit


Rscript --vanilla ./src/DEG_RT.R -d "./data/kallisto/siDDX/" -o "./results/siDDX/" -t 4 -f 4
 - Formatting input files
 - DEG analysis
converting counts to integer mode
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
   

|Log2FoldChange| Threshold > 4
Number of genes with Readthrough : 3396
 - Plots
 
 Rscript --vanilla ./src/DEG_RT.R -d "./data/kallisto/siGL2/" -o "./results/siGL2/" -t 4 -f 4
 - Formatting input files
 - DEG analysis
converting counts to integer mode
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
   

|Log2FoldChange| Threshold > 4
Number of genes with Readthrough : 3353
 - Plots
 
 
```

## Step 5b : Differentiel entre les RT d'une condition (eg siDDX) vs Ctrl (eg siGL2) 


```
Rscript --vanilla ./src/DEG_2cond.R --help 
Usage: ./src/DEG_2cond.R [options]


Options:
        -a CHARACTER, --cond1=CHARACTER
                csv file with raw counts for the condition 1

        -b CHARACTER, --cond2=CHARACTER
                csv file with raw counts for the condition 2 (control condition)

        -o CHARACTER, --out=CHARACTER
                output repertory [default= ./]

        -f NUMBER, --lfch=NUMBER
                log2FoldChange Threshold [default= 0.4]

        -h, --help
                Show this help message and exit
                

Rscript --vanilla ./src/DEG_2cond.R -a "./results/siDDX/rawcounts.csv" -b "./results/siGL2/rawcounts.csv" -o "./results/siDDX-vs-siGL2/" -f 0.4
 - Formatting input files
 - DEG analysis
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
   

|Log2FoldChange| Threshold > 0.4
Number of genes with Readthrough : 1345
 - Plots
 - Done

 
 
``` 

N.B : On retrouve significatifs les marqueurs vaidés dans [Terrone et al, 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9458439/), Figure 1B : 
- NCS1 
- PPARD 
- RAB36 
- NCKAP5L 

Le gène SH3TC1 n'a pas de RT significativement variable entre siDDX vs siGL2

/!\ créer un argument pour ajouter un coldata (meme si valeur pa defaut)
/!\ sortie finale avec colonnes moyenne d'expression du last exon sans RT ( basemean pour chaque condition, dc 2 colonnes )
/!\ sortie finale avec colonnes moyenne d'expression ( basemean pour chaque condition, dc 2 colonnes ) du gene total sans RT (mais pas dans l'analyse stat) 



## Step 6 : Selection des genes diff sig & comptage des reads / bases 

1 - selection des genes DE sig, 
2 - BamCoverage sur les RT pour un vecteur de couverture sur les RT
3 - Add Xavier Grand as dev




## Step 7 : Bi-segmentation (fin du RT) 


## Step 8 : Ratio & formattage final 



