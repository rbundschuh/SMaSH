# SMaSH
A tool for sample swap identification in high throughput sequencing studies

### Introduction ###

Sample swaps are a real concern in high throughput sequencing studies.  SMaSH helps detect such sample swaps by integrating the information from over 6000 carefully selected single nucleotide polymorphism (SNP) sites from across the human genome to identify which samples in a group of sequencing data sets are derived from the same human individual.  Importantly, SMaSH is able to verify sample identity between different data types, such as RNA-Seq, exome, and MethylCap-Seq data.

### Usage ###

SMaSH is written in python 2.7 and is run by simply typing

    SMaSH.py -bam ALL

This requires the file snps_hg19.vcf to be stored in the current directory and will automatically compare all bam files in the current directory.  The bam files must be indexed and must result from alignment to human genome hg19.  The output will be written to the file pval\_out.txt in the form of a matrix of p-values for all pairwise comparisons of the samples in the directory.  Low p-values indicate high probability that the corresponding samples are derived from the same individuals; p-values close to one indicate that the corresponding samples are derived from different individuals.

Various command line options are explained by running

    SMaSH.py

by itself and allow choosing specific bam files, providing different SNP index files, such as, e.g., the file snps_GRCh38.vcf for data aligned to GRCh38, and/or choosing a different output file.
