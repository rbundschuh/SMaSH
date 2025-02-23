# SMaSH
A tool for sample swap identification in high throughput sequencing studies

### Introduction ###

Sample swaps are a real concern in high throughput sequencing studies.  SMaSH helps detect such sample swaps by integrating the information from over 6000 carefully selected single nucleotide polymorphism (SNP) sites from across the human genome to identify which samples in a group of sequencing data sets are derived from the same human individual.  Importantly, SMaSH is able to verify sample identity between different data types, such as RNA-Seq, exome, and MethylCap-Seq data.

### Usage ###

SMaSH is written in Python and is run by simply typing:

    SMaSH.py sample1.bam sample2.bam [sample3.bam [etc...]]

This requires the file snps_hg19.vcf to be stored in the directory of the SMaSH.py.  The bam files must be indexed and must result from alignment to human genome hg19 (for GRCh38 see below). Sam files will be very slow if they are not very small.

The output will be written to the file pval\_out.txt in the form of a matrix of p-values for the null hypothesis that the samples are from different individuals for all pairwise comparisons of the samples in the directory.  Thus, low p-values indicate high probability that the corresponding samples are derived from the same individuals; p-values close to one indicate that the corresponding samples are derived from different individuals.  P-values for samples to be from the same individual are written to the file issameindividual\_pvalue\_out.txt. The latter correspond to the p-values described in the paper presenting SMaSH.  While these are just one minus the p-values reported in pval\_out.txt, reporting them separately allows distinguishing between very small p-values, where reporting 1-p results in rounding to one.

    SMaSH.py ALL

Looks for *.bam and *.cram in the current directory.  

Various command line options are explained by running:

    SMaSH.py --help

by itself and allow choosing specific sam, bam, or cram files, providing different SNP index files, e.g. `snps_GRCh38.vcf` for data aligned to GRCh38, and/or choosing a different output file.

```
usage: SMaSH.py [-h] [-bam BAM_COMMA_LIST] [-i INFILE] [-o OUTNAME]
                [-chr_index CHR_INDEX] [-pos_index POS_INDEX]
                [-ref_index REF_INDEX] [-alt_index ALT_INDEX] [-regenerate]
                [-output_dir OUTPUT_DIR] [-include_rgid]
                [bam [bam ...]]

positional arguments:
  bam                   BAM/SAM/CRAM files to check. Note BAMs must end in .bam and
                        be indexed

optional arguments:
  -h, --help            show this help message and exit
  -bam BAM_COMMA_LIST, --bam BAM_COMMA_LIST
                        [deprecated] input BAM/SAM/CRAM files to be tested (comma
                        separated list) or 'ALL' to use all BAMs and CRAMs in
                        current dir.
  -i INFILE, --positions INFILE
                        input locations file
  -o OUTNAME, --output_file OUTNAME
                        output file name [pval_out.txt]
  -chr_index CHR_INDEX, --chr_index CHR_INDEX
                        index of chromosome column in locations file [0]
  -pos_index POS_INDEX, --pos_index POS_INDEX
                        index of position column in locations file [1]
  -ref_index REF_INDEX, --ref_index REF_INDEX
                        reference allele index
  -alt_index ALT_INDEX, --alt_index ALT_INDEX
                        alternate allele index
  -regenerate, --regenerate
                        regenerate SNP read counts and all calculations
  -output_dir OUTPUT_DIR, --output_dir OUTPUT_DIR
                        The directory to save output files. [default: ./]
  -include_rgid, --include_rgid
                        include BAM's Read Group ID value in output and only
                        print the BAM's basename, not full path

```

### Dependencies ###


SMaSH requires the following python libraries:
    pysam
    scipy
    numpy

### Citation ###

If you use SMaSH in your project, please cite

Maximillian Westphal, David Frankhouser, Carmine Sonzone, Peter G. Shields, Pearlly Yan, and Ralf Bundschuh. SMaSH: Sample matching using SNPs in humans. _BMC Genomics_ **20**(Suppl 12): 1001 (2019). https://doi.org/10.1186/s12864-019-6332-7.

