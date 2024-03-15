# Metagenomic analysis of coral associated microbiome using DADA2

This repository contain R scripts used to reproduce the results of this paper: 
Nhung, D. T., & Van Ngoc, B. (2020). Bioinformatic approaches for analysis of coral-associated bacteria using R programming language. Vietnam Journal of Biotechnology, 18(4), 733-743.

### Pipeline

The standard **DADA2** and **phyloseq** libraries and their associated pipelines were applied to analyse the composition and diversity of microbiome associated with coral reefs in certain ecosystems.

In-detail intruction for the DADA2 package can be found at: https://benjjneb.github.io/dada2/tutorial_1_8.html

16S rRNA sequencing files (originally in .SRA format) were downloaded from NCBI according to the instructions: [](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/)
After downloading, all files were converted to FASTQ format as input for the DADA2 pipeline using the following command:

```
dir > list.txt # append all file names into a text file
for i in $(cat list.txt); do echo $i; date; fasterq-dump -S $i; done
```

Details of R and the used 16S rRNA sequences as well as the locations and samples will be updated later
