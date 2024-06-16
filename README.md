# FishCODE
FishCODE: a web-based information platform for comprehensive omics data exploration in fish research.  
The script should be run in the Linux environment.  

## Epigenomics upstream script  
### Linux environment  
1.Before the process begins, you need to make sure that the following dependent software in your Linux environment is installed successfully and can be used properly.    
FASTX Toolkit 0.0.13, BSMAP 2.9, samtools 1.3.1, Python 2.7.15 [sys, time, os, array, optparse]  

### Quality control
$ bash filter_fastx.sh 555065 SRR9697472  
Note: The first parameter of the quality control script is the project name, and the second parameter is the prefix name of fastq.gz. But before you run the script, you need to replace "@SRR" on line 20 of "filter_fastx.sh" with the special flag string of your sequencing data. What are special identifiers? Take NCBI's original sequencing file as an example. The special string of SRRXXXX.fastq.gz is "@SRR". You can zless open the fastq.gz file to view the first line and replace the "NCBI's" in line 20 of filter_fastx.sh. @SRR".  

### Count
$ bash methy.sh 555065 SRR9697472 ./zebrafish/GCF_000002035.6_GRCz11_genomic.fna  
Note: The first parameter of the quality control script is the project name, the second is the prefix name of fastq.gz, and the third parameter is the absolute path of the genome of the corresponding species.  

**Before you start the process, here are a few things you need to know**  
1. The output result of the third step will be in 02methyPos/01methratio/555065 under the same directory.  
2. All directory structures and demonstration data examples are included in the compressed package, allowing users to try it easily.  
3. *methratio.txt is the final output result, where *methratio is the methylation status of all sites, and the final result is the filtered site file with coverage.  
4. Considering that methylation computing resources consume a lot, we recommend that your running server has at least 100G of RAM, 200G of storage space, and 10 CPU cores. (actually depends on your genome size).  

**Getting started (Linux)**  
$ cd ./SNP_calling_pipeline  
$ snakemake -s SNP_calling_pipeline_snakemake.py --cores 2

## Mendelian test
The process includes the Mendelian test.  

**Before you start the snakemake process, here are a few things you need to do**  
1.gunzip ./Mendelian_test/example_data/multiple_sort_vcf.vcf.gz

**Getting started (Linux)**   
$ cd Mendelian_test  
$ snakemake -s mendelian_test_snakemake.py --cores 1

## Hardy Weinberg test
The process includes the Hardy Weinberg test.  
**Getting started (Linux)**  
$ cd Hardy_Weinberg_test  
$ snakemake -s Hardy_Weinberg_test/hardy_pipeline_snakemake.py --cores 1
