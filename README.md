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

## Transcriptome upstream script
### Linux environment  
1.Before the process begins, you need to make sure that the following dependent software in your Linux environment is installed successfully and can be used properly.    
gffread, salmon 1.4.0, bedtools v2.26.0, mashmap, GSQCToolkit_v2.3.3 (IlluQC.pl), salmon 1.4.0  
2. The pp.conf file records the absolute paths of the dependent software mentioned below. Before you start you must configure this file and replace the execution path of the software with the absolute path to the software on your server!  

### Create database index for analysis
$ bash First_bulid_index.sh ./build_faidx_example/GCF_901000725.2_fTakRub1.2_genomic.fna ./build_faidx_example/GCF_901000725.2_fTakRub1.2_genomic.gff 30  
Note: 1. Your genome and corresponding annotation files must be in the same folder.   
2. All the following working directories default to the decompressed directory.  
3. The first parameter is the genome location, the second parameter is the gff file location, and the third parameter is the number of threads.  

### Quality control
$ bash filter.sh 792045 SRR17641338  
Note: The first parameter is the name of the upper folder where fastq.gz is located, or the project name, and the second parameter is the prefix of fastq.gz. For the specific file structure, please refer to the compressed package. The user must create a project under the parent directory 00_data, and store the sequencing files to be processed under the project name. "./00_data/projectname/*fastq.gz".  

### Count
$ bash salmon.sh 792045 SRR17641338 ./build_faidx_example/GCF_901000725.2_fTakRub1.2_genomic.fna ./build_faidx_example/GCF_901000725.2_fTakRub1.2_genomic.gff  
Note: The first parameter is the project name of the second step, the second parameter is the prefix of fastq.gz, the third parameter is the genome location, and the fourth parameter is the location of the gff file (note that the location of this genome and annotation file must Strictly the input genome and gff file path for the first step)  

**Before you start the process, here are a few things you need to know**  
1. The output result of the third step will be in 01_result/792045/SRR17641338_quant/ under the same directory.  
2. All directory structures and demonstration data examples are included in the compressed package, allowing users to try it easily.  
3. In order to facilitate users to generate transcriptome data corresponding to the species included in FishCODE, we provide a download package of the genome and annotation files of FishCODE's 35 fish species in the "Help->Usage->gene information" section (http://bioinfo.ihb.ac.cn/fishcode/help).  

## Genome script  
please refer to https://github.com/hliihbcas/FishSNP

## Analysis and visualization code
### Dependent R packages
1.Before the process begins, you need to make sure that the following dependent software in your Linux environment is installed successfully and can be used properly.    
please refer FishCODE.yaml.  
2.If conda works under your Linux environment, you canimport the corresponding environment.  
$ create -f FishCODE.yaml via conda env create -f  
3. Enter the corresponding directory and grant execute permission  
$ cd R_analysis_visualization  
$ chmod 755 ./*R  

### Transcriptome analysis
$ ./Rversion1.R de_analyse_batch 11 ./data/example_trans/*.data 0 ./data/example_trans/feature > ./data/example_trans/err.txt 2>&1  

### Time-series transcriptome analysis
$ ./Rversion1.R de_analyse_batch_time 11 ./data/example_time/*.data 0 ./data/example_time/feature > ./data/example_time/err.txt 2>&1  

### Evolutionary tree construction  
$ Rscript ./evolution_tree.R ./data/example_evo/common_aln.phy

### SNP annotation  
$ ./run_snp2.py -sh ./snpEff_pip1_1.sh -f ./data/example_snp/common.data -i new_grasscarp -c --address xxxxxxxxx@xxx.com -w ./data/example_snp -o ./data/example_snp/out  
Note: Considering that SNP annotations need to be built in advance, and that SNP annotations mainly rely on SNPEff software (https://pcingola.github.io/SnpEff), the usage is also relatively simple. This script is mainly used for review and reference.

### Methylation analysis  
$ Rscript ./methylkit.R SRR17851372,SRR17851373 SRR17851372,SRR17851373 3 1000 1000 0.01 25 all ./data/example_meth/result_1.txt 1,0 10 ./data/example_meth no 0.2 0.01 PRJNA802599,PRJNA802599  
Note: Differential methylation analysis.  

$ ./find_genes_by_loc.py --gene-file=./Cyprinus_carpio.txt --gene-name=Gene --gene-chr=Chr --gene-start=Start --gene-end=End --loc-chr=chr --loc-start=start --loc-end=stop -i ./data/example_meth/result_1.txt -o ./data/example_meth/SRR17851372,SRR17851373.meth.gene.CG.txt  
Note: Annotation of genes directly associated with differentially methylated regions.  

$ ./find_genes_nearby_loc.py --gene-file=./Cyprinus_carpio.txt --gene-name=Gene --gene-chr=Chr --gene-start=Start --gene-end=End --loc-chr=chr --loc-start=start --loc-end=stop  --adjacent=2000 -i ./data/example_meth/SRR17851372,SRR17851373.meth.gene.CG.txt -o ./data/example_meth/SRR17851372,SRR17851373.meth.2000.CG.txt  
Note：Annotation of (~2000bp) genes indirectly associated with differentially methylated regions.  
