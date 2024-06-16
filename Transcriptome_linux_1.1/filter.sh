#!/usr/bin/env bash
source ./pp.conf
##PATH_prefetch=/home/yduan/soft/bio/web/sratoolkit.2.9.2-centos_linux64/bin
##PATH_Aspera=/home/wtzhang/.aspera/connect/bin
### input file, a list file of SRA accessions.

PROJECT=$1 ##eg.cf_rnaseq, gc_bsseq
SRR=$2
wkdir=`pwd`

mkdir ${wkdir}/00_data
mkdir ${wkdir}/00_data/${PROJECT}

### Filter data
N=`ls ${wkdir}/00_data/${PROJECT}/${SRR}*fastq.gz |wc -l `
echo $N
if [[ $N -eq 2 ]];then
  IlluQC.pl -pe ${wkdir}/00_data/${PROJECT}/${SRR}_1.fastq.gz ${wkdir}/00_data/${PROJECT}/${SRR}_2.fastq.gz 2 A -z g
elif [[ $N -eq 1 ]];then
  IlluQC.pl -se ${wkdir}/00_data/${PROJECT}/${SRR}*.gz  2 A -z g
elif [[ $N -eq 3 ]];then
  IlluQC.pl -pe ${wkdir}/00_data/${PROJECT}/${SRR}_1.fastq.gz ${wkdir}/00_data/${PROJECT}/${SRR}_2.fastq.gz 2 A -z g
fi 
## Remove raw fastq data
FILTERED_SIZE=`du ${wkdir}/00_data/${PROJECT}/IlluQC_Filtered_files/${SRR}*_filtered.gz | awk '{print $1}' |head -1`
if [[ ${FILTERED_SIZE} -ge 3072 ]] ;then
  rm ${wkdir}/00_data/${PROJECT}/${SRR}*fastq.gz 
fi 
