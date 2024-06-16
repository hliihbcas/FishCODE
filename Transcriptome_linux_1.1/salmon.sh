# !/usr/bin/env bash
source ./pp.conf

PROJECT=$1 ##eg.cf_rnaseq, gc_bsseq
SRR=$2
wkdir=`pwd`

#GENOME_INDEX=$3 ##eg.GENOME_cf, GENOME_sm, GENOME_gc, GENOME_nt
#GFF_INDEX=$4 ##eg.GFF_cf,...
GENOME=$3
GFF=$4
echo "i am here"
echo ${GENOME}
echo ${GFF}
mkdir ${wkdir}/01_result
mkdir ${wkdir}/01_result/${PROJECT}
threads=20

##decoy sequences
################# cf,sm,nt. [error] In FixFasta, two references with the same name but different sequences
######## check
##grep '^>' ${GFF%.*}_decoy/gentrome.fa |cut -d' ' -f1 |sort |uniq -c |awk '{if($1!=1)print}' |wc -l
##awk '{if($3!="mRNA" && $3!="transcript" && $3!="exon" && $3!="CDS" && $3!="gene" )print $3}' ${GFF} |sort |uniq -c
#####C_gene_segment
#####V_gene_segment
############ log end
##awk '{if($3=="C_gene_segment" || $3=="V_gene_segment")$3="mRNA";print $0}' ${GFF} > ${GFF%.*}_modify.gff 
################ cf,sm,nt ,end

### 
#gffread ${GFF%.*}_modify.gff -g ${GENOME} -w ${GFF%.*}_modify_exons.fa 
#generateDecoyTranscriptome.sh -a ${GFF%.*}_modify.gff  -g ${GENOME} -t ${GFF%.*}_modify_exons.fa -o ${GFF%.*}_decoy -j ${threads} 

##salmon:index
#${salmon} index -p ${threads} -t ${GFF%.*}_decoy/gentrome.fa -i ${GFF%.*}_merged -k 23 --keepDuplicates -d ${GFF%.*}_decoy/decoys.txt 

##2.3.3 count by salmon(v1.0)
#cat sample.list|while read m
#do
N=`ls ${wkdir}/00_data/${PROJECT}/IlluQC_Filtered_files/${SRR}*.gz |wc -l `
if [[ $N -ge 2 ]];then 
  ${salmon} quant --gcBias -i ${GFF%.*}_merged -l A -1 ${wkdir}/00_data/${PROJECT}/IlluQC_Filtered_files/${SRR}_1.fastq.gz_filtered.gz -2 ${wkdir}/00_data/${PROJECT}/IlluQC_Filtered_files/${SRR}_2.fastq.gz_filtered.gz  -o ${wkdir}/01_result/${PROJECT}/${SRR}_quant  -p ${threads} -g ${GFF%.*}_modify.gff --mimicBT2 --useEM  ##-g ${GFF%.*}_decoy/gentrome.fa
elif [[ $N -eq 1 ]];then
  ${salmon} quant --gcBias -i ${GFF%.*}_merged -l A -r ${wkdir}/00_data/${PROJECT}/IlluQC_Filtered_files/${SRR}.fastq.gz_filtered.gz -o ${wkdir}/01_result/${PROJECT}/${SRR}_quant  -p ${threads} -g ${GFF%.*}_modify.gff --mimicBT2 --useEM ##
fi 
mv ${wkdir}/00_data/${PROJECT}/IlluQC_Filtered_files/${SRR}* /home/hli/mybiodb/GBTP/00_data/collect_IlluQC/${PROJECT}
echo ${SRR}" is done"
#done
