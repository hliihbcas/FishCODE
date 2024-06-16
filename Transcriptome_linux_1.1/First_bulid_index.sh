# !/usr/bin/env bash
source ./pp.conf 

##GENOME_INDEX=$1 ##eg.cf, sm, gc, nt
#GENOME=${GENOME_Lepisosteus_oculatus}
#GFF=${GFF_Lepisosteus_oculatus}

GENOME=$1
GFF=$2

#cpu=30
cpu=$3

##decoy sequences
################# cf,sm,nt. [error] In FixFasta, two references with the same name but different sequences
######## check
##grep '^>' ${GFF%.*}_decoy/gentrome.fa |cut -d' ' -f1 |sort |uniq -c |awk '{if($1!=1)print}' |wc -l
##awk '{if($3!="mRNA" && $3!="transcript" && $3!="exon" && $3!="CDS" && $3!="gene" )print $3}' ${GFF} |sort |uniq -c
#####C_gene_segment
#####V_gene_segment
############ log end
awk 'BEGIN{FS="\t";OFS="\t"}{if($3=="C_gene_segment" || $3=="V_gene_segment" || $3=="D_gene_segment" || $3=="J_gene_segment")$3="mRNA";print $0}' ${GFF} > ${GFF%.*}_modify.gff
################ cf,sm,nt ,end
gffread ${GFF%.*}_modify.gff -g ${GENOME} -w ${GFF%.*}_modify_exons.fa 
./generateDecoyTranscriptome.sh -a ${GFF%.*}_modify.gff  -g ${GENOME} -t ${GFF%.*}_modify_exons.fa -o ${GFF%.*}_decoy -j $cpu 

##salmon:index
${salmon} index -p 20 -t ${GFF%.*}_decoy/gentrome.fa -i ${GFF%.*}_merged -k 23 --keepDuplicates -d ${GFF%.*}_decoy/decoys.txt 

##2.3.3 count by salmon(v1.0)
#gtffile="merged.annotated.gtf"
#gtffile="merged.annotated_split.gtf"

#cat sample.list|while read m
#do
#salmon quant -i ${GFF%.*}_merged -l A -1 filter_data/${m}_1.filter.fq.gz -2 filter_data/${m}_2.filter.fq.gz -o ${m}_quant -g $gtffile -p 20 --mimicBT2 --useEMoe
#done
