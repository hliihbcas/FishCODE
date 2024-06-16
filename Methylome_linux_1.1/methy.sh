#usage methy.sh <sample>  
SAMPLE=$2    # /home/yxjiang/mybiodb/WGBS/13_2.clean.R1.fastq.gz
BIOPROJECT=$1
GENOME=$3

METHRATIO='python ./methratio.py'
mkdir 01bsmap/${BIOPROJECT}
mkdir 02methyPos
#mkdir 03diff/${BIOPROJECT}
mkdir 02methyPos/01methratio/${BIOPROJECT}
#mkdir 02methyPos/03mCs/${BIOPROJECT}
wkdir=`pwd`
########## bsmap  methritio
N=`ls ${wkdir}/00_data/${BIOPROJECT}/${SAMPLE}*fastx |wc -l `
echo $N
echo "test"
if [[ $N -eq 2 ]];then
#echo "bsmap -a ./00_data/${BIOPROJECT}/${SAMPLE}_1_fastx -b ./00_data/${BIOPROJECT}/${SAMPLE}_2_fastx -d ${GENOME}  -o ./01bsmap/${BIOPROJECT}/${SAMPLE}.bam -v 0.1 -x 500 -m 50 -g 2 -R -u -n 0 -p 8"

#bsmap -a ./00_data/${BIOPROJECT}/${SAMPLE}_1_fastx -b ./00_data/${BIOPROJECT}/${SAMPLE}_2_fastx -d ${GENOME}  -o ./01bsmap/${BIOPROJECT}/${SAMPLE}.bam -v 0.1 -x 500 -m 50 -g 2 -R -u -n 0 -p 8
bsmap -a ./00_data/${BIOPROJECT}/${SAMPLE}_1_fastx -b ./00_data/${BIOPROJECT}/${SAMPLE}_2_fastx -d ${GENOME}  -o ./01bsmap/${BIOPROJECT}/${SAMPLE}.bam -R -p 8
elif [[ $N -eq 1 ]];then
echo "1"
#bsmap -a ./00_data/${BIOPROJECT}/${SAMPLE}_fastx -d ${GENOME}  -o ./01bsmap/${BIOPROJECT}/${SAMPLE}.bam -v 0.1 -x 500 -m 50 -g 2 -R -u -n 0 -p 8
bsmap -a ./00_data/${BIOPROJECT}/${SAMPLE}_fastx -d ${GENOME}  -o ./01bsmap/${BIOPROJECT}/${SAMPLE}.bam -R -p 8

fi
#wait
echo "bsmap is done"
samtools sort ./01bsmap/${BIOPROJECT}/${SAMPLE}.bam -o ./02methyPos/01methratio/${BIOPROJECT}/${SAMPLE}.sorted.bam
samtools index ./02methyPos/01methratio/${BIOPROJECT}/${SAMPLE}.sorted.bam
${METHRATIO} -d ${GENOME} -o ./02methyPos/01methratio/${BIOPROJECT}/${SAMPLE}.methratio   ./02methyPos/01methratio/${BIOPROJECT}/${SAMPLE}.sorted.bam
#methratio.py -d ${GENOME} -o ./02methyPos/01methratio/${BIOPROJECT}/${SAMPLE}.methratio ./02methyPos/01methratio/${BIOPROJECT}/${SAMPLE}.sorted.bam
sed -i '1s/#chr/chr/'   ./02methyPos/01methratio/${BIOPROJECT}/${SAMPLE}.methratio
#awk 'BEGIN{FS=OFS="\t";print "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT"}$9!="0"&&FNR!=1{if($3="+"){FLA="F"}else{FLA="R"};print $1"."$2,$1,$2,FLA,$9,$7/$9,1-$7/$9}' ./02methyPos/01methratio/${BIOPROJECT}/${SAMPLE}.methratio > ./02methyPos/01methratio/${BIOPROJECT}/${SAMPLE}.methratio.txt
awk 'BEGIN{FS=OFS="\t";print "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT"}$9!="0"&&FNR!=1{if($3="+"){FLA="F"}else{FLA="R"};print $1"."$2,$1,$2,FLA,$9,$7/$9*100,(1-$7/$9)*100}' ./02methyPos/01methratio/${BIOPROJECT}/${SAMPLE}.methratio > ./02methyPos/01methratio/${BIOPROJECT}/${SAMPLE}.methratio.txt
#Rscript methy_mCs.R ./02methyPos/01methratio/${BIOPROJECT}/${SAMPLE}.methratio ./02methyPos/02lambdaDNA/${BIOPROJECT}/${SAMPLE}.BSratio ./02methyPos/03mCs/${BIOPROJECT}/${SAMPLE}.pvalue

