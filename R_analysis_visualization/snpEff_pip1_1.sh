SPECIES=$1 ##eg. GrassCarp, commoncarp
POSITIVE=$2
PARA=$3
PARA2=`echo ${PARA}|sed 's/,/ \\-/g'|sed 's/\[//g'|sed 's/\]//g'`
echo ${PARA2}

#OUT_=`pwd`/snpEffout
#OUT=${OUT_}/$3
OUT=$4 #可以直接引用自动建立的目录，临时存储结果文件，
cd ${OUT}
# mkdir ${OUT}
###1.在~/snpEff/目录中，创建一个文件夹：data
###2.在~/snpEFF/data目录下，创建两个文件夹
###    AT_10/   genomes/
###    这两个文件夹中，分别放置了GTF文件和基因组文件
###    genes.gtf sequences.fa
###3.编辑~/snpEff/snpEff.config文件
###    在文件的最后一行添加信息：
###    AT_10.genome: AT

#java -jar ~/snpEff/snpEff.jar build -c ~/snpEff/snpEff.config -gtf22 -v AT_10

###参数说明
###java -jar: Java环境下运行程序
###-c snpEff.config配置文件路径
###-gtf22 设置输入的基因组注释信息是gtf2.2格式
###-gff3 设置输入基因组注释信息是gff3格式
###-v 设置在程序运行过程中输出的日志信息
###最后的AT_10参数 设置输入的基因组版本信息，和~/snpEff/snpEff.config配置文件中添加的信息一致
# snpEff eff -c /home/lzhang/biosoft/snpEff/snpEff.config ${SPECIES} ${POSITIVE} > Effresult.snp.eff.vcf -csvStats Effresult.csv -stats Effresult.html


java -Xmx8g -jar /home/lzhang/biosoft/snpEff/snpEff.jar  -c /home/lzhang/biosoft/snpEff/snpEff.config  ${SPECIES} ${POSITIVE} \-${PARA2} > Effresult.snp.eff.vcf -csvStats Effresult.csv -stats Effresult.html
mkdir Result
mv Effresult* Result
zip -r Effresult.zip Result #会保留原始的文件

# rm -rf Result Effresult.zip





#java -Xmx8g -jar snpEff.jar GRCh37.75 examples/test.chr22.vcf > test.chr22.ann.vcf
#java -Xmx8g -jar ~/soft/snpEff/snpEff.jar eff -c /home/Xiaxq/soft/snpEff/snpEff.config ${SPECIES} ${POSITIVE} > positive.snp.eff.vcf -csvStats positive.csv -stats positive.html
###最终会产生四个文件 positive.snp.eff.vcf positive.html positive.csv positive.genes.txt
###可以在positive.snp.eff.vcf文件中，分析自己的后续基因位点了。
#java -jar /home/wtzhang/bin/snpEff.jar  -v commoncarp 100children2parents_snp.pass.vcf
