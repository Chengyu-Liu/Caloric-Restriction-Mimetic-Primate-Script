#!/bin/bash
#This is a pipeline for ssRNA-seq analysis
#Current version 1: J.Q.Z, 2021-10-27 modified by JXY 2022-05-05

array=( "$@" )

#Usage
if [ ! -n "$1" ]
then
  echo "********************************************************************************************"
  echo "*                     PipeRiboseq: pipeline for RNA-seq analysis.                         *"
  echo "*                            Version 2, 2021-10-27, Q.Z                                    *"
  echo "* Usage: `basename $0`                                                                     *"
  echo "*        Required (paired end data):                                                       *"
  echo "*                  -F [fastq files1]                                                       *"
  echo "*                  -f [fastq files2]                                                       *"
  echo "*                  -n [The prefix samplename]                                              *"
  echo "*                  -o [The outputpath]                                                     *"
  echo "*        Optional: -p [Number of CPUs, default=24]                                         *"
  echo "*                  -g [The file of Gene GTF files]                                         *"
  echo "*                  -G [The file of REs GTF files]                                          *"
  echo "*                  -i [The file of hisat2 index files]                                     *"
  echo "*                  -I [The file of STAR index files]                                       *"
  echo "* Inputs: The raw fastq files                                                              *"
  echo "* Run: Default to run trim_galore,fastqc,hisat2,sam2bam,featureCounts & deeptools          *"
  echo "*      Figures will be generated in /plots folder, and bigWig files in /tracks folder      *"
  echo "* Outputs: All output files will be generated in the same folder as the pipeline submitted *"
  echo "* This pipeline requires 'conda activate daily'                                            *"
  echo "********************************************************************************************"
  exit 1
fi

#Get parameters
fq1="unassigned"
fq2="unassigned"
sample="unassigned"
out="unassigned"
CPU=24
gtf="unassigned"
REgtf="unassigned"
index="unassigned"
STARindex="unassigned"

for arg in "$@"
do
 if [[ $arg == "-F" ]]
  then
    fq1=${array[$counter+1]}
    echo '   raw fastq1: '$fq1
 elif [[ $arg == "-f" ]]
  then
    fq2=${array[$counter+1]}
    echo '   raw fastq2: '$fq2
 elif [[ $arg == "-n" ]]
  then
    sample=${array[$counter+1]}
    echo '   Sample name: '$sample
 elif [[ $arg == "-o" ]]
  then
    out=${array[$counter+1]}
    echo '   The output of all files: '$out
 elif [[ $arg == "-p" ]]
  then
    CPU=${array[$counter+1]}
    echo '   CPU: '$CPU
 elif [[ $arg == "-g" ]]
  then
    gtf=${array[$counter+1]}
    echo '   GTF: '$gtf
 elif [[ $arg == "-G" ]]
  then
    REgtf=${array[$counter+1]}
    echo '   REGTF: '$REgtf
 elif [[ $arg == "-i" ]]
  then
    index=${array[$counter+1]}
    echo '   Index: '$index
 elif [[ $arg == "-I" ]]
  then
    STARindex=${array[$counter+1]}
    echo '   STARIndex: '$STARindex
 fi
  let counter=$counter+1
done

echo "*****************************************************************************"
echo "1. Generate the all output files"

cd $out
mkdir -p 02_trim 03_mapping 04_count
trim=$out/02_trim
map=$out/03_mapping
count=$out/04_count

for path in `ls $out | grep 0`
do
mkdir -p $path/logs
done

echo "*****************************************************************************"
echo "2. Run trim_galore"

trim_log=$trim/logs/${sample}.log
trim_out=$trim/$sample
mkdir -p $trim_out
trim_galore --fastqc --path_to_cutadapt /share/home/jiangxiaoyu/anaconda3/bin/cutadapt --stringency 3 --paired --output_dir $trim_out $fq1 $fq2 2>$trim_log

echo "*****************************************************************************"
echo "3. Mapping using STAR"

star_log=$map/logs/${sample}.star.log
map_out=$map/$sample
mkdir -p $map_out
# input ###
clean1=$trim_out/*1.fq.gz
clean2=$trim_out/*2.fq.gz
STARout=$map_out/STAR
mkdir -p $STARout

~/anaconda3/bin/STAR --runMode alignReads --readFilesCommand zcat --genomeDir $STARindex --runThreadN 24 --readFilesIn $clean1 $clean2 --outFileNamePrefix $STARout/${sample}.Gene. --alignSJDBoverhangMin 1 --outReadsUnmapped Fastx --outFilterMismatchNoverLmax 0.05 --outFilterScoreMinOverLread 0.90 --outFilterMatchNminOverLread 0.90 --alignIntronMax 1000000 --outSAMtype BAM SortedByCoordinate 2>$star_log
samtools index $STARout/${sample}.Gene.Aligned.sortedByCoord.out.bam 2>>$star_log

echo "*****************************************************************************"
echo "4. assigning sequence reads to genomic features"

count_log=$count/logs/${sample}.log
count_out=$count/$sample
mkdir -p $count_out

featureCounts -o $count_out/${sample}.rna.STAR.txt -Q 10 -p -T $CPU -F GTF -a $gtf -t exon -g gene_id -s 2 -B -C --minOverlap 60 $STARout/${sample}.Gene.Aligned.sortedByCoord.out.bam 2>>$count_log
awk '{print $1,$7}' $count_out/${sample}.rna.STAR.txt | grep -v "#" > $count_out/${sample}.STAR_final.txt
