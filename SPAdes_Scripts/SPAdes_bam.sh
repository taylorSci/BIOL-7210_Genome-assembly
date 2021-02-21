#!/bin/bash


mkdir $outputDir/SPAdes/bam

files=($(ls $outputDir/fastp/fastq_output/))

cd $outputDir/SPAdes/extra

for j in "${files[@]}";
do
  bwa index $outputDir/SPAdes/contigs/${j}_SPAdes.fasta
  bwa mem $outputDir/SPAdes/contigs/${j}_SPAdes.fasta $files/$j/${j}_1_fp.fq.gz $files/$j/${j}_2_fp.fq.gz > ${j}.sam
  samtools fixmate -O bam ${j}.sam ${j}.bam
  samtools sort -O bam -o ${j}_sorted.bam -T temp ${j}.bam
  cp ${j}_sorted.bam $outputDir/SPAdes/bam/${j}_SPAdes.bam
done