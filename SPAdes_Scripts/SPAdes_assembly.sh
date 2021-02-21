#!/bin/bash

mkdir $outputDir/SPAdes
mkdir $outputDir/SPAdes/contigs

files=($(ls $outputDir/fastp/fastq_output/))

cd $toolsDir/SPAdes/SPAdes-3.15.0/bin/
mkdir SPAdes_total

for i in "${files[@]}";
do 
  ./spades.py --isolate -1 $files/$i/${i}_1_fp.fq.gz -2 $files/$i/${i}_2_fp.fq.gz -o SPAdes_total/$i -t 4
  cp SPAdes_total/$i/contigs.fasta $outputDir/SPAdes/contigs/${i}_SPAdes.fasta
done

