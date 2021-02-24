#!/bin/bash

mkdir $outputDir/SPAdes
mkdir $outputDir/SPAdes/contigs
mkdir $outputDir/SPAdes/extra

files=($(ls $outputDir/fastp/fastq_output/))

cd $toolsDir/SPAdes/SPAdes-3.15.0/bin/
mkdir SPAdes_total

for i in "${files[@]}";
do 
  ./spades.py --isolate -1 $files/$i/${i}_1_fp.fq.gz -2 $files/$i/${i}_2_fp.fq.gz -o $outputDir/SPAdes/extra/$i -t 4
  cp $outputDir/SPAdes/extra/$i/contigs.fasta $outputDir/SPAdes/contigs/${i}_SPAdes.fasta
done

