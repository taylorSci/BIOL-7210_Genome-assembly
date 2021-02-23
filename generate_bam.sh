#!/usr/bin/env bash

inputDir=${@:$OPTIND:1}
outputDir=($(dirname $inputDir)/output)

mkdir temp
cd temp
mkdir $outputDir/bam_files

file_prefix=($(ls $outputDir/read-QC/fastp))

for i in "${file_prefix[@]}";
do
    bwa index $outputDir/assemblies/skesa/${i}.skesa.fa
    bwa mem $outputDir/assemblies/skesa/${i}.skesa.fa $outputDir/read-QC/fastp/$i/${i}_1_fp.fq.gz $outputDir/read-QC/fastp/$i/${i}_2_fp.fq.gz > ./${i}.skesa.sam
    samtools fixmate -O bam ./${i}.skesa.sam ./${i}.skesa.bam
    samtools sort -O bam -o $outputDir/bam_files/${i}_sorted.skesa.bam -T temp ./${i}.skesa.sam

    bwa index $outputDir/assemblies/SPAdes/contigs/${i}_SPAdes.fasta
    bwa mem $outputDir/assemblies/SPAdes/contigs/${i}_SPAdes.fasta $outputDir/read-QC/fastp/$i/${i}_1_fp.fq.gz $outputDir/read-QC/fastp/$i/${i}_2_fp.fq.gz > ./${i}_SPAdes.sam
    samtools fixmate -O bam ./${i}_SPAdes.sam ./${i}_SPAdes.bam
    samtools sort -O bam -o $outputDir/bam_files/${i}_sorted.SPAdes.bam -T temp ./${i}_SPAdes.sam
done

cd ..
rm -r temp