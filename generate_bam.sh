#!/usr/bin/env bash

inputDir=${@:$OPTIND:1}
outputDir=($(dirname $inputDir)/output)

mkdir temp
cd temp
mkdir $outputDir/bam_files

file_prefix=($(ls $outputDir/read_QC/fastp))

for i in "${file_prefix[@]}";
do
    bwa index $outputDir/assemblies/SKESA/${i}_SKESA.fasta
    bwa mem $outputDir/assemblies/SKESA/${i}_SKESA.fasta $outputDir/read_QC/fastp/$i/${i}_1_fp.fq.gz $outputDir/read_QC/fastp/$i/${i}_2_fp.fq.gz > ./${i}_SKESA.sam
    samtools fixmate -O bam ./${i}_SKESA.sam ./${i}_SKESA.bam
    samtools sort -O bam -o $outputDir/bam_files/${i}_sorted.SKESA.bam -T temp ./${i}_SKESA.sam
    mv $outputDir/assemblies/SKESA/*.amb $outputDir/assemblies/SKESA/extra
    mv $outputDir/assemblies/SKESA/*.ann $outputDir/assemblies/SKESA/extra
    mv $outputDir/assemblies/SKESA/*.bwt $outputDir/assemblies/SKESA/extra
    mv $outputDir/assemblies/SKESA/*.pac $outputDir/assemblies/SKESA/extra
    mv $outputDir/assemblies/SKESA/*.sa $outputDir/assemblies/SKESA/extra

    bwa index $outputDir/assemblies/SPAdes/contigs/${i}_SPAdes.fasta
    bwa mem $outputDir/assemblies/SPAdes/contigs/${i}_SPAdes.fasta $outputDir/read_QC/fastp/$i/${i}_1_fp.fq.gz $outputDir/read_QC/fastp/$i/${i}_2_fp.fq.gz > ./${i}_SPAdes.sam
    samtools fixmate -O bam ./${i}_SPAdes.sam ./${i}_SPAdes.bam
    samtools sort -O bam -o $outputDir/bam_files/${i}_sorted.SPAdes.bam -T temp ./${i}_SPAdes.sam
    mv $outputDir/assemblies/SPAdes/contigs/*.amb $outputDir/assemblies/SPAdes/extra
    mv $outputDir/assemblies/SPAdes/contigs/*.ann $outputDir/assemblies/SPAdes/extra
    mv $outputDir/assemblies/SPAdes/contigs/*.bwt $outputDir/assemblies/SPAdes/extra
    mv $outputDir/assemblies/SPAdes/contigs/*.pac $outputDir/assemblies/SPAdes/extra
    mv $outputDir/assemblies/SPAdes/contigs/*.sa $outputDir/assemblies/SPAdes/extra

# NEED AEKANSH'S REAPR CONTIG FILE REPAIR CODE HERE

    bwa index $outputDir/assemblies/ABySS/contigs_reapr/${i}_repaired_ABySS.fasta
    bwa mem $outputDir/assemblies/ABySS/contigs_reapr/${i}_repaired_ABySS.fasta $outputDir/read_QC/fastp/$i/${i}_1_fp.fq.gz $outputDir/read_QC/fastp/$i/${i}_2_fp.fq.gz > ./${i}_ABySS.sam
    samtools fixmate -O bam ./${i}_ABySS.sam ./${i}_ABySS.bam
    samtools sort -O bam -o $outputDir/bam_files/${i}_sorted.ABySS.bam -T temp ./${i}_ABySS.sam
    mv $outputDir/assemblies/ABySS/contigs_reapr/*.amb $outputDir/assemblies/ABySS/extra
    mv $outputDir/assemblies/ABySS/contigs_reapr/*.ann $outputDir/assemblies/ABySS/extra
    mv $outputDir/assemblies/ABySS/contigs_reapr/*.bwt $outputDir/assemblies/ABySS/extra
    mv $outputDir/assemblies/ABySS/contigs_reapr/*.pac $outputDir/assemblies/ABySS/extra
    mv $outputDir/assemblies/ABySS/contigs_reapr/*.sa $outputDir/assemblies/ABySS/extra

done
cd ..
rm -r temp