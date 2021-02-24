#!/usr/bin/env bash

#remove next 4 lines
c=6
m=10
inputDir=${@:$OPTIND:1}
outputDir=($(dirname $inputDir)/output)

#hard link no longer needed since skesa will be installed usign conda
#ln $inputDir/toolsDir/skesa/SKESA/skesa .run_skesa
mkdir $outputDir/assemblies/SKESA/
mkdir $outputDir/assemblies/SKESA/extra

skesa_in=($(ls $outputDir/read_QC/fastp))

for i in "${skesa_in[@]}";
do
    skesa --reads $outputDir/read_QC/fastp/$i/${i}_1_fp.fq.gz,$outputDir/read_QC/fastp/$i/${i}_2_fp.fq.gz --cores $c -- memory $m > $outputDir/assemblies/SKESA/${i}_SKESA.fasta
done