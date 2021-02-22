#!/usr/bin/env bash

#replace the pathing with $outputDir for final script. I'm keeping the pathing here as I test things locally 
cd ~/class/7210/genome_assembly/
if [ ! -d skesa ]
then
    mkdir ~/class/7210/genome_assembly/skesa
    mkdir ~/class/7210/genome_assembly/skesa/contigs
    mkdir ~/class/7210/genome_assembly/skesa/bam_files
fi

files=($(ls $HOME/class/7210/genome_assembly/fastp/))
cd ~/class/7210/genome_assembly/skesa

for i in "${files[@]}";
do
# be sure to replace --cores with --$cores and --memory with --$memory
    skesa --reads $HOME/class/7210/genome_assembly/fastp/$i/${i}_1_fp.fq.gz,$HOME/class/7210/genome_assembly/fastp/$i/${i}_2_fp.fq.gz --cores 6 -- memory 10 > $HOME/class/7210/genome_assembly/skesa/contigs/$i.skesa.fa
done


