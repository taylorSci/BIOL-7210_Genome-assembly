#!/usr/bin/env bash

#will need to change this path to reflect the environment of the server
if [ ! -d skesa_output ]
then
    mkdir skesa_output
fi

#change this to PATH TO fastp_output for final version
cd $HOME/class/7210/genome_assembly/assemblies/fastp_output

ls *_1_* > reads1.txt
ls *_2_* > reads2.txt
paste -d, reads1.txt reads2.txt > reads.csv
sed 's/_1_fp.fq.gz/.skesa.fa/g' reads1.txt > out.txt

cat reads.csv | while read pairs
do
    skesa --reads $pairs --cores 6 --memory 10 > ../skesa_output/$pairs
done
rm *.txt
rm *.csv

#change this to PATH TO skesa_output for final version
cd $HOME/class/7210/genome_assembly/assemblies/skesa_output
for file in *
do
    mv "$file" "${file/_1_fp.fq.gz*/.skesa.fa}"
done
#hello
