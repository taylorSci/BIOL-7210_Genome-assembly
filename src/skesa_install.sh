#!/usr/bin/env bash

#skesa installation
cd $toolsDir
mkdir skesa
cd skesa
git clone https://github.com/ncbi/SKESA.git
cd SKESA
conda install -c conda-forge/label/gcc7 boost
make -f Makefile.nongs

#bwa and samtools installation for generating bam files
conda install -c bioconda bwa
conda install -c bioconda samtools

