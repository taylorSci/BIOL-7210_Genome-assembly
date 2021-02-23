#!/usr/bin/env bash

# skesa installation
conda install -c bioconda skesa

# bwa and samtools installation for generating bam files
conda install -c bioconda bwa
conda install -c bioconda samtools

# all dependencies are installed alongside the tools thanks to the magic of conda
