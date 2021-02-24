#!/bin/bash

helpMessage="
USAGE
	genome_assembly [OPTIONS...] INPUT_READS_DIRECTORY

DESCRIPTION
A script to install and run a pipeline which assembles sets of paired-end reads in FASTQ format into  genome assemblies with multiple tools.
Script can be run with or without tool installation option (for pipeline reuse).
Preprocesses input reads, provides assembly quality metrics, and attempts to reconcile assemblies into meta-assemblies.
Developed on Illumina bridge amplification dye sequencing data.

PREREQUISITES:
	git
	conda

TOOLS INSTALLED/INVOKED:
	Read quality:
		FastQC
	Genome assembly:
		ABySS
		SKESA
		SPAdes
	Assembly quality control:
		REAPR
		dnAQET
	Assembly reconciliation:
		GAM-NGS

OPTIONS
	-h 	display help
	-i 	install pipeline
	-o 	[OUTPUT_FOLDER] (defaults to sibling ('../output') of input reads directory)
	-t 	[TOOLS_FOLDER] Directory where pipeline tools will be installed (defaults to sibling ('../tools') of input reads directory)
"



install_=false
outputDir=""
toolsDir=""


while getopts "hio:t:" option
do
	case $option in
		h) echo $helpMessage;;
		i) install_=true
			echo "Installing tools";;
		o) outputDir=$OPTARG;;
		t) toolsDir=$OPTARG;;
		*) echo "UNKNOWN OPTION $OPTARG PROVIDED"
			exit;;
	esac
done

inputDir=${@:$OPTIND:1}
if [ -z $outputDir ]; then
	outputDir=($(dirname $inputDir)/output)
fi
mkdir $outputDir
if [ -z $toolsDir ]; then
	toolsDir=($(dirname $inputDir)/toolsDir)
fi
mkdir $toolsDir

if [ $install_ == true ]
then
	echo 'Installing fastp ...'
	mkdir -p $toolsDir/fastp && cd $_
	conda install -c bioconda fastp
	echo 'Installing multiqc...'
	mkdir -p $toolsDir/fastp && cd $_
	conda install -c bioconda -c conda-forge multiqc

fi
mkdir $outputDir/fastp
cd $inputDir
for dir in *
	do
		mkdir ${outputDir}/fastp/${dir}_fp
		cd $(dirname $inputDir)
		fastp -i $inputDir/${dir}/${dir}_1.fq.gz -I $inputDir/${dir}/${dir}_2.fq.gz -o ${dir}_1_fp.fq.gz -O ${dir}_2_fp.fq.gz -f 5 -t 5 -5 -3 -M 28 -W 20 -e 28 -c
		mv *fp* $outputDir/fastp/${dir}_fp
		mv fastp* $outputDir/fastp/${dir}_fp
		cd $inputDir
done

cd $(dirname $inputDir)
cd ${outputDir}/fastp
multiqc .



