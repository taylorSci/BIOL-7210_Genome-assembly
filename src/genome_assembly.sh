#!/usr/bin/env bash

# TODO Finalize help message
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
	-c	[NUMER_OF_CORES] Number of cores that will be used to run the pipeline
	-m	[MEMORY] Amount of memory to allocate to the pipeline
"

# Parse user command line arguments  # TODO add assembly parameters
install_=false
outputDir=""
toolsDir=""
while getopts "hi:o:t:c:m:" option
do
	case $option in
		h) 	echo $helpMessage;;
		i) 	install_=true;;
		o) 	outputDir=$OPTARG;;
		t) 	toolsDir=$OPTARG;;
		c) 	cores=$OPTARG;;
		m) 	mem=$OPTARG;;
		*) 	echo "UNKNOWN OPTION $OPTARG PROVIDED"
			exit;;
	esac
done
inputDir=${@:$OPTIND:1}
if [ -z $outputDir ]; then
do
	outputDir=($(dirname $inputDir)/output)
done
if [ -z $toolsDir ]; then
do
	toolsDir=($(dirname $inputDir)/toolsDir)
done

#Install tools
if $install_
	if ! (type git && type conda) &> /dev/null
	then
		echo "git and conda are required to run installations."
		exit
	fi

	echo "Installing FastQC..."

	# Install ABySS and its dependencies
	echo "Installing ABySS..."
	conda install -c bioconda abyss

	# Install SKESA and its dependencies
	echo "Installing SKESA..."
	conda install -c bioconda skesa

	echo "Installing SPAdes..."
	
	# Install bwa and its dependencies
	echo "Installing bwa..."
	conda install -c bioconda bwa

	# Install samtools and its dependencies
	echo "Installing samtools"
	conda install -c bioconda samtools

	echo "Installing REAPR..."

	echo "Installing dnAQET..."

	if ! (type gam-create && type gam-merge) &> /dev/null
	then
		echo "Installing GAM-NGS..."
		for dependency in gxx_linux_64 cmake zlib boost google-sparsehash bwa samtools
		do
			conda install $dependency
		done
		git clone https://github.com/vice87/gam-ngs $toolsDir/gam-ngs
		mkdir $toolsDir/gam-ngs/build
		boostPath=$CONDA_PREFIX/include/boost
		sparsehashPath=$CONDA_PREFIX/include/sparsehash
		cmake -DBOOST_ROOT=$boostPath -DBoost_NO_BOOST_CMAKE=TRUE -DSPARSEHASH_ROOT=$sparsehashPath $toolsDir/gam-ngs/
		make $toolsDir/gam-ngs/build
	fi
fi

assemblers="ABySS SKESA SPAdes"

# Preprocess reads
echo "Analyzing and trimming reads..."


# Assemble genomes
echo "Assembling with ABySS..."


echo "Assembling with SKESA..."
mkdir $outputDir/assemblies/SKESA/
mkdir $outputDir/assemblies/SKESA/extra
skesa_in=($(ls $outputDir/read_QC/fastp))
for i in "${skesa_in[@]}";
do
    skesa --reads $outputDir/read_QC/fastp/$i/${i}_1_fp.fq.gz,$outputDir/read_QC/fastp/$i/${i}_2_fp.fq.gz --cores $c -- memory $m > $outputDir/assemblies/SKESA/${i}_SKESA.fasta
done
echo "Contigs generated for SKESA"

echo "Assembling with SPAdes..."

# generate BAM files for post Assembly QC
echo "Generating BAM files from ABySS, SKESA, and SPAdes contigs"
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

# Reconcile assemblies
for assembler in assemblers
do
	echo "Analyzing $assembler output..."

done

echo "Reconciling assemblies..."


echo "Analyzing meta-assembly output..."

