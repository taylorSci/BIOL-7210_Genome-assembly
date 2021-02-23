#!/bash/bin

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

	echo "Installing ABySS..."

	echo "Installing SKESA..."

	echo "Installing SPAdes..."

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
#TAYLOR, I AM MAKING THE HARD LINK TO SKESA HIDDEN. PLEASE FEEL FREE TO CHANGE THIS IF IT SHOULD BE UNHIDDEN.
ln $inputDir/$toolsDir/skesa/SKESA/skesa .run_skesa
mkdir $outputDir/assemblies/skesa/
skesa_in=($(ls $outputDir/read-QC/fastp))

for i in "${skesa_in[@]}";
do
    ./.run_skesa --reads $outputDir/read-QC/fastp/$i/${i}_1_fp.fq.gz,$outputDir/read-QC/fastp/$i/${i}_2_fp.fq.gz --cores $cores -- memory $mem > $outputDir/assemblies/skesa/$i.skesa.fa
done
echo "Contigs generated for SKESA"

echo "Assembling with SPAdes..."

echo "Generating BAM files from ABySS, SKESA, and SPAdes contigs"
# TAYLOR, I'VE IMPLEMENTED THE BAM FILE GENERATION FOR SPADES AND SKESA AND HAVE TESTED IT. I NEED TO TALK TO YOU ABOUT HOW TO DO THIS FOR ABYSS (DEADLY WHITESPACES LOL)
# I'VE ALSO EDITED THE SCRIPT LAYOUT WORD DOC SLIGHTLY TO INCLUDE A "BAM_FILES" DIRECTORY IN OUTPUTS WHERE YOU CAN GRAB EVERYTHING FROM
# ALSO, I'VE ADDED -C AND -M FLAGS FOR CORES AND MEMORY INPUTS
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


# Reconcile assemblies
for assembler in assemblers
do
	echo "Analyzing $assembler output..."

done

echo "Reconciling assemblies..."


echo "Analyzing meta-assembly output..."

