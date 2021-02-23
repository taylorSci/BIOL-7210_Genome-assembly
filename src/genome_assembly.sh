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
"

# Parse user command line arguments  # TODO add assembly parameters
install_=false
outputDir=""
toolsDir=""
while getopts "hi:o:t:" option
do
	case $option in
		h) 	echo $helpMessage;;
		i) 	install_=true;;
		o) 	outputDir=$OPTARG;;
		t) 	toolsDir=$OPTARG;;
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
	conda install -c bioconda spades

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


echo "Assembling with SPAdes..."
mkdir $outputDir/assemblies/SPAdes
mkdir $outputDir/assemblies/SPAdes/contigs
mkdir $outputDir/assemblies/SPAdes/extra

files=($(ls $outputDir/read_QC/fastp/))


for i in "${files[@]}";
do 
  spades.py -1 $outputDir/read_QC/fastp/$i/${i}_1_fp.fq.gz -2 $outputDir/read_QC/fastp/$i/${i}_2_fp.fq.gz -o $outputDir/assemblies/SPAdes/extra/$i -t 4
  cp $outputDir/assemblies/SPAdes/extra/$i/contigs.fasta $outputDir/assemblies/SPAdes/contigs/${i}_SPAdes.fasta
done

# Reconcile assemblies
for assembler in assemblers
do
	echo "Analyzing $assembler output..."

done

echo "Reconciling assemblies..."


echo "Analyzing meta-assembly output..."

