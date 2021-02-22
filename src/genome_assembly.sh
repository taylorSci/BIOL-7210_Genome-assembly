#!/bash/bin

# TODO Finalize help message
helpMessage="
USAGE
	genome_assembly [OPTIONS...] <INPUT_READS_DIRECTORY>

DESCRIPTION
A script to install and run a pipeline which assembles sets of paired-end reads in FASTQ format into  genome assemblies with multiple tools.
Script can be run with or without tool installation option (for pipeline reuse).
Preprocesses input reads, provides assembly quality metrics, and attempts to reconcile assemblies into meta-assemblies.
Sequence alignment maps are generated on assemblies to assist with QC and reconciliation.
Developed on Illumina bridge amplification dye sequencing data.

PREREQUISITES:
	git
	conda

TOOLS INSTALLED/INVOKED:
	Read quality:
		Fastp
		MultiQC
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
	-h 						display help
	-i 						install pipeline
	-o 	<OUTPUT_FOLDER> 	(defaults to sibling ('../output') of input reads directory)
	-t 	<TOOLS_FOLDER>	 	directory where pipeline tools will be installed (defaults to sibling ('../tools') of input reads directory)
	-p 						do NOT write TOOLS_FOLDER to PATH and modify startup file (eg .bash_profile) accordingly
"

# Parse optional arguments  # TODO add assembly parameters
install_=false
outputDir=""
toolsDir=""
pathFlag=false
while getopts "hi:o:t:p" option
do
	case $option in
		h) 	echo $helpMessage;;
		i) 	install_=true;;
		o) 	outputDir=$OPTARG;;
		t) 	toolsDir=$OPTARG;;
		p) 	pathFlag=true;;
		*) 	echo "UNKNOWN OPTION $OPTARG PROVIDED"
			exit;;
	esac
done

# Parse positional arguments
inputDir=${@:$OPTIND:1}

# Make required directories
if [ -z $outputDir ]; then
do
	outputDir=($(dirname $inputDir)/output)
done
mkdir outputDir
if [ -z $toolsDir ]; then
do
	toolsDir=($(dirname $inputDir)/toolsDir)
done
mkdir toolsDir

# Modify PATH variable and startup file
if install_ && ! pathFlag
then
	if ! [[ $PATH =~ $toolsDir ]]
	then
		export PATH=$PATH:$toolsDir
		for startupFile in .bash_profile .bash_login .profile
		do
			if [ -f ~/$startupFile ]
				printf "\n\nif ! [[ \$PATH =~ $toolsDir ]]\nthen\n\texport PATH=\$PATH:$toolsDir\nfi\n" >> ~/$startupFile
				break
			fi
		done
	fi
fi

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
		mkdir -p $toolsDir/gam-ngs/build
		boostPath=$CONDA_PREFIX/include/boost
		sparsehashPath=$CONDA_PREFIX/include/sparsehash
		cmake -DBOOST_ROOT=$boostPath -DBoost_NO_BOOST_CMAKE=TRUE -DSPARSEHASH_ROOT=$sparsehashPath $toolsDir/gam-ngs/
		make $toolsDir/gam-ngs/build
		ln $toolsDir/gam-ngs/bin/gam-create $toolsDir/gam-create
	fi
fi

assemblers="ABySS SKESA SPAdes"

# Preprocess reads
echo "Analyzing and trimming reads..."


# Assemble genomes & generating sequence alignment maps
echo "Assembling with ABySS..."


echo "Assembling with SKESA..."


echo "Assembling with SPAdes..."


# Reconcile assemblies
for assembler in assemblers
do
	echo "Analyzing $assembler output..."

done

echo "Reconciling assemblies..."


echo "Analyzing meta-assembly output..."

