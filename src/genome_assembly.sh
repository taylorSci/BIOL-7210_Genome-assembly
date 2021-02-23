#!/usr/bin/bash

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
	Assembly reconciliation:
		GAM-NGS

OPTIONS
	-h 										display help
	-i 										install pipeline
	-o 	<OUTPUT_FOLDER> 					(defaults to sibling ('../output') of input reads directory)
	-t 	<TOOLS_FOLDER>	 					directory where pipeline tools will be installed (defaults to sibling ('../tools') of input reads directory)
	-p 										do NOT write TOOLS_FOLDER to PATH and modify startup file (eg .bash_profile) accordingly
	-b 	<MIN_BLOCK_SIZE>			10		GAM-NGS parameter (default taken from https://doi.org/10.1186/1471-2105-14-S7-S6)
	-c 	<BLOCK_COVERAGE_THRESHOLD> 	0.75	GAM-NGS parameter (default taken from https://doi.org/10.1186/1471-2105-14-S7-S6)
"

# Parse optional arguments  # TODO add assembly parameters
install_=false
outputDir=""
toolsDir=""
pathFlag=false
Bmin=10
Tc=0.75
while getopts "hi:o:t:p" option
do
	case $option in
		h) 	echo $helpMessage;;
		i) 	install_=true;;
		o) 	outputDir=$OPTARG;;
		t) 	toolsDir=$OPTARG;;
		p) 	pathFlag=true;;
		b) 	Bmin=$OPTARG;;
		c) 	Tc=$OPTARG;;
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
mkdir $outputDir
if [ -z $toolsDir ]; then
do
	toolsDir=($(dirname $inputDir)/toolsDir)
done
mkdir $toolsDir

# Modify PATH variable and startup file
if $install_ && ! $pathFlag
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

	# TODO Add dependencies for all tools here
	echo "Installing required dependencies..."
	for dependency in gxx_linux_64 cmake zlib boost google-sparsehash bwa samtools
	do
		conda install $dependency
	done

	echo "Installing FastQC..."

	echo "Installing ABySS..."

	echo "Installing SKESA..."

	echo "Installing SPAdes..."

	echo "Installing REAPR..."

	echo "Installing dnAQET..."

	if ! (type gam-create && type gam-merge) &> /dev/null
	then
		echo "Installing GAM-NGS..."
		git clone https://github.com/vice87/gam-ngs $toolsDir/gam-ngs
		mkdir -p $toolsDir/gam-ngs/build
		boostPath=$CONDA_PREFIX/include/boost
		sparsehashPath=$CONDA_PREFIX/include/sparsehash
		cmake -DBOOST_ROOT=$boostPath -DBoost_NO_BOOST_CMAKE=TRUE -DSPARSEHASH_ROOT=$sparsehashPath $toolsDir/gam-ngs/
		make $toolsDir/gam-ngs/build
		ln $toolsDir/gam-ngs/bin/gam-create $toolsDir/gam-create
	fi
fi

isolates=$(ls $inputDir | grep -o -E "CGT[0-9]{4}")
assemblers="ABySS SKESA SPAdes"
ABySSBams=$outputDir/abyss/abyss_output/bwa_mem_corrected/\${PATTERN}_sorted.bam
ABySSContigs=$outputDir/abyss/abyss_output/corrected_fasta/abyss_\$PATTERN.fa
ABySSReapr=$outputDir/REAPR/QA_ABYSS/\${PATTERN}_QA/05.summary.report.txt
SKESABams=$outputDir/skesa/bam_mem/\${PATTERN}_sorted_mem.bam
SKESAContigs=$outputDir/skesa/skesa_outputs/\$PATTERN.skesa_contig.fa
SKESAReapr=$outputDir/REAPR/QA_SKESA/\${PATTERN}_QA/05.summary.report.txt
SPAdesBams=$outputDir/SPAdes/SPAdes_bam_mem/\${PATTERN}_sorted.bam
SPAdesContigs=$outputDir/SPAdes/SPAdes_Output/\$PATTERN/contigs.fasta
SPAdesReapr=$outputDir/REAPR/QA_SPAdes/\${PATTERN}_QA/05.summary.report.txt
reconDir=$outputDir/assembly-reconciliation/

# Preprocess reads
echo "Analyzing and trimming reads..."


# Assemble genomes & generating sequence alignment maps
echo "Assembling with ABySS..."


echo "Assembling with SKESA..."


echo "Assembling with SPAdes..."

# Assembly QC with REAPR
for assembler in assemblers
do
	echo "Analyzing $assembler output..."
	
done

# Reconcile assemblies
echo "Reconciling assemblies..."
tmpDir=$reconDir/tmp
mkdir $tmpDir

# Tests whether the first assembly is better than the second
# Compares by number of errors, then number of warnings, then N50
# numErrs1 numErrs2 numWarns1 numWarns2 n50_1 n50_2
compare_assemblies {
	if $1 < $2
	then
		echo true
	elif $1 > $2
		echo false
	else
		if $3 < $4
		then
			echo true
		elif $3 > $4
		then
			echo false
		else
			if $5 > $6
			then
				echo true
			else
				echo false
			fi
		fi
	fi
}
for PATTERN in isolates
do
	echo "Merging $PATTERN..."
	for assembler in assemblers
	do
		# Index BAMs, make BAM list files
		eval "samtools index $(eval "echo \$${assembler}Bams")"
		eval "printf $(eval "echo \$${assembler}Bams\n")" > $reconDir/${assembler}_$PATTERN.bamlist
		printf "200 800" >> $reconDir/${assembler}_$PATTERN.bamlist
		#eval "ln $(eval "echo \$${assembler}Contigs") $recondir/${assembler}_$PATTERN"

		# Extract quality metrics
		declare ${assembler}NumErrs=$(eval "grep -E '[0-9]+ errors:' $(eval "echo \$${assembler}Reapr") | grep -E -o '[0-9]+'")
		declare ${assembler}NumWarns=$(eval "grep -E '[0-9]+ warnings:' $(eval "echo \$${assembler}Reapr") | grep -E -o '[0-9]+'")
		declare ${assembler}N50=$(eval "grep -E 'N50' $(eval "echo \$${assembler}Reapr") | head -n1 | sed -E 's/N50 = ([0-9]+).*/\1/'")
	done

	# Determine assembly order
	ret1=$(compare_assemblies $ABySSNumErrs $SKESANumErrs $ABySSNumWarns $SKESANumWarns $ABySSN50 $SKESAN50)
	ret2=$(compare_assemblies $SPAdesNumErrs $SKESANumErrs $SPAdesNumWarns $SKESANumWarns $SPAdesN50 $SKESAN50)
	ret3=$(compare_assemblies $ABySSNumErrs $SPAdesNumErrs $ABySSNumWarns $SPAdesNumWarns $ABySSN50 $SPAdesN50)
	if $ret1
	then
		if $ret2
		then
			if $ret3
			then
				first=ABySS
				second=SPAdes
				third=SKESA
			else
				first=SPAdes
				second=ABySS
				third=SKESA
			fi
		else
			first=ABySS
			second=SKESA
			third=SPAdes
		fi
	else
		if $ret2
		then
			first=SPAdes
			second=SKESA
			third=ABySS
			fi
		else
			if $ret3
			then
				first=SKESA
				second=ABySS
				third=SPAdes
			else
				first=SKESA
				second=SPAdes
				third=ABySS
			fi
		fi
	fi

	# Merge first two
	eval "gam-create $(eval "echo --master-bam \$${first}Bams --slave-bams \$${second}Bams --min-block-size $Bmin --output $tmpDir/intermediate")"
	eval "gam-merge $(eval "echo --master-bam \$${first}Bams --slave-bams \$${second}Bams --blocks-file $tmpDir/intermediate.blocks --master-fasta \$${first}Contigs --slave-fasta \$${second}Contigs --output $tmpDir/intermediate")"

	# Prep partially merged assembly
	bwa index $tmpDir/intermediate.gam.fasta
	bwa mem $tmpDir/intermediate.gam.fasta $inputDir/$PATTERN/$PATTERN_1.fq.gz $inputDir/$PATTERN/$PATTERN_2.fq.gz > $tmpDir/intermediate.sam
	samtools fixmate -O bam $tmpDir/intermediate.sam $tmpDir/intermediate.bam
	samtools sort -O bam -o $tmpDir/intermediate_sorted.bam -T temp $tmpDir/intermediate.bam
	samtools index $tmpDir/intermediate_sorted.bam
	printf $tmpDir/intermediate_sorted.bam\n > $tmpDir/intermediate.bamlist
	printf "200 800" >> $tmpDir/intermediate.bamlist

	#Complete merging
	eval "gam-create $(eval "echo --master-bam $tmpDir/intermediate_sorted.bam --slave-bams \$${third}Bams --min-block-size $Bmin --output $tmpDir/intermediate")"
	eval "gam-merge $(eval "echo --master-bam $tmpDir/intermediate_sorted.bam --slave-bams \$${third}Bams --blocks-file $tmpDir/intermediate.blocks --master-fasta $tmpDir/intermediate.gam.fasta --slave-fasta \$${third}Contigs --output $outputDir/$PATTERN_meta.fa")"
	printf $PATTERN\n$first:$second:$third >> merging-order.txt
done

echo "Analyzing meta-assembly output..."

