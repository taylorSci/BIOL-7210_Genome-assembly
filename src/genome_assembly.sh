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
	-n	[NUMER_OF_CORES] 			6		Number of cores that will be used to run the pipeline
	-m	[MEMORY] 					10		Amount of memory to allocate to the pipeline (GB)
"

# Parse optional arguments  # TODO add assembly parameters
install_=false
outputDir=""
toolsDir=""
pathFlag=false
Bmin=10
Tc=0.75
cores=6
mem=10
while getopts "hio:t:pb:c:m:n:" option
do
	case $option in
		h) 	echo $helpMessage
			exit;;
		i) 	install_=true;;
		o) 	outputDir=$OPTARG;;
		t) 	toolsDir=$OPTARG;;
		p) 	pathFlag=true;;
		b) 	Bmin=$OPTARG;;
		c) 	Tc=$OPTARG;;
		n) 	cores=$OPTARG;;
		m) 	mem=$OPTARG;;
		*) 	echo "UNKNOWN OPTION $OPTARG PROVIDED"
			exit;;
	esac
done

# Parse positional arguments
inputDir=${@:$OPTIND:1}

# Make required directories
if [ -z $outputDir ]
then
	outputDir=($(dirname $inputDir)/output)
fi
mkdir -p $outputDir
if [ -z $toolsDir ]
then
	toolsDir=($(dirname $inputDir)/toolsDir)
fi
mkdir -p $toolsDir

# Modify PATH variable and startup file
if $install_ && ! $pathFlag
then
	if ! [[ $PATH =~ $toolsDir ]]
	then
		export PATH=$PATH:$toolsDir
		for startupFile in .bash_profile .bash_login .profile
		do
			if [ -f ~/$startupFile ]
			then
				printf "\n\nif ! [[ \$PATH =~ $toolsDir ]]\nthen\n\texport PATH=\$PATH:$toolsDir\nfi\n" >> ~/$startupFile
				break
			fi
		done
	fi
fi

# Install tools
if $install_
then
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

	# Install ABySS and its dependencies
	echo "Installing ABySS..."
	conda install -c bioconda abyss

	# Install SKESA and its dependencies
	echo "Installing SKESA..."
	conda install -c bioconda skesa

	echo "Installing SPAdes..."
	conda install -c bioconda spades

	echo "Installing REAPR..."
	if ! type reapr &> /dev/null
	then
		cd $toolsDir
		mkdir REAPR
		cd REAPR
		wget ftp://ftp.sanger.ac.uk/pub/resources/software/reapr/Reapr_1.0.18.tar.gz
		tar -xzf Reapr_1.0.18.tar.gz -C $toolsDir/REAPR
		./configure
		make
		make install
		ln $toolsDir/REAPR/src/reapr.pl $toolsDir reapr
	fi

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

isolates=$(ls $fastpDir)
assemblies="ABySS SKESA SPAdes"

fastpDir=$outputDir/read_QC/fastp
contigPath=$outputDir/assemblies/\$assembly/contigs/\${PATTERN}_$assembly.fasta
bamPath=$outputDir/bam_files/\${PATTERN}_sorted.\$assembly.bam
qaPath=$outputDir/REAPR/unreconciled/QA_\$assembly/\${PATTERN}_QA/05.summary.report.txt
reconDir=$outputDir/assembly-reconciliation

ABySSDir=$outputDir/assemblies/ABySS
SKESADir=$outputDir/assemblies/SKESA
SPAdesDir=$outputDir/assemblies/SPAdes

ABySSContigs=$ABySSDir/contigs/\${PATTERN}_ABySS.fasta
SKESAContigs=$SKESADir/contigs/\${PATTERN}_SKESA.fasta
SPAdesContigs=$SPAdesDir/contigs/\${PATTERN}_SPAdes.fasta

# Preprocess reads
echo "Analyzing and trimming reads..."


# Assemble genomes & generating sequence alignment maps
echo "Assembling with ABySS..."
mkdir -p $ABySSDir/contigs
mkdir -p $ABySSDir/extra

for i in $isolates;
do 
	abyss-pe k=70 in="$fastpDir/$i/${i}_1_fp.fq.gz $fastpDir/$i/${i}_2_fp.fq.gz" name=$ABySSDir/extra/$i
	# Repair contig names to something REAPR likes
	reapr facheck $ABySSDir/extra/$PATTERN/${PATTERN}-contigs.fa $ABySSDir/contigs/${PATTERN}_ABySS.fasta
done

echo "Assembling with SKESA..."
mkdir -p $SKESADir/contigs
mkdir -p $SKESADir/extra

for i in $isolates;
do
	skesa --reads $outputDir/read_QC/fastp/$i/${i}_1_fp.fq.gz,$outputDir/read_QC/fastp/$i/${i}_2_fp.fq.gz --cores $cores -- memory $mem > $outputDir/assemblies/SKESA/${i}_SKESA.fasta
	ln $SKESADir/extra/$i/${i}-contigs.fa $SKESADir/contigs/${i}_SKESA.fasta
done

echo "Assembling with SPAdes..."
mkdir -p $SPAdesDir/contigs
mkdir -p $SPAdesDir/extra

for i in $isolates;
do 
	spades.py -1 $fastpDir/$i/${i}_1_fp.fq.gz -2 $fastpDir/$i/${i}_2_fp.fq.gz -o $SPAdesDir/extra/$i -t 4
	ln $SPAdesDir/extra/$i/contigs.fasta $SPAdesDir/contigs/${i}_SPAdes.fasta
done

# Post-assembly QC
echo "Analyzing $assembler output..."

for PATTERN in $isolates;
do
	for assembly in $assemblies
	do
		# Generate and index BAM files
		eval "bwa index $(echo $contigPath)"
		eval "bwa mem $(echo "$contigPath") $fastpDir/$PATTERN/${PATTERN}_1_fp.fq.gz $fastpDir/$PATTERN/${PATTERN}_2.fq.gz" > $outputDir/bam_files/${PATTERN}.$assembly.sam
		samtools fixmate -O bam $outputDir/bam_files/$PATTERN.$assembly.sam $outputDir/bam_files/$PATTERN.$assembly.bam
		eval "samtools sort -O bam -o $(echo $bamPath) -T temp $outputDir/bam_files/$PATTERN.$assembly.bam"
		eval "samtools index $(echo $bamPath)"

		# REAPR QA
		eval "reapr pipeline $(echo $contigPath $bamPath $(dirname $qaPath))"
	done

done

# Reconcile assemblies
echo "Reconciling assemblies..."
tmpDir=$reconDir/tmp
mkdir -p $tmpDir
rm -f $reconDir/merging-order.txt

# Tests whether the first assembly is better than the second
# Compares by number of errors, then number of warnings, then N50
# numErrs1 numErrs2 numWarns1 numWarns2 n50_1 n50_2
compare_assemblies() {
	if [[ $1 < $2 ]]
	then
		echo true
	elif [[ $1 > $2 ]]
	then
		echo false
	else
		if [[ $3 < $4 ]]
		then
			echo true
		elif [[ $3 > $4 ]]
		then
			echo false
		else
			if [[ $5 > $6 ]]
			then
				echo true
			else
				echo false
			fi
		fi
	fi
}
for PATTERN in $isolates
do
	echo "Merging $PATTERN..."
	for assembler in ABySS SKESA SPAdes
	do
		# Make BAM list files
		eval "echo $(echo $bamPath)" > $reconDir/${assembler}_$PATTERN.bamlist
		echo "200 800" >> $reconDir/${assembler}_$PATTERN.bamlist

		# Extract quality metrics
		declare ${assembler}NumErrs=$(eval "grep -E '[0-9]+ errors:' $(echo $qaPath) | grep -E -o '[0-9]+'")
		declare ${assembler}NumWarns=$(eval "grep -E '[0-9]+ warnings:' $(echo $qaPath) | grep -E -o '[0-9]+'")
		declare ${assembler}N50=$(eval "grep -E 'N50' $(eval echo $qaPath) | head -n1 | sed -E 's/N50 = ([0-9]+).*/\1/'")
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
	echo "$PATTERN: Performing initial merge..."
	gam-create --master-bam $reconDir/${first}_$PATTERN.bamlist --slave-bam $reconDir/${second}_$PATTERN.bamlist --min-block-size $Bmin --output $tmpDir/intermediate
	eval "gam-merge $(eval "echo --master-bam $reconDir/${first}_$PATTERN.bamlist --slave-bam $reconDir/${second}_$PATTERN.bamlist --blocks-file $tmpDir/intermediate.blocks --master-fasta \$${first}Contigs --slave-fasta \$${second}Contigs --output $tmpDir/intermediate")"

	# Prep partially merged assembly
	echo "$PATTERN: Prepping second merge..."
	bwa index $tmpDir/intermediate.gam.fasta
	bwa mem $tmpDir/intermediate.gam.fasta $fastpDir/$PATTERN/${PATTERN}_1_fp.fq.gz $fastpDir/$PATTERN/${PATTERN}_2.fq.gz > $tmpDir/intermediate.sam
	samtools fixmate -O bam $tmpDir/intermediate.sam $tmpDir/intermediate.bam
	samtools sort -O bam -o $tmpDir/intermediate_sorted.bam -T temp $tmpDir/intermediate.bam
	samtools index $tmpDir/intermediate_sorted.bam
	echo $tmpDir/intermediate_sorted.bam > $tmpDir/intermediate.bamlist
	echo "200 800" >> $tmpDir/intermediate.bamlist

	#Complete merging
	echo "$PATTERN: Performing second merge..."
	gam-create --master-bam $tmpDir/intermediate.bamlist --slave-bam $reconDir/${third}_$PATTERN.bamlist --min-block-size $Bmin --output $tmpDir/intermediate
	eval "gam-merge $(eval "echo --master-bam $tmpDir/intermediate.bamlist --slave-bam $reconDir/${third}_$PATTERN.bamlist --blocks-file $tmpDir/intermediate.blocks --master-fasta $tmpDir/intermediate.gam.fasta --slave-fasta \$${third}Contigs --output $reconDir/${PATTERN}_meta")"
	echo -e "$PATTERN\t$first:$second:$third" >> $reconDir/merging-order.txt

	echo "Analyzing meta-assembly output..."
	# Align BAMs for meta-QC
	bwa index $reconDir/${PATTERN}_meta.gam.fasta
	bwa mem $reconDir/${PATTERN}_meta.gam.fasta $fastpDir/$PATTERN/${PATTERN}_1_fp.fq.gz $inputDir/$PATTERN/${PATTERN}_2_fp.fq.gz > $reconDir/${PATTERN}_meta.sam
	samtools fixmate -O bam $reconDir/${PATTERN}_meta.sam $reconDir/${PATTERN}_meta.bam
	samtools sort -O bam -o $reconDir/${PATTERN}_meta_sorted.bam -T temp $reconDir/${PATTERN}_meta.bam

	# QC of meta-assembly
	reapr pipeline $reconDir/${PATTERN}_meta.bam $reconDir/${PATTERN}_meta.gam.fasta $outputDir/REAPR/reconciled/$PATTERN
done

# TODO Provide folder with clean output
