# Team 1 Genome Assembly

Team 1 Genome Assembly repository. This repository includes the script genome_assembly which inputs a directory of data (oraginzed in sub-directories includeing pair end reads of the format 1_fq.gz and 2_fq.gz for each isolate) and goes through the entire genome assembly pipeline.

## Part 1: Pre-assembly QC using fastp and consolidation of all analysis filed using multiqc.
	
Fastp is an all-in-one trimming and filtering tool for raw fastq files. Outputs are quality controlled reads denoted by <isolatename>_fp_1.fq.gz and <isolatename>_fp_2.fq.gz, and fastp.html and fastp.json files used for visualization of filter results.


The fastp command used was as
`fastp -i *_1.fq.gz -I *_2.fq.gz -o ${dir}_1_fp.fq.gz -O ${dir}_2_fp.fq.gz -f 5 -t 5 -5 -3 -M 28 -W 20 -e 28 -c`

Options used:
-i, -I, -o, -O indicate input and output file names for both reads, 
-f 5;	trim of the first 5 bases for both reads
-t 5;	trim of the last 5 bases for both reads
-5 -3;	cut front and cut tail sliding windows will be used
-M 28;	the mean quality cutoff for the sliding window (phred score >= 28)
-W 20;	number of bases in sliding window
-e 28;	get rid of any reads that have an everage quality score below 28
-c;	base correction


## Part2: De-novo assembly using SKESA, ABySS, and SPAdes assemblers
SKESA, SPAdes and ABySS are devono assemlers which employ the de bruijn graph method of assembling reads into contigs. This method involves breaking the reads into small fragments of size k, known as k-mers, and then traversing the Eulerian path through the k-mers to goin reads into larger contigs. SPAdes and SKESA automatucally determine and optimize whcih k-mers to use to assemble the genomes. ABySS requires input of a sinlge k-mer size which we optimized to 21 using k-mer optimizartion software, KmerGenie. Below are the usage cases for the three assemblers:
```
abyss-pe k=21 name=CGTxxxx in=‘CGTxxxx_1.fq.gz CGTxxxx_2.fq.gz’
spades.py --careful –1 <CGTxxx_1.fq.gz> -2 <CGTxxx_2.fq.gz> -o <output directory> -t <# of cores> -k <k-mer sizes>
skesa --reads CGTxxxx_1.fq.gz,CGTxxxx_2.fq.gz --cores <cores> --memory <memory> > CGTxxxx_SKESA.fasta
```
## Part 3: Post Assembly QC and Meta-Assembly

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
	-h	display help
	-i	install pipeline
	-o	[OUTPUT_FOLDER] (defaults to sibling ('../output') of input reads directory)
	-t	[TOOLS_FOLDER] Directory where pipeline tools will be installed (defaults to sibling ('../tools') of input reads directory)
"
