#!/bash/bin
echo "Analyzing SKESA Output..."
mkdir $outputDir/QA_SKESA

ls $outputDir/assemblies/SKESA/*.fasta > $outputDir/QA_SKESA/contigs.txt
ls $outputDir/bam_files/*_SKESA.bam > $outputDir/QA_SKESA/bam.txt

cd $outputDir/bam_files
ls *_SKESA.bam > $outputDir/QA_SKESA/temp.txt
cd $outputDir/QA_SKESA/
cat temp.txt | tr -d "_sorted.SKESA.bam" > output.txt
rm temp.txt

cat contigs.txt | while read line
do reapr facheck $line
done

paste -d, contigs.txt bam.txt output.txt > temp.csv
cat temp.csv | tr "," " " > contig_bam_output.csv
rm temp.csv

cat contig_bam_output.csv | while read line
do reapr pipeline $line
done
