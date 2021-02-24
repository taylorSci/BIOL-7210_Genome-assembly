#!/bash/bin
echo "Analyzing SPAdes Output..."
mkdir $outputDir/QA_SPAdes

ls $outputDir/assemblies/SPAdes/contigs/*.fasta > $outputDir/QA_SPAdes/contigs.txt
ls $outputDir/bam_files/*_SPAdes.bam > $outputDir/QA_SPAdes/bam.txt

cd $outputDir/bam_files
ls *_SPAdes.bam > $outputDir/QA_SPAdes/temp.txt
cd $outputDir/QA_SPAdes/
cat temp.txt | tr -d "_SPAdes.bam" > output.txt
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
