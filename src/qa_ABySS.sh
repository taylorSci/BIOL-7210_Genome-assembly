#!/bash/bin
echo "Analyzing ABYSS Output..."
mkdir $outputDir/QA_ABySS

ls $outputDir/assemblies/ABySS/contigs_reapr/*_repaired_ABySS.fasta > $outputDir/QA_ABySS/contigs.txt
ls $outputDir/bam_files/*_ABySS.bam > $outputDir/QA_ABySS/bam.txt
cd $outputDir/bam_files
ls *_ABySS.bam > $outputDir/QA_ABySS/temp.txt
cd $outputDir/QA_ABYSS/
cat temp.txt | tr -d "ABySS.bam" > output.txt
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
