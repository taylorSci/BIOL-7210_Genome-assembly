#!/bash/bin

mkdir $outputDir/assemblies/ABySS/contigs_reapr
ls $outputDir/assemblies/ABYSS/*.fasta > $outputDir/assemblies/ABySS/contigs_reapr/contigs.txt
cd $outputDir/assemblies/ABySS/contigs_reapr

cat contigs.txt | while read line
do reapr facheck $line $line_repaired_ABySS.fasta
done
