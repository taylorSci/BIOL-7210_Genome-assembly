#!/usr/bin/python

import os

samples = os.listdir("/home/team1/genome_assembly/output/fastp/")



for each in samples:
  os.system("bwa index /home/team1/genome_assembly/output/abyss/abyss_output/contigs/abyss_" + each + "_fastp-contigs.fa")
  os.system("bwa mem /home/team1/genome_assembly/output/abyss/abyss_output/contigs/abyss_" + each + "_fastp-contigs.fa /home/team1/genome_assembly/output/fastp/" + each + "/" + each + "_1_fp.fq.gz /home/team1/genome_assembly/output/fastp/" + each + "/" + each + "_2_fp.fq.gz > " + each + ".sam")
  os.system("samtools fixmate -O bam " + each + ".sam " + each + ".bam")
  os.system("samtools sort -O bam -o " + each + "_sorted.bam -T temp " + each + ".bam") 
  os.system("cp " + each + "_sorted.bam /home/team1/genome_assembly/output/abyss/abyss_output/bam_mem/")
