#!/bin/bash
# Search for RNAseq data on https://www.ncbi.nlm.nih.gov/sra/, and fill in the SRR IDs, 
# cell line ID and a P (paired-end) or S (single-end) separated by commas 
# (SRR8615345,NAMALWA,P) in MapSamples.txt, 
# and run the script (./MapRNAseqHumanSRA.sh on the command line)
# The first time, make the bash script executable (chmod 755 MapRNAseqHumanSRA.sh)
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2019,
# info: m.f.derooij@amsterdamumc.nl
mkdir pairedEnd
mkdir singleEnd
for line in `cat MapSamples.txt`
do 
	s="$(echo $line | cut -d"," -f1)"
	c="$(echo $line | cut -d"," -f2)"
	p="$(echo $line | cut -d"," -f3)"
	echo ID:$s Cell:$c Reads:$p
	prefetch $s
	mv ~/ncbi/sra/$s.sra $PWD 
	if [ $p = P ] 
	then
		fastq-dump -v --split-files $s.sra
		seqtk trimfq ${s}_1.fastq > ${s}_1.fq
		seqtk trimfq ${s}_2.fastq > ${s}_2.fq
		hisat2 -x ~/HumanGenome/hg38\
		--known-splicesite-infile ~/HumanGenome/hg38_splicesites.txt\
	 	-1 ${s}_1.fq -2 ${s}_2.fq -S ${s}_${c}.sam\
	 	--summary-file ${s}_${c}.summmary.txt
	 	samtools view -S -b ${s}_${c}.sam > ${s}_${c}.bam
	 	mv ${s}_${c}.sam pairedEnd/${s}_${c}.sam
		rm ${s}_1.fastq
		rm ${s}_2.fastq
		rm ${s}_1.fq 
		rm ${s}_2.fq
	else
		fastq-dump -v $s.sra
		seqtk trimfq $s.fastq > $s.fq
		hisat2 -x ~/HumanGenome/hg38\
		--known-splicesite-infile ~/HumanGenome/hg38_splicesites.txt\
	 	-U $s.fq -S ${s}_${c}.sam --summary-file ${s}_${c}.summmary.txt
	 	samtools view -S -b ${s}_${c}.sam > ${s}_${c}.bam
	 	mv ${s}_${c}.sam singleEnd/${s}_${c}.sam
		rm $s.fastq 
		rm $s.fq 
	fi
	samtools sort ${s}_${c}.bam > ${s}_${c}.sort.bam
	samtools index ${s}_${c}.sort.bam
	rm ${s}_${c}.bam
done
