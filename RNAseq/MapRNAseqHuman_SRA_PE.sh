#!/bin/bash

for s in `cat MapFiles.txt`
do 
	echo $s
	prefetch $s
	mv ~/ncbi/sra/$s.sra $PWD 
	fastq-dump -v --split-files $s.sra
	seqtk trimfq ${s}_1.fastq > ${s}_1.fq
	seqtk trimfq ${s}_2.fastq > ${s}_2.fq

	hisat2 -x ~/HumanGenome/hg38\
	--known-splicesite-infile ~/HumanGenome/hg38_splicesites.txt\
 	-1 ${s}_1.fq -2 ${s}_2.fq -S $s.sam --summary-file $s.summmary.txt
	samtools view -S -b $s.sam > $s.bam
	samtools sort $s.bam > $s.sort.bam
	samtools index $s.sort.bam
	rm $s.sam
	rm $s.bam
	rm ${s}_1.fastq
	rm ${s}_2.fastq
	rm ${s}_1.fq 
	rm ${s}_2.fq
done
