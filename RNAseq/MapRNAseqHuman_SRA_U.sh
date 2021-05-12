#!/bin/bash

for s in `cat MapFiles.txt`
do 
	echo $s
	prefetch $s
	mv ~/ncbi/sra/$s.sra $PWD 
	fastq-dump -v $s.sra
	seqtk trimfq $s.fastq > $s.fq

	hisat2 -x ~/HumanGenome/hg38\
	--known-splicesite-infile ~/HumanGenome/hg38_splicesites.txt\
 	-U $s.fq -S $s.sam --summary-file $s.summmary.txt
	samtools view -S -b $s.sam > $s.bam
	samtools sort $s.bam > $s.sort.bam
	samtools index $s.sort.bam
	rm $s.bam
	rm $s.fastq
	rm $s.fq
done
