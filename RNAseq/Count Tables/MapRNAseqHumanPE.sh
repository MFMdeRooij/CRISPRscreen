#!/bin/bash
# Paired-end RNAseq mapping from fastq files:
# To install the Linux tools (cutadapt,seqtk,samtools,hisat2) follow the instructions in GeneExpressionRNAseqHisat2.R
# The first time, make this bash script executable (chmod 755 MapRNAseqHumanPE.sh)
# List samples in MapSamples.txt 
# Run the script (./MapRNAseqHumanPE.sh) on the command line
# Afterwards you can view the mapping in IGV viever, or make a count table with the R script (GeneExpressionRNAseqHisat2.R)
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2021
# info: m.f.derooij@amsterdamumc.nl
for s in `cat MapSamples.txt`
do 
	echo $s

	# Remove Illumina's TruSeq adapters
	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA\
		 -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\
		 -o ${s}_1.ca.fastq -p ${s}_2.ca.fastq\
		    ${s}_1.fastq.gz ${s}_2.fastq.gz

	# Trim low quality bases
	seqtk trimfq ${s}_1.ca.fastq > ${s}_1.fq
	seqtk trimfq ${s}_2.ca.fastq > ${s}_2.fq

	# Map to genome
	hisat2 -x ~/HumanGenome/hg38\
	--known-splicesite-infile ~/HumanGenome/hg38.105_spliceSites.txt\
 	-1 ${s}_1.fq -2 ${s}_2.fq -S ${s}.sam --summary-file ${s}.summmary.txt
 	
	# Create, sort, and index BAM file for IGV-viewer
	samtools view -S -b ${s}.sam > ${s}.bam
	samtools sort ${s}.bam > ${s}.sort.bam
	samtools index ${s}.sort.bam
	rm ${s}.bam
	
	#rm ${s}_1.fastq.gz
	#rm ${s}_2.fastq.gz
	rm ${s}_1.ca.fastq
	rm ${s}_2.ca.fastq
	rm ${s}_1.fq
	rm ${s}_2.fq
done
