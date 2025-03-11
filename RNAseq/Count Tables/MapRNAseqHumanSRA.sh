#!/bin/bash
# To install the Linux tools (cutadapt,seqtk,samtools,hisat2) follow the instructions in GeneExpressionRNAseqHisat2.R
# The first time, make this bash script executable (chmod 755 MapRNAseqHumanSRA.sh)
# Search for RNAseq data on https://www.ncbi.nlm.nih.gov/sra/, and fill in the SRR IDs, cell line ID and a P (paired-end) or S (single-end) separated by commas (SRR8615345,NAMALWA,P) in MapSamples.txt 
# Run the script (./MapRNAseqHumanSRA.sh) on the command line
# Afterwards you can view the mapping in IGV viever, or make a count table with the R script (GeneExpressionRNAseqHisat2.R)
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2021
# info: m.f.derooij@amsterdamumc.nl

mkdir -p pairedEnd
mkdir -p singleEnd
for line in `cat MapSamples.txt`
do 
	s="$(echo $line | cut -d ',' -f 1)"
	c="$(echo $line | cut -d ',' -f 2)"
	p="$(echo $line | cut -d ',' -f 3 | sed $'s/\r//')"
	echo ID:$s Cell:$c Reads:$p

	# Download from NCBI
	prefetch $s
	mv $PWD/${s}/${s}.sra $PWD
	rm -r $PWD/${s}
	while : 
	do
		if [ $p = P ] 
		then
			# Convert sra to fastq
			fastq-dump -v --split-files ${s}.sra

			# Remove Illumina's TruSeq adapters
			cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA\
				 -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\
				 -o ${s}_1.ca.fastq -p ${s}_2.ca.fastq\
			 	    ${s}_1.fastq ${s}_2.fastq

			# Trim low quality bases
			seqtk trimfq ${s}_1.ca.fastq > ${s}_1.fq
			seqtk trimfq ${s}_2.ca.fastq > ${s}_2.fq

			# Map to genome
			hisat2 -x ~/HumanGenome/hg38\
			--known-splicesite-infile ~/HumanGenome/hg38.105_spliceSites.txt\
			-1 ${s}_1.fq -2 ${s}_2.fq -S ${s}_${c}.sam\
			--summary-file ${s}_${c}.summmary.txt

			# Create BAM file for IGV-viewer
			samtools view -S -b ${s}_${c}.sam > ${s}_${c}.bam

			mv ${s}_${c}.sam pairedEnd/
			rm ${s}_1.fastq
			rm ${s}_2.fastq
			rm ${s}_1.ca.fastq
			rm ${s}_2.ca.fastq
			rm ${s}_1.fq 
			rm ${s}_2.fq
			break
		elif [ $p = S ] 
		then
			# Convert sra to fastq
			fastq-dump -v ${s}.sra

			# Remove Illumina's TruSeq adapters
			cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA\
				 -o ${s}.ca.fastq ${s}.fastq

			# Trim low quality bases
			seqtk trimfq ${s}.ca.fastq > ${s}.fq

			# Map to genome
			hisat2 -x ~/HumanGenome/hg38\
			--known-splicesite-infile ~/HumanGenome/hg38.105_spliceSites.txt\
			-U $s.fq -S ${s}_${c}.sam --summary-file ${s}_${c}.summmary.txt

			# Create BAM file for IGV-viewer
			samtools view -S -b ${s}_${c}.sam > ${s}_${c}.bam

			mv ${s}_${c}.sam singleEnd/
			rm ${s}.fastq
			rm ${s}.ca.fastq 
			rm ${s}.fq 
			break
		else
			echo What means $p?, Is $s paired-end \(P\) or single-end \(S\)?
			read p
		fi
	done

	# Sort and index BAM file for IGV-viewer
	samtools sort ${s}_${c}.bam > ${s}_${c}.sort.bam
	samtools index ${s}_${c}.sort.bam
	rm ${s}_${c}.bam
	rm ${s}.sra
done