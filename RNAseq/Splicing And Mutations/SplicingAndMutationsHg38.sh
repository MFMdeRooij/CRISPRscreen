#!/bin/bash
# First Download a recent GTP file from Ensembl download page, and make a file with known splice sites with the Hisat2 
# delivered python script (python hisat2_extract_splice_sites.py hg38.105.gtf > ~/HumanGenome/hg38.105_spliceSites.txt)
# The first time, make this bash script executable (chmod 755 SplicingAndMutationsHg38.sh)
# Search for RNAseq data on https://www.ncbi.nlm.nih.gov/sra/, and fill in the SRR IDs, 
# cell line ID and a P (paired-end) or S (single-end) separated by commas 
# (SRR8615282,HCT116,P) in MapSamples.txt,
# and run the script (./SplicingAndMutationsHg38.sh on the command line)
# Afterwards you can view the mapping in IGV viever, or plot the splicing  and mutationswith the R
# script (SplicingAndMutations.R)
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2022,
# info: m.f.derooij@amsterdamumc.nl

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
			
			rm ${s}_${c}.sam
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
			
			rm ${s}_${c}.sam
			rm $s.fastq
			rm $s.ca.fastq
			rm $s.fq
			break
		else
			echo What means $p?, Is $s paired-end \(P\) or single-end \(S\)?
			read p
		fi
	done
	samtools sort ${s}_${c}.bam > ${s}_${c}.sort.bam
	samtools index ${s}_${c}.sort.bam
	rm ${s}_${c}.bam
	rm $s.sra
	
	mkdir ${s}_${c}
	mv ${s}_${c}.summmary.txt ${s}_${c}

	samtools view -b ${s}_${c}.sort.bam "chr10:44370165-44385092" > ${s}_${c}/${s}_${c}.hcxcl12.bam	
	samtools faidx ~/HumanGenome/hg38.fa "chr10:44370165-44385092" > hcxcl12.fa
	samtools view -b ${s}_${c}.sort.bam "chr11:35139171-35232402" > ${s}_${c}/${s}_${c}.hcd44.bam
	samtools faidx ~/HumanGenome/hg38.fa "chr11:35139171-35232402" > hcd44.fa
	samtools view -b ${s}_${c}.sort.bam "chr7:5527147-5530601" > ${s}_${c}/${s}_${c}.hactb.bam
	samtools faidx ~/HumanGenome/hg38.fa "chr7:5527147-5530601" > hactb.fa

	cd ${s}_${c}
	for bamfile in *.bam; do
		samtools view -H $bamfile > temp.sam
		samtools view $bamfile | awk '($6 ~ /N/)' >> temp.sam
		samtools view -bS temp.sam > splice_$bamfile
		samtools index $bamfile
		samtools index splice_$bamfile
	done
	rm temp.sam

	cd ..
	rm ${s}_${c}.sort.bam
	rm ${s}_${c}.sort.bam.bai
done