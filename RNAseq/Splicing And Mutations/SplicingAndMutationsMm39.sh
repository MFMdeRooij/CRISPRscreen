#!/bin/bash
# First Download a recent GTP file from Ensembl download page, and make a file with known splice sites with the Hisat2 
# delivered python script (python hisat2_extract_splice_sites.py mm39.105.gtf > ~/MouseGenome/mm39.105_spliceSites.txt)
# The first time, make this bash script executable (chmod 755 SplicingAndMutationsMm39.sh)
# Search for RNAseq data on https://www.ncbi.nlm.nih.gov/sra/, and fill in the SRR IDs, 
# cell line ID and a P (paired-end) or S (single-end) separated by commas 
# (SRR6145108,AEC,S) in MapSamples.txt,
# and run the script (./SplicingAndMutationsMm39.sh on the command line)
# Afterwards you can view the mapping in IGV viever, or plot the splicing and mutations with the R
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
			hisat2 -x ~/MouseGenome/mm39\
			--known-splicesite-infile ~/MouseGenome/mm39.105_splicesites.txt\
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
			hisat2 -x ~/MouseGenome/mm39\
			--known-splicesite-infile ~/MouseGenome/mm39.105_splicesites.txt\
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

	samtools view -b ${s}_${c}.sort.bam "chr6:117145496-117158328" > ${s}_${c}/${s}_${c}.mcxcl12.bam	
	samtools faidx ~/MouseGenome/mm39.fa "chr6:117145496-117158328" > mcxcl12.fa
	samtools view -b ${s}_${c}.sort.bam "chr2:102641486-102731970" > ${s}_${c}/${s}_${c}.mcd44.bam
	samtools faidx ~/MouseGenome/mm39.fa "chr2:102641486-102731970" > mcd44.fa
	samtools view -b ${s}_${c}.sort.bam "chr5:142888870-142892509" > ${s}_${c}/${s}_${c}.mactb.bam
	samtools faidx ~/MouseGenome/mm39.fa "chr5:142888870-142892509" > mactb.fa

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