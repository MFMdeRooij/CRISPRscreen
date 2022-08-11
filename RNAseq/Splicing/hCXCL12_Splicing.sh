#!/bin/bash
# First Download a recent GTP file from Ensembl download page, and make a file with known splice sites with the Hisat2 
# delivered python script (python hisat2_extract_splice_sites.py hg38.105.gtf > ~/HumanGenome/hg38.105_spliceSites.txt)
# The first time, make this bash script executable (chmod 755 hCXCL12_Splicing.sh)
# Search for RNAseq data on https://www.ncbi.nlm.nih.gov/sra/, and fill in the SRR IDs, 
# cell line ID and a P (paired-end) or S (single-end) separated by commas 
# (SRR8615345,NAMALWA,P) in MapSamples.txt,
# and run the script (./hCXCL12_Splicing.sh on the command line)
# Afterwards you can view the mapping in IGV viever, or plot the splicing with the R
# script (Splicing.R)
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2022,
# info: m.f.derooij@amsterdamumc.nl

for line in `cat MapSamples.txt`
do 
	s="$(echo $line | cut -d"," -f1)"
	c="$(echo $line | cut -d"," -f2)"
	p="$(echo $line | cut -d"," -f3)"
	echo ID:$s Cell:$c Reads:$p
	prefetch $s
	mv ~/ncbi/sra/$s.sra $PWD 
	while : 
	do
		if [ $p = P ] 
		then
			fastq-dump -v --split-files $s.sra
			seqtk trimfq ${s}_1.fastq > ${s}_1.fq
			seqtk trimfq ${s}_2.fastq > ${s}_2.fq
			hisat2 -x ~/HumanGenome/hg38\
			--known-splicesite-infile ~/HumanGenome/hg38.105_spliceSites.txt\
			-1 ${s}_1.fq -2 ${s}_2.fq -S ${s}_${c}.sam\
			--summary-file ${s}_${c}.summmary.txt
			samtools view -S -b ${s}_${c}.sam > ${s}_${c}.bam
			rm ${s}_${c}.sam
			rm ${s}_1.fastq
			rm ${s}_2.fastq
			rm ${s}_1.fq 
			rm ${s}_2.fq
			break
		elif [ $p = S ] 
		then
			fastq-dump -v $s.sra
			seqtk trimfq $s.fastq > $s.fq
			hisat2 -x ~/HumanGenome/hg38\
			--known-splicesite-infile ~/HumanGenome/hg38.105_spliceSites.txt\
			-U $s.fq -S ${s}_${c}.sam --summary-file ${s}_${c}.summmary.txt
			samtools view -S -b ${s}_${c}.sam > ${s}_${c}.bam
			rm ${s}_${c}.sam
			rm $s.fastq 
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

	samtools view -b ${s}_${c}.sort.bam "chr10:44370165-44385092" > ${s}_${c}/${s}_${c}.cxcl12.bam
	samtools depth -d 0 ${s}_${c}/${s}_${c}.cxcl12.bam > ${s}_${c}/${s}_${c}.cxcl12.bam.txt
		
	samtools view -b ${s}_${c}.sort.bam "chr7:5527147-5530601" > ${s}_${c}/${s}_${c}.actb.bam
	samtools depth -d 0 ${s}_${c}/${s}_${c}.actb.bam > ${s}_${c}/${s}_${c}.actb.bam.txt

	cd ${s}_${c}
	for bamfile in *.bam; do
		samtools view -H $bamfile > temp.sam
		samtools view $bamfile | awk '($6 ~ /N/)' >> temp.sam
		samtools view -bS temp.sam > splice_$bamfile
		samtools depth -d 0 splice_$bamfile > splice_$bamfile.txt
	done
	rm temp.sam
	samtools index ${s}_${c}.cxcl12.bam
	cd ..
	
	rm ${s}_${c}.sort.bam
	rm ${s}_${c}.sort.bam.bai
done