#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#####################################################################################################
# We used the screen design of Jastrzebski et al Methods Mol Biol 2016
# Open in Spyder, adjust the settings, and run the script
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2019, info: m.f.derooij@amsterdamumc.nl
#####################################################################################################
# Import python modules
import os
import gzip
import csv
#####################################################################################################
#                                            SETTINGS

# Workdirectory
os.chdir("H:/BioWin/Screens/")

# Screen files
screens = ("TGACCA_Namalwa_Adhesion.fastq.gz",) 

# Library file
library = "LibraryKinomeBrunelloLentiGuide.csv"

# Library type: 1 = CRISPR, 2 = shRNA
lib = 1
if lib == 1: 
    # Location CRISPR seq-1, length CRISPR, and upstream nucleotides (lentiCRISPRv2/LentiGuide: 41, 20, "CACCG")
    loc,CRISPRsize,upseq  = 41,20,"CACCG"
if lib == 2: 
    # Location shRNA seq-1 and length shRNA, and upstream nucleotides (pKLO.1 TRC1/1.5 and TRC2: 42, 21, "ACCGG")
    loc,CRISPRsize,upseq = 42,21,"ACCGG"

# Barcodes (all of the same length)
barcode = ("CGTGAT", "ACATCG", "GCCTAA", "TGGTCA", "CACTGT", "ATTGGC", 
	       "GATCTG", "TCAAGT", "CTGATC", "AAGCTA", "GTAGCC", "TACAAG")

# Number of allowed mismaches in barcodes, and max indels in constant part of the reads (PCR-primer mutations) (default = 1,3)
BCmm,indel = 1,3

######################################################################################################
# Make dictionary with barcode numbers
BCnum = {}
for i in range(len(barcode)):
    BCnum[barcode[i]]=i+1

# Make dictionary with barcode pieces 
BCpieces = {}
BCsize = len(barcode[0])
index = 0
for bc in barcode:
  index+=1
  for nt in range(BCsize):
    BCpieces[str(index)+"."+str(nt+1)] = bc[nt]

# Load library and make count table precursor
libraryGuides,libraryGenes,librarySeqs = [],[],[]
guides = open(library, 'r')
line = guides.readline()
while line:
    listLine = line.split(",")
    libraryGuides.append(listLine[0])
    libraryGenes.append(listLine[1])
    librarySeqs.append(listLine[2].strip("\n"))
    line = guides.readline()
guides.close()    
# List of lists (multidimensional list)
countTablePre = [libraryGuides, libraryGenes, librarySeqs]

for screen in screens:
    # Make multidimensional list with read sequences 
    readTable = []
    for i in range(len(barcode)):
        readTable.append([])
    
    # Read Fastq file and write reads in the right barcode file
    screenData = screen
    FASTQ = gzip.open(screenData, 'r')
    while True:
        name = FASTQ.readline()
        if not name:
            break
        seq = FASTQ.readline()
        plusline = FASTQ.readline()
        qual = FASTQ.readline()
        seq = seq.decode("utf-8")
        
        # Determine barcode
        BCnumber = 0
        BCreads = seq[0:BCsize]
        if BCreads in barcode:
            # When barcode has no mismatches
            BCnumber = BCnum[BCreads]
        else:
            if BCmm > 0:
                  BCread = {}
                  for nt in range(BCsize):
                      BCread["s"+str(nt+1)] = seq[nt]  
            # Calculate similarities with barcode hash
            for ibc in range(len(barcode)):
                m1 = 0
                for i in range(BCsize):
                    if BCread["s"+str(i+1)] == BCpieces[str(ibc+1)+"."+str(i+1)]:
                        m1+=1
                if m1 >= BCsize-BCmm:
                    BCnumber = ibc+1
                    break
    
        # Determine CRISPR sequence
        if BCnumber > 0:
           seqtag = seq[loc:loc+CRISPRsize]
           # When there is an indel in constant region, change CRISPR seq location 
           if indel > 0:
               upseqlength = len(upseq)
               seqcheck = seq[loc-upseqlength:loc]   
               if upseq != seqcheck: 
                   for i in range(-indel,indel):
                       preseq = seq[loc-upseqlength+i:loc-upseqlength+upseqlength+i]
                       if preseq == upseq:
                           seqtag = seq[loc+i:loc+CRISPRsize+i]
                           break
				
           # Make multidimensional array with read sequences
           readTable[BCnumber-1].append(seqtag)
    FASTQ.close() 
    
    # Start new count table
    readCounts = []
    for i in range(len(barcode)):
        readCounts.append([])
        # Produce count tables per barcode and add to multidimensional list
        for j in librarySeqs:
            rc = readTable[i].count(j)
            readCounts[i].append(rc)
            # Name BCs
        readCounts[i][0] = "BC"+str(i+1)
    countTable = countTablePre + readCounts
    transpose = [*zip(*countTable)]
    
    # Write count table to csv file
    cell = screen.split("/") [-1]
    myCsv = open("CountTable_"+cell+".csv", "w", newline='')
    csvWriter = csv.writer(myCsv, delimiter=',')
    csvWriter.writerows(transpose)
    myCsv.close() 
