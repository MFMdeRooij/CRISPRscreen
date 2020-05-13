#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#####################################################################################################
# We used the screen design of Jastrzebski et al Methods Mol Biol 2016
# Open in Spyder, adjust the settings, and run the script
# If you don't have Biopython installed comment out line 15, (or install biopython: command line/Windows Powershell: conda install biopython)
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2019, info: m.f.derooij@amsterdamumc.nl
#####################################################################################################
# Import python modules
import os
import numpy as np
import pandas as pd
import gzip
import sys
from Bio import pairwise2
#####################################################################################################
#                                            SETTINGS
# Workdirectory
os.chdir("H:/BioWin/Screens/")

# Screen files
screens = ("test.fastq.gz",) 

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
barcode = ["CGTGAT", "ACATCG", "GCCTAA", "TGGTCA", "CACTGT", "ATTGGC", 
		   "GATCTG", "TCAAGT", "CTGATC", "AAGCTA", "GTAGCC", "TACAAG"]

# Number of allowed mismaches in barcodes, and max indels in constant part of the reads (PCR-primer mutations) (default = 1,3)
BCmm,indel = 1,3

# Print information about non-perfect barcode and guide sequences: 0 = No, 1 = yes
printUnID = 1

# Use biopython pairwise2 alignment of barcodes and guides: 0 = No, 1 = yes
biopythonBarcode = 0
biopythonGuide = 0 # Not recommanded (mutations could be inactivating the guide)

######################################################################################################
    
# Make dictionary with barcode numbers
BCnum = dict(zip(barcode, np.arange(12)+1))

# Make dictionary with barcode pieces 
BCpieces = {}
BCsize = len(barcode[0])
for index, bc in enumerate(barcode, start=1):
  for nt in range(BCsize):
    BCpieces[str(index)+"."+str(nt+1)] = bc[nt]

# Load library and make count table precursor
Library = pd.read_csv(library)
Library = pd.DataFrame([["Unidentified"]*len(Library.columns)], columns=Library.columns.tolist()).append(Library,ignore_index=True)
LibrarySeq = Library["Sequence"].to_numpy()

def barcodeDetermination(seq):
    # Determine barcode
    BCreads = seq[0:BCsize]
    if printUnID == 1:
        global BCmismatch
        global BCalign
    try:
        BCnumber = BCnum[BCreads]
    except KeyboardInterrupt:
        sys.exit()
    except:
        BCnumber = 0
        if BCmm > 0:
              BCread = {}
              for nt in range(BCsize):
                  BCread["s"+str(nt+1)] = seq[nt]  
        # Calculate similarities with barcode dictionary
        for ibc in range(len(barcode)):
            m1 = 0
            for i in range(BCsize):
                if BCread["s"+str(i+1)] == BCpieces[str(ibc+1)+"."+str(i+1)]:
                    m1+=1
            if m1 >= BCsize-BCmm:
                BCnumber = ibc+1
                if printUnID == 1:
                    print("This barcode (%s --> %s) had %d mismatch" % (BCreads, barcode[BCnumber-1], BCsize-m1)) 
                    BCmismatch+=1
                    print("Barcodes with mismatches mapped: %d" % (BCmismatch))
                break  
           
    if BCnumber == 0:
        if biopythonBarcode == 1:
            BCreads = seq[0:BCsize+1]
            bcarr = np.array([])
            for bc in barcode:
                alignments = pairwise2.align.localms(BCreads, bc, 5, -3, -5, -5, score_only=True)
                bcarr=np.append(bcarr, int(alignments)) 
            if bcarr.max()>=15:    
                BCnumber = 1+bcarr.argmax() 
                if printUnID == 1:
                    print("Barcode (%s --> %s) aligned by biopython" % (BCreads, barcode[BCnumber-1]))
                    BCalign+=1
                    print("Barcodes which are aligned by biopython: %d" % (BCalign))
            else:
                BCnumber = 0
    if printUnID == 1:
        if BCnumber == 0:
            print(seq+"This barcode cannot be identified: "+BCreads)
    return(BCnumber)
   
def guideDetermination(seq):
    guide = seq[loc:loc+CRISPRsize]
    if printUnID == 1:
        global total
        global indels
        global align
        total+=1
        print("\nTotal reads analyzed: %d" % (total))
    # When there is an indel in constant region, change CRISPR seq location 
    if indel > 0:
        upseqlength = len(upseq)
        seqcheck = seq[loc-upseqlength:loc]   
        if upseq != seqcheck: 
            for i in range(-indel,indel):
                preseq = seq[loc-upseqlength+i:loc-upseqlength+upseqlength+i]
                if preseq == upseq:
                    guide = seq[loc+i:loc+CRISPRsize+i]
                    if printUnID == 1:
                        print("This read (%s) had %d indel" % (seq, i))
			indels+=1
                        print("Reads with indels mapped: %d" % (indels))
                    break
                    
    guideNumber=np.where(LibrarySeq==guide)[0]    
    if guideNumber.size == 0: 
        if biopythonGuide == 1:
            guidereads = seq[loc-3:]
            guidearr = np.array([])
            for guides in LibrarySeq:
                alignments = pairwise2.align.localms(guidereads, guides, 5, -3, -5, -5, score_only=True)
                guidearr=np.append(guidearr, int(alignments)) 
            if guidearr.max()>=75:    
                guideNumber = guidearr.argmax() 
                if printUnID == 1:
                    print("Guide (%s --> %s) aligned by biopython" % (guidereads, LibrarySeq[guideNumber])) 
                    align+=1
                    print("Reads which are aligned by biopython: %d" % (align))  
            else:
                guideNumber = 0
        else:
            guideNumber = 0       
    if printUnID == 1:
        if guideNumber==0:     
            print(seq+"This guide cannot be identified: "+guide)
    return guideNumber

# Count the reads
for screen in screens:
    if printUnID == 1:
        total = 0
        indels = 0
        align = 0
        BCmismatch = 0
        BCalign = 0
    countTableToFillZero = np.zeros((len(Library), len(barcode)+1))
    # Read Fastq file and put counts in the count table
    screenData = screen
    with gzip.open(screenData, 'r') as FASTQ:
	    while True:
    		FASTQ.readline()
    		seq = FASTQ.readline().decode("utf-8")
    		if not seq:
    		    break
    		FASTQ.readline()
    		FASTQ.readline()
    		countTableToFillZero[guideDetermination(seq),barcodeDetermination(seq)]+=1

    #Write count table to file
    if printUnID == 1:
        print("%d reads were analyzed, in which %d (%d%%) were mapped" % (int(countTableToFillZero.sum()), 
                                int(countTableToFillZero[1:,1:].sum()), int(countTableToFillZero[1:,1:].sum())/
                                                                            int(countTableToFillZero.sum())*100))
    countTableFilled = pd.concat([Library, pd.DataFrame(countTableToFillZero, columns=["Unidentified"]+barcode)], axis=1)
    countTableFilled.to_csv("CountTable_"+screen+".csv", index=False)
