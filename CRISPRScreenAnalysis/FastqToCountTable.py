#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#####################################################################################################
# We used the screen design of Jastrzebski et al Methods Mol Biol 2016
# Nextera PCR primers can be found in 'PCR design Nextera-NovaSeq.xlsx'
# Apart from anaconda3, install Biopython: command line/Windows Powershell: conda install biopython
# Open in Spyder, adjust the settings, and run the script
# To produce count tables keep the 4 variables after line 47 on 0
# To have a nice overview of the quality of your reads, keep the 4 variables after line 47 on 1, run the script in Spyder
# for a few minutes, and subsequently run the read quality overview lines (select line 307-327 and press F9)
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2020, info: m.f.derooij@amsterdamumc.nl
#####################################################################################################
# Import python modules
import os
import numpy as np
import pandas as pd
import gzip
import sys
import re
from Bio import pairwise2
#####################################################################################################
#                                            SETTINGS
# Workdirectory
os.chdir("C:/BioWin/CRISPRscreen/")

# FASTQ files
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
BCmm,indel = 1,3 # for barcodes of 6 nucleotides do not use more than 1 mismatch (to see multiple matches comment out 'break' in line 105, and use printUnID = 1 in line 50)

# FOR FAST MAPPING KEEP NEXT 4 VARIABLES ON 0

# Count and print information about the identification of barcodes and guides: 0 = No, 1 = yes
printUnID = 1
# From line 115 and 200 you can comment out any block if you want to skip it, especially biopython's alignment is very slow

# Use regex and Biopython's pairwise2 alignments for barcodes and guides identification: 0 = No, 1 = yes
biopythonBarcode = 1
biopythonGuide = 1 # (Keep in mind that mutations in the guide sequence could be non-functional or aspecific)

# Remove unidentified counts from count table (recommended for further analysis) 0 = yes, 1 = no !!
RemoveUnID = 1

# Reverse-complement of the template part of the reverse primer 
revCompRevPrimerSeq = "GGCTCCGGTGCCCGTCAGTG"
# Primerdimer (dependent on primer design, high fidelity polymerases have proofreading (exonuclease activity), in which not nessesary 
# the 3'ends have to be complementary, so remove the first few nucleotides if no primerdimers can be found)

# Downstream nucleotides after guide sequence (CRISPR: "GTTTT", shRNA (TRC): "CTCGAG")
downseq = "GTTTT"
######################################################################################################   
# Make dictionary with barcode numbers
BCnum = dict(zip(barcode, np.arange(12)+1))
# Make dictionary with barcode pieces 
BCpieces = {}
BCsize = len(barcode[0])
for index, bc in enumerate(barcode, start=1):
  for nt in range(BCsize):
    BCpieces[str(index)+"."+str(nt+1)] = bc[nt]

# Load library
Library = pd.read_csv(library)
Library = pd.concat((pd.DataFrame([["Unidentified"]*len(Library.columns)], columns=Library.columns.tolist()),Library), ignore_index=True)
LibrarySeq = Library["Sequence"].to_numpy()

def barcodeDetermination(seq):
    '''Identify barcodes'''
    BCreads = seq[0:BCsize]
    if printUnID == 1:
        global BCcorrect, BCmismatch, BCdel, BCins, BCalign
    try:
        BCnumber = BCnum[BCreads]
        if printUnID == 1 and BCnumber.size != 0:         
            BCcorrect+=1
            print("Correct barcodes: %d" % (BCcorrect))  
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
                        print("This barcode (%s --> %s) has %d mismatch(es)" % (BCreads, barcode[BCnumber-1], BCsize-m1)) 
                        BCmismatch+=1
                        print("Barcodes identified with max %d mismatch(es): %d" % (BCmm,BCmismatch))
                    break  
                         
            if biopythonBarcode == 1:
            # You can comment out any block if you want to skip it
            
                if BCnumber == 0:
                    # Deletion
                    BCreads = seq[0:BCsize-1]
                    matchSeq = []
                    for i in range(len(BCreads)+1):
                        pattern = re.compile(r''+str(BCreads[0:i])+'.'+str(BCreads[i:len(BCreads)]))
                        for bc in barcode:   
                            matches = pattern.findall(bc)
                            for match in matches:
                                matchSeq.append(bc)
                    if len(matchSeq) != 0:
                        BCnumber = BCnum[max(set(matchSeq), key=matchSeq.count)]
                        if printUnID == 1:
                            print("This barcode (%s --> %s) is found by regex" % (BCreads, barcode[BCnumber-1]))          
                            BCdel+=1
                            print("Barcodes identified by regex (a deletion): %d" % (BCdel))  

                if BCnumber == 0:
                    # Insertion
                    BCreads = seq[0:BCsize+1]
                    matchSeq = []
                    for i in range(len(BCreads)):
                        pattern = re.compile(r''+str(BCreads[0:i])+str(BCreads[i+1:len(BCreads)]))
                        for bc in barcode:   
                            matches = pattern.findall(bc)
                            for match in matches:
                                matchSeq.append(bc)
                    if len(matchSeq) != 0:
                        BCnumber = BCnum[max(set(matchSeq), key=matchSeq.count)]
                        if printUnID == 1:
                            print("This barcode (%s --> %s) is found by regex" % (BCreads, barcode[BCnumber-1]))          
                            BCins+=1
                            print("Barcodes identified by regex (an insertion): %d" % (BCins))  
        
                if BCnumber == 0:
                    # Alignment
                    BCreads = seq[0:BCsize+1]
                    bcarr = np.array([])
                    for bc in barcode:
                        alignments = pairwise2.align.localms(BCreads, bc, 5, -3, -5, -5, score_only=True)
                        bcarr=np.append(bcarr, int(alignments)) 
                    if bcarr.max() > BCsize*2/3*5:    
                        BCnumber = 1+bcarr.argmax() 
                        if printUnID == 1:
                            print("This barcode (%s --> %s) is aligned by biopython" % (BCreads, barcode[BCnumber-1]))
                            BCalign+=1
                            print("Barcodes identified by biopython aligment (more mutations): %d" % (BCalign))

    if printUnID == 1 and BCnumber == 0:
        print(seq+"This barcode cannot be identified: "+seq[0:BCsize])
    return BCnumber
   
def guideDetermination(seq):
    '''Identify guides'''
    guide = seq[loc:loc+CRISPRsize]   
    if biopythonGuide == 1:
        guideForInsert = seq[loc:loc+CRISPRsize+1]
    if printUnID == 1:
        global total, correct, indels, mismatch, deletion, insertion, align, dimer, otherGuide 
        total+=1
        print("\nTotal reads analyzed: %d" % (total))
    # When there is an indel in constant region, change CRISPR seq location
    if indel > 0:
        upseqlength = len(upseq)
        seqcheck = seq[loc-upseqlength:loc]   
        if upseq != seqcheck: 
            for i in range(-indel,indel+1):
                preseq = seq[loc-upseqlength+i:loc-upseqlength+upseqlength+i]
                if preseq == upseq:
                    guide = seq[loc+i:loc+CRISPRsize+i]
                    if biopythonGuide == 1:
                        guideForInsert = seq[loc+i:loc+CRISPRsize+i+1]
                    if printUnID == 1:
                        print("This read (%s) has %d indel" % (seq, i))
                        indels+=1
                        print("Reads with indels of max %dnt in upstream region of the guide: %d" % (indel,indels))
                    break                 
    guideNumber=np.where(LibrarySeq==guide)[0]
    if printUnID == 1 and guideNumber.size != 0:
        correct+=1
        print("Correct guides: %d" % (correct))  
        
    if biopythonGuide == 1:
    # You can comment out any block if you want to skip it
    
        if guideNumber.size == 0: 
            # Deletion
            Guidereads = guide[0:CRISPRsize-1]
            matchSeq = []
            for i in range(len(Guidereads)+1):
                pattern = re.compile(r''+str(Guidereads[0:i])+'.'+str(Guidereads[i:len(Guidereads)]))
                for seqG in LibrarySeq:   
                    matches = pattern.findall(seqG)
                    for match in matches:
                        matchSeq.append(seqG)
            if len(matchSeq) != 0:
                guideNumber=np.where(LibrarySeq==max(set(matchSeq), key=matchSeq.count))[0]  
                if printUnID == 1:
                    print("This guide (%s --> %s) is found by regex" % (Guidereads, LibrarySeq[guideNumber]))          
                    deletion+=1
                    print("Guides identified by regex (a deletion): %d" % (deletion))  

        if guideNumber.size == 0:         
            # Mismatch
            Guidereads = guide
            matchSeq = []
            for i in range(len(Guidereads)):
                pattern = re.compile(r''+str(Guidereads[0:i])+'.'+str(Guidereads[i+1:len(Guidereads)]))
                for seqG in LibrarySeq:   
                    matches = pattern.findall(seqG)
                    for match in matches:
                        matchSeq.append(seqG)
            if len(matchSeq) != 0:
                guideNumber=np.where(LibrarySeq==max(set(matchSeq), key=matchSeq.count))[0]   
                if printUnID == 1:
                    print("This guide (%s --> %s) is found by regex" % (Guidereads, LibrarySeq[guideNumber]))          
                    mismatch+=1
                    print("Guides identified by regex (a mismatch): %d" % (mismatch))  
                         
        if guideNumber.size == 0: 
            # Insertion
            Guidereads = guideForInsert
            matchSeq = []
            for i in range(len(Guidereads)):
                pattern = re.compile(r''+str(Guidereads[0:i])+str(Guidereads[i+1:len(Guidereads)]))
                for seqG in LibrarySeq:   
                    matches = pattern.findall(seqG)
                    for match in matches:
                        matchSeq.append(seqG)
            if len(matchSeq) != 0:
                guideNumber=np.where(LibrarySeq==max(set(matchSeq), key=matchSeq.count))[0]   
                if printUnID == 1:
                    print("This guide (%s --> %s) is found by regex" % (Guidereads, LibrarySeq[guideNumber]))          
                    insertion+=1
                    print("Guides identified by regex (an insertion): %d" % (deletion))  

        if guideNumber.size == 0: 
            # Alignment
            Guidereads = seq[loc-3:]
            guidearr = np.array([])
            for guides in LibrarySeq:
                alignments = pairwise2.align.localms(Guidereads, guides, 5, -3, -5, -5, score_only=True)
                guidearr=np.append(guidearr, int(alignments)) 
            if guidearr.max() > CRISPRsize*2/3*5:    
                guideNumber = guidearr.argmax() 
                if printUnID == 1:
                    print("This guide (%s --> %s) is aligned by biopython" % (Guidereads, LibrarySeq[guideNumber])) 
                    align+=1
                    print("Guides identified by biopython alignment (more indels and/or mismatches): %d" % (align))  

        if guideNumber.size == 0 and printUnID == 1: 
            # primerdimer or other CRISPR guides
            Guidereads = seq
            if re.search(str(revCompRevPrimerSeq), Guidereads):
                print("This read is a primerdimer: %s" % (Guidereads))
                dimer+=1
                print("Primerdimers: %d" % (dimer)) 
            else:
                pattern = re.compile(r''+upseq+'(\w{'+str(CRISPRsize)+'})'+str(downseq))
                match = pattern.search(Guidereads)
                if match != None:
                    print("This guide (%s --> %s) is found by regex, but cannot be found in the library" % (Guidereads, match.group(1)))          
                    otherGuide+=1
                    print("Unidentified guides (contamination?): %d" % (otherGuide)) 

    if guideNumber.size == 0:
        guideNumber = 0       
        if printUnID == 1:    
            print(seq+"This guide cannot be identified: "+guide)   
    return guideNumber

# Count the reads
for screen in screens:
    if printUnID == 1:
        total=correct=indels=mismatch=deletion=insertion=align=dimer=otherGuide = 0
        BCcorrect=BCmismatch=BCdel=BCins=BCalign = 0
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

    # Overview read quality
    print("\nSUMMARY (%s):\n%d reads were analyzed, in which %d (%d%%) were mapped" % (screen, int(countTableToFillZero.sum()), 
        int(countTableToFillZero[1:,1:].sum()), int(countTableToFillZero[1:,1:].sum())/int(countTableToFillZero.sum())*100))
    if printUnID == 1:
        print("Correct barcodes: %d" % (BCcorrect))  
        if BCmm > 0:
            print("Barcodes identified with max %d mismatch(es): %d" % (BCmm,BCmismatch))
        if biopythonBarcode == 1:
            print("Barcodes identified by regex (a deletion): %d" % (BCdel))
            print("Barcodes identified by regex (an insertion): %d" % (BCins))
            print("Barcodes identified by biopython aligment (more mutations): %d" % (BCalign))
        if indel > 0:
            print("Reads with indels of max %dnt in upstream region of the guide: %d" % (indel,indels))
        if biopythonGuide == 1:
            print("Primerdimers: %d" % (dimer))  
            print("Correct guides: %d" % (correct)) 
            print("Guides identified by regex (a mismatch): %d" % (mismatch))
            print("Guides identified by regex (a deletion): %d" % (deletion))  
            print("Guides identified by regex (an insertion): %d" % (insertion))  
            print("Guides identified by biopython alignment (more indels and/or mismatches): %d" % (align))        
            print("Unidentified guides (contamination?): %d" % (otherGuide)) 
        print("Check also the 'countTableFilled' table in the next step, to see the distribution of unidentified counts per barcode/guide")
            
    # Write count table to file
    countTableFilled = pd.concat([Library, pd.DataFrame(countTableToFillZero, columns=["Unidentified"]+barcode)], axis=1)
    if RemoveUnID == 0:
       countTableFilled = countTableFilled.drop([0]).drop(["Unidentified"], axis=1)
    countTableFilled.to_csv("CountTable_"+screen+".csv", index=False)
