# -*- coding: utf-8 -*-
"""
Quick Sanger Sequence Analysis of cloning a sgRNA or shRNA insert into a plasmid:
    - Install the Biopython package (command line: conda install -c anaconda biopython)
    - Add the reference sequences (like the Gibson oligo sequences) to Sanger_RefsCloningOligo.csv
    - Add the ab1 files and the Sanger_RefsCloningOligo.csv file to the given folder, and run the script in Spyder.
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2021, info: m.f.derooij@amsterdamumc.nl
"""
##################################################################################################################################
#                                                           SETTINGS

folder = 'H:/BioWin/SangerSeqs'
##################################################################################################################################
import pandas as pd
import glob
import time
from Bio import SeqIO
from Bio import pairwise2

start = time.time()
files = glob.glob(folder+'/*.ab1')
files = sorted(files, key=lambda x: int("".join([i for i in x if i.isdigit()])))
open(folder+'/sangerfa.fa', 'w').close()
for num, seq in enumerate(files):
    with open(seq, "rb") as input_handle:
        with open(folder+'/sangerfa.fa', 'a') as output_handle:
            sequence = SeqIO.read(input_handle, "abi")
            SeqIO.write(sequence, output_handle, "fasta")
    print("Converted file %d" % (num+1)) 
print("\n")  
input_handle.close()
output_handle.close()

ref = pd.read_csv(folder+"/Sanger_RefsCloningOligo.csv")
ref['Sequence'] = ref['Sequence'].str.upper() 
open(folder+'/SangerSeqIDs.txt', 'w').close()
with open(folder+'/sangerfa.fa', "r") as in_handle:
    with open(folder+'/SangerSeqIDs.txt', "a") as out_handle:
        fasta_sequences = SeqIO.parse(in_handle,'fasta')
        for fasta in fasta_sequences:
             name, sequence = fasta.id, str(fasta.seq)   
             ID=0
             for refseq in ref['Sequence']:
                 alignments = pairwise2.align.localms(refseq, sequence, 5, -3, -5, -5, score_only=False)
                 if alignments[0].score == len(refseq)*5:               
                     print("%s is %s" % (name, ref["ID"][ref["Sequence"]==refseq].values[0]))
                     out_handle.write("%s is %s\n" % (name, ref["ID"][ref["Sequence"]==refseq].values[0]))
                     ID=+1
                 elif alignments[0].score > len(refseq)*4.5:
                     print("%s could be %s -> check peeks in Snapgene Viewer" % (name, ref["ID"][ref["Sequence"]==refseq].values[0]))
                     out_handle.write("%s could be %s -> check peeks in Snapgene Viewer\n" % (name, ref["ID"][ref["Sequence"]==refseq].values[0]))
                     ID=+1
                 else:
                    pass
             if ID == 0:
                print("%s cannot be identified" % name)
                out_handle.write("%s cannot be identified\n" % name)
in_handle.close()
out_handle.close()
stop = time.time()
print("\nThe Sanger sequences were analyzed in %d minutes and %d seconds\n" % ((stop-start)/60, (stop-start)%60))
