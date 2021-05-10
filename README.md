# CRISPR screen analysis from fastq to robust statistics in Linux or Microsoft Windows.
With my scripts you can analyze your own CRISPR screen completely and in a robust way with little knowledge of R, Perl, and Python. 
The scripts are compatible with Linux as well as Microsoft Windows (with a small adjustment in the Perl scripts).
Open the scripts (as text file) for information how to use them.
- Install Perl, R, R studio, and Anaconda3 (during installation of Anaconda3, add Python to your path (check the box)).

CRISPRScreenAnalysis:
- From fastq to count table use FastqToCountTable.pl with a library file and your fastq files (Perl). 
- If you don't like Perl, you can also use the Python script (FastqToCountTable.py) in Spyder, but this one is much slower than the Perl script (Python).
- From count table to robust statistics use CRISPRScreenAnalysis.R with CRISPRScreenAnalysisLibraries.csv and your count tables (R).
- To get information about the quality of the read mapping (e.g. how much reads map to the library, contaminations, primerdimers), use Quality.pl + Quality.R (Perl/R). To obtain more information (e.g. mutations in barcodes/guides) use FastqToCountTable.py (Python).

MAviewer:
- To interactively view the MA plots, add the output from CRISPRScreenAnalysis.R to the files/DataCRISPR03New folder of MAviewer, and open MAviewer with the sh file (Linux) or bat file (Windows), and discover your screen results (You can add a shortcut to these files for on your desktop. As an example I included the published screens of B cell lymphoma cell lines of Phelan et al, Nature 2018. You can also check the gene expression (RNAseq) of your top hits in a small B cell line panel (derived from the public SRA database, in which most is performed by the Broad Institute). The RNAseq data are in TPM, which are subseqently normalized by median of ratios (DESeq2) between cell lines (Python).

CRISPRScreenPlots:
- To produce publishing-grade MA plots, split-volcano plots, and guide fold change plots use the MAplot_GeneOfInterest.py, VolcanoSplit_GeneOfInterest.py, and guideFCPlot.py scripts in Spyder (Python).

RNAseq:
- To look at gene expression of your CRISPR targets, you can download RNAseq data from most cell lines from the NCBI-SRA database, and analyze it with GeneExpressionRNAseqHisat2.R. Because this script uses various Linux apps, this does not work on Windows (R/Linux). With the Bash scripts, you can map the RNAseq data with Hisat2 outside R, which is much faster. Afterwards you can continue in R to produce count tables, or view the read mapping in IGV viewer (e.g. to look at mutations, differential splicing between cell lines, or check CRISPR induced deletions).

SangerSeq:
- After performing the CRISPR screen, you might validate the hits. With the Sanger script, you can analyze your Sanger sequences from your CRISPRguide cloning in a few seconds. Apart from Anaconda3, you need to install the Biopython package (Python). 

Synergy:
- After a synthetic lethality CRISPR screen in which a nice drug is available for the best hits, you can check for synergy with the initial drug. When you have nice S-curves you can use the 4 component fit script. When that one is not working, try the Cubic fit script (R).

OncoPrint:
  - To summarize the CRISPR screen hits for multiple cell lines, you can use the OncoPrint script (R). 

TestData (for CRISPR screen analysis):
- In the TestData folder is a test file available (test.fastq.gz), which contains the first 750k reads of a loss-of-adhesion CRISPR screen. Barcodes 1-3 corresponds to the preadhesion replicates, Barcodes 7-9 to the PMA-induced adhesion replicates, and Barcodes 4-6 to the anti-IgM-induced adhesion replicates (Note that this is the same design as a synthetic lethality screen (T0, Control, Treatment). The easiest way to analyze this, is to run the FastqToCountTable.py in Spyder to produce the count table (Know that the Perl script is much faster, so that is recommended for a complete FASTQ file), and subsequently the R script CRISPRScreenAnalysis.R in R studio to perform the statistics. For the TestData analysis you only have to adjust the workdirectories. You should get the same data as in the TestData/Output folder. This all takes only a few minutes (Windows or Linux). After performing the DESeq2 analysis succesfully, you can also test the MA plot, Volcano plot and GuideFC plot scripts in Spyder.
- Good luck!, Martin de Rooij, The Netherlands

For more information, suggestions, or bugs: m.f.derooij@amsterdamumc.nl
