# CRISPR screen analysis from fastq to robust statistics in Linux or Microsoft Windows.
With my scripts you can analyze your own CRISPR screen completely and in a robust way with little knowledge of R, Perl, and Python. 
The scripts are compatible with Linux as well as Microsoft Windows (with a small adjustment in the Perl scripts).
Open the scripts (as text file) for information how to use them.
- Install Perl, R, R studio, and Anaconda3 (during installation of Anaconda3, add Python to your path (check the box)).
- From fastq to count table use FastqToCountTable.pl with a library file and your fastq files (Perl). 
If you don't like Perl, you can also use the Python script (FastqToCountTable.py) in Spyder, but this one is much slower than the Perl  version.
- From count table to robust statistics use CRISPRScreenAnalysis.R with CRISPRScreenAnalysisLibraries.csv and your count tables (R).
- To get some information about the quality of the read mapping use in addition to FastqToCountTable.pl also Quality.pl and your fastq files, and subsequently Quality.R (Perl/R). With FastqToCountTable_Quality.py you can also get some insights in the unmappable reads (Python). Our experience is that up to 30% of the library reads contain mismatches in the guide sequences.
- To interactively view the MA plots, add the output from CRISPRScreenAnalysis.R to the files/DataCRISPR03New folder of MAviewer, and open MAviewer with the sh file (linux) or bat file (Windows), and discover your screen results (You can add a shortcut to these files for on your desktop. As an example I included the published screens of B cell lymphoma cell lines of Phelan et al, Nature 2018 (Python).
- To make high quality MA plots, guide fold change plots, and split-volcano plots, use the MAplot_GeneOfInterest.py, guideFCPlot.py, and VolcanoSplit_GeneOfInterest.py scripts in Spyder (Python).
- To look at gene expression of your CRISPR targets, you can download RNAseq data from most cell lines from the NCBI-SRA database, and analyze it with GeneExpressionRNAseqHisat2.R. Because this script uses various Linux apps, this does not work on Windows (R/Linux).
- In the TestData folder is a test file available (test.fastq.gz), which contains the first 750k reads of the adhesion CRISPR screen. Barcode 1-3 corresponds to the preadhesion replicates, Barcode 7-9 to the PMA-induced adhesion replicates, and Barcode 4-6 to the anti-IgM-induced adhesion replicates (Note that this is the same design as a synthetic lethality screen (T0, Control, Treatment). The easiest way to analyze this, is to run the FastqToCountTable.py in Spyder to produce the count table (Know that the Perl script is much faster, so that is recommended for a complete FASTQ file), and subsequently the R script CRISPRScreenAnalysis.R in R studio to perform the statistics. For the TestData analysis you only have to adjust the workdirectories. You should get the same data as in the TestData/Output folder. This takes only a few minutes (Windows or Linux).
- Good luck!, Martin de Rooij

For more information, suggestions, or bugs: m.f.derooij@amsterdamumc.nl
