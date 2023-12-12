# CRISPR screen analysis from fastq to robust statistics in Linux or Microsoft Windows.
With my scripts you can analyze your own CRISPR screen completely and in a robust way with little knowledge of R, Perl, and Python. 
The scripts are compatible with Linux as well as Microsoft Windows with a small adjustment in the Perl scripts, and with the exception of the RNAseq scripts and the MAGeCK option in the CRISPRScreenAnaysis.
Open the scripts (as text file) for information how to use them.
- Install Perl, R, R studio, and Anaconda3 (during installation of Anaconda3, add Python to your path (check the box)).

CRISPRScreenAnalysis:
- From fastq to count table use FastqToCountTable.pl with a library file and your fastq files (Perl). 
- If you don't like Perl, you can also use the Python script (FastqToCountTable.py) in Spyder, but this one is much slower than the Perl script (Python).
- From count table to robust statistics use CRISPRScreenAnalysis.R with CRISPRScreenAnalysisLibraries.csv and your count tables (R).
- GSEA_CRISPRscreen.R can be used to obtain pathway information about the CRISPR screen results.
- To get information about the quality of the read mapping (e.g. how much reads map to the library, contaminations, primerdimers), use Quality.pl + Quality.R (Perl/R). To obtain more information (e.g. mutations in barcodes/guides) use FastqToCountTable.py (Python).

CRISPRScreenPlots:
- To produce publishing-grade MA plots, split-volcano plots, and guide fold change plots use the MAplot_GeneOfInterest.py, VolcanoSplit_GeneOfInterest.py, and guideFCPlot.py scripts in Spyder (Python). MAplot_GenesOfInterest.R can be used to make MA plots in Rstudio, and is much faster than the Python script. With the SyntheticLethality_GenesOfInterest.R or .py scripts, you can get a nice visualization of your synthetic lethal hits (R or Python). 

CRISPRScreenPublic:
- To check the quality of lethality screens, of which only CRISPR scores are published, by comparing the distibution of essential and non-essential genes, use CRISPRscreenAnalysisFromScores.R script (R). 
- To plot CRISPRscreen or RNAseq values from DepMap in a PDF file, download the data files from DepMap.org and use CRISPRscores_DepMap.R or RNAseq_DepMap.R (R).

MAviewer:
- To interactively view the MA plots, add the output from CRISPRScreenAnalysis.R to the files/DataCRISPR03New folder of MAviewer, and open MAviewer with the sh file (Linux) or bat file (Windows), and discover your screen results (You can add a shortcut to these files for on your desktop. As an example I included the published screens of B cell lymphoma cell lines of Phelan et al, Nature 2018. You can also check the gene expression (RNAseq) of your top hits in a small B cell line panel (derived from the public SRA database, in which most is performed by the Broad Institute). The RNAseq data are in TPM, which are subseqently normalized by median of ratios (DESeq2) between cell lines (Python).

OncoPrint:
  - To summarize the CRISPR screen hits for multiple cell lines, you can use the OncoPrint.R script (R). 

RNAseq:
- To look at gene expression of your CRISPR targets, you can download RNAseq data from most cell lines from the NCBI-SRA database, and analyze it with GeneExpressionRNAseqHisat2.R and the Bash files. Because this script uses various Linux apps, it does not work on Windows (R/Linux). With the sorted bam files you can also view the read mapping in IGV viewer (e.g. to look at mutations, differential splicing between cell lines, or check CRISPR induced deletions). From the generated count tables you can perform differential expression, gene set enrichment analysis, make gene level bar/PCA plots, or gene family pie plots (R). You can also look at (predicted) cell-cell interactions between different cell types or autocrine factors. Furthermore, you can look at differential splicing of a particular gene.

SangerSeq:
- After performing the CRISPR screen, you might validate the hits. With the Sanger script, you can analyze your Sanger sequences from your CRISPRguide cloning in a few seconds. Apart from Anaconda3, you need to install the Biopython package (Python). 

Synergy:
- After a synthetic lethality CRISPR screen in which a nice drug is available for the best hits, you can check for synergy with the initial drug. When you have nice S-curves you can use the 4 component fit script. When that one is not working, try the Cubic fit script (R). For these analyses you use the single drugs and the diagonal combinations. You can use BlissSynergy.R in Rstudio (R) or BlissSynergy.py in Spyder (Python) to calculate the Bliss synergy of a matrix-design.  

TestData (for CRISPR screen analysis):
- In the TestData folder is a test file available (test.fastq.gz), which contains the first 750k reads of a loss-of-adhesion CRISPR screen. Barcodes 1-3 corresponds to the preadhesion replicates, Barcodes 7-9 to the PMA-induced adhesion replicates, and Barcodes 4-6 to the anti-IgM-induced adhesion replicates (Note that this is the same design as a synthetic lethality screen (T0, Control, Treatment). The easiest way to analyze this, is to run the FastqToCountTable.py in Spyder to produce the count table (Know that the Perl script is much faster, so that is recommended for a complete FASTQ file), and subsequently the R script CRISPRScreenAnalysis.R in R studio to perform the statistics. For the TestData analysis you only have to adjust the workdirectories. You should get the same data as in the TestData/Output folder. This all takes only a few minutes (Windows or Linux). After performing the DESeq2 analysis succesfully, you can also test the MA plot, Volcano plot and GuideFC plot scripts in Spyder. Optionally you can download the complete data set (https://www.ncbi.nlm.nih.gov/sra/?term=SRR16971271), and compare your output with the data in our publication (de Rooij et al Nature Communications 2022; the count table and other output tables are available in the Source Data file). There can be little differences due to updates of the used packages.
- Good luck!, Martin de Rooij, The Netherlands

Publications:
- Martin F.M. de Rooij, Yvonne J. Thus, Nathalie Swier, Roderick L. Beijersbergen, Steven T. Pals, and Marcel Spaargaren. A loss-of-adhesion CRISPR-Cas9 screening platform to identify cell adhesion-regulatory proteins and signaling pathways. Nature Communications 13(1):2136 (2022). https://doi.org/10.1038/s41467-022-29835-y
- Yvonne J. Thus, Martin F.M. de Rooij, Roderick L. Beijersbergen, and Marcel Spaargaren. An Unbiased CRISPR-Cas9 Screening Method for the Identification of Positive and Negative Regulatory Proteins of Cell Adhesion. Bio-protocol Journal 12(21):e4545 (2022). https://doi.org/10.21769/BioProtoc.4545 
- Yvonne J. Thus, Martin F.M. de Rooij, Nathalie Swier, Roderick L. Beijersbergen, Jeroen E.J. Guikema, Marie Jose Kersten, Eric Eldering, Steven T. Pals, Arnon P. Kater, and Marcel Spaargaren. Inhibition of casein kinase 2 sensitizes mantle cell lymphoma to venetoclax through MCL-1 downregulation. Haematologica 108(3):797-810 (2023). https://doi.org/10.3324/haematol.2022.281668  
- Marthe Minderman, Avital Amir, Willem Kraan, Esther J.M. Schilder-Tol, Monique E.C.M. Oud, Cornelis G. Scheepstra, Arnold L. Noorduyn, Philip M. Kluin, Marie Jose Kersten, Marcel Spaargaren, and Steven T. Pals. Immune evasion in primary testicular and central nervous system lymphomas: HLA loss rather than 9p24.1/PD-L1/PD-L2 alterations. Blood. 138(13):1194-1197 (2021). https://doi.org/10.1182/blood.2021011366  

For more information, suggestions, or bugs: m.f.derooij@amsterdamumc.nl

