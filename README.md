# CRISPRscreen analysis from fastq to robust statistics in Linux or Microsoft Windows.
With my scripts you can analyze your own CRISPR screen completely and in a robust way with little knowledge of R, Perl, and Python. 
The scripts are compatible with linux as well as Microsoft Windows (with a small adjustment in the Perl script).
Open the scripts (as text file) for information how to use them.
- Install Perl, R, R studio, and Anaconda3 (during installation of Anaconda3, add Python to your path (check the box)).
- From fastq to count table use FastqToCountTable.pl with a library file and your fastq files (Perl).
- From count table to robust statistics use CRISPRScreenAnalysis.R with CRISPRScreenAnalysisLibraries.csv and your count tables (R).
- To get some information about the quality of the read mapping use in addition to FastqToCountTable.pl also Quality.pl and your fastq files, and subsequently Quality.R (Perl/R). 
- To interactively view the MA plots, add the output from CRISPRScreenAnalysis.R to the files/DataCRISPR03New folder of MAviewer, and open MAviewer with the sh file (linux) or bat file (Windows), and discover your screen results (You can add a shortcut to these files for on your desktop. As an example I included the published screens of B cell lymphoma cell lines of Phelan et al, Nature 2018 (Python).
- To make an high quality MA plot, use the MAplot_GeneOfInterest_v1.1.py script (Python).

For more info, suggestions, or bugs: M.F.deRooij@amsterdamumc.nl
