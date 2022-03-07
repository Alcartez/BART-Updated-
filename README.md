# README #

### What is this repository for? ###

* BART (Bioinformatics Array R Tool) is a shiny web application which performs microarray analysis. 
* BART automates the download and analysis process across diverse microarray platforms.
* BART can work with a GSE accession ID from GEO or an expression table uploaded by the user. 
* BART performs differential expression testing and generates hierarchically clustered heat maps and PCA plots to facilitate sample comparison, as well as volcano plots and pathway enrichment analysis to summarize differential expression testing.

### How do I get set up? ###

* BART is available at bart.salk.edu (recommended)

* To install and use BART from command line:
	
	1) `Git clone https://bitbucket.org/Luisa_amaral/bart.git`
	
	2) `cd bart`

	3) Open `R`

	4) `>install.packages(“shiny”)`        

	5) `>library(“shiny”)`

	6) `>runApp(".")` 	

	7) Allow BART to attempt to automatically install packages (may request some permissions)

	8) BART will appear in an interactive browser window for use
	
### How does BART work? ###

* BART automatically parses any experimental data available on GEO and suggests groupings to use for differential expression testing
* If CEL files are used as input, BART detects whether to use the Affy (Gautier et al., 2004) or Oligo (Carvalho et al., 2010) platforms and performs RMA normalization that attempts to remove local biases across samples in order to enable meaningful differential expression testing
* If starting from an expression table, BART automatically detects whether a log2 transformation is needed and normalizes the data
* The hclust R function is used to perform clustering of the top 1000 expressed normalized genes
* Principal Component Analysis is used on the top 1000 expressed genes to transform and visualize samples as a function of the top two independent descriptive variables, which show the most sample variance
* BART leverages the LIMMA bioinformatics package (Ritchie et al., 2015) to perform differential expression testing
* BART automatically generates volcano plots for each differential expression comparison and highlights genes that are differentially expressed with adjusted p-values less than 0.05
* BART performs functional enrichment analysis by using the GAGE and WebGestAlt R package to determine which KEGG pathways are enriched
* All data and plots can then be downloaded for additional analysis and validation

### Who do I talk to? ###

* Repo owner or admin