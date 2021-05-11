# MasterThesis: Epigenetic regulation and transcriptional effects of full-length L1HS retrotransposons in human neural progenitor cells
Streamlined Snakemake pipelines and Rscripts created to define the epigenetic status of full-lenght LINE-1 and HERV elements in human neuron progenitor cells. 
- Author: Ninoslav Pandiloski
## Table of Contents:
- [Background][1]
- [Tools][2]
- [Usage][3]

## Background:
Three Snakemake pipelines were created to analyze the dataset in the project: Epigenetic.sn - defining the epigenetic status around FL-TEs, Bulk-RNA.sn - capturing the transcriptional activity of the elements, and DESeq.sn - counting reads and preparing the data for downstream R analysis and visualization. The Rscript asigns unique gene names to the counts and performs differential expression analysis of Human vs Chimpanzee datasets using DESeq2. Furthermore, configuration files were created such as wildcards.yaml - containing all the wildcards used throughout as well as lunarc_config - with all configurations to run the Snakemake files onto the Lunarc-LSENS server.

## Tools:
A variety of tools were used throughout the project, for quality assesment, file transformation, and extraction of information
- Tools used: 
	- FastQC(0.11.8) and MultiQC(1.2) - Quality assesment of raw data;
	- Bowtie2(2.3.4.2) and STAR(2.6.0c) - Mapping the epigenetics CUT&RUN data and the Bulk-RNAseq data to hg38, respectively;
	- deepTools(2.5.4) - Computing signal matrices over genomic region and visualizing them into heatmaps;
	- R(3.6.0) and SEACR(1.3) - R was used as a prerequisite module for SEACR's peak calling of the CUT&RUN data;
	- Subread(1.6.3) - Subread's featureCounts was used to count the number of aligned reads to genomic features;
	- R(4.0.0) - Used for further analysis of the featureCounts output and visualization;
	- Snakemake (5.4.2) - Used to build the pipelines used for our analysis

## Usage
The pipelines were optimized to be run onto the Lunarc-LSENS cluster using the configurations files using. They were ran using the following command:
- ```snakemake -n -s Snakefile.sn -j 6 --cluster-config /projects/fs1/nino/MastersThesis/bin/config_files/lunarc_config.json --cluster "sbatch -A {cluster.account} -p {cluster.partition} --tasks-per-node {cluster.tasks-per-node} -t {cluster.time} -o {cluster.o} -e {cluster.e} -J {cluster.J} -N {cluster.N}" --latency-wait 120 --rerun-incomplete```


---- 

[1]:	#Background
[2]:	#Tools
[3]:	#Usage
