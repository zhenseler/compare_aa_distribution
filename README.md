# compare_aa_distribution.py

## Purpose
	Compare how site-specific amino acid frequencies differ in a particular gene between
	two user-supplied groups of samples

## Dependencies
	- diamond (https://anaconda.org/bioconda/diamond)
	- mafft
	- dit
	
## Usage
	python compare_aa_distribution.py <database_path> <sample_fastq_path> <nMA_path>


## Input
	1. A protein database of the gene you want to examine (dmnd)
	2. Quality metagenomic reads (fasta)

## Workflow
	1. Creates a Master Alignment(MA) from protein database
	2. Number positions in MA (can also remove sites with low coverage)
	3. Pull out metagenomic reads mapping to database via DIAMOND
	4. From best hit, associate each position in each read with its corrosponding
		position in the numbered MA
	5. Construct an amino acid frequency table for each sample
	6. Compare amino acid frequency tables between groups

## Outputs
	- A count table showing the counts of each amino acid at each position within
		the numbered master alignment
	- A frequency table of the count table
	

-Z
