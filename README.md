# imsig
Immune Cell Gene Signatures for Profiling the Microenvironment of Solid Tumours

ImSig is a set of gene signatures that can be used to estimate the relative abundance of immune cells in tissue transcriptomics data.

# Requirement
*Expression matrix* of transcriptomics data, with HGNC gene symbols as rows and samples as columns. The Gene symbols need to be set as rownames and duplicates are not allowed. The expression matrix can be imported into R similar to the basic syntax below.

`exp = read.table('exp.txt', header = T, row.names = 1, sep = '\t')`


# Usage
