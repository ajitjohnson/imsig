# imsig: Immune Cell Gene Signatures for Profiling the Microenvironment of Solid Tumours

ImSig is a set of gene signatures that can be used to estimate the relative abundance of immune cells in tissue transcriptomics data especially in cancer datasets.

# Requirements
### Mandatory
**Expression matrix** of transcriptomics data, with HGNC gene symbols as rows and samples as columns. The Gene symbols need to be set as rownames and duplicates are not allowed. Missing values are also not allowed within the expression matrix. The expression matrix can be imported into R using something similar to the basic syntax below.

`exp = read.table('exp.txt', header = T, row.names = 1, sep = '\t')`

### Optional
**Survival data** (event and time to event) for the set of patients is required if the effect of immune infiltration in patient prognosis is to be determined. 

# Usage

#### Install package

`install.packages("imsig")`

#### Get an idea of the imported data

`gene_stat (exp = exp, r = 0.7)`

This function returns a table of basic stats of the data. 
