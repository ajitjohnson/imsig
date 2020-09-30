# *imsig*: Immune Cell Gene Signatures for Profiling the Microenvironment of Solid Tumours

*ImSig* is a set of gene signatures that can be used to estimate the relative abundance of immune cells in tissue transcriptomics data especially in cancer datasets.

The basic **principle** behind *ImSig* analysis is that for a given immune cell type to be called as present in a dataset, it is not sufficient for the signature genes to be expressed but also need to be co-expressed. *ImSig* genes were designed to be co-expressed in tissue transcriptomic data and so failure to co-express in the user's dataset may indicate its absence. Once the user has identified the cell types present in their dataset, they can compute its relative abundance and carry out survival analysis with this package.

### Citation
Immune Cell Gene Signatures for Profiling the Microenvironment of Solid Tumors <br>
Ajit J. Nirmal, Tim Regan, Barbara B. Shih, David A. Hume, Andrew H. Sims and Tom C. Freeman <br>
Cancer Immunol Res November 1 2018 (6) (11) 1388-1400; DOI: 10.1158/2326-6066.CIR-18-0342 <br>

# Requirements
### Mandatory
**Expression matrix** of transcriptomics data, with HGNC gene symbols in rows and samples in columns. The Gene symbols need to be set as *rownames* and duplicates are not allowed. Missing values are also not allowed within the expression matrix. After installing and loading *ImSig*, type `head(example_data)` in the R console to view an example expression file.

The expression matrix can be imported into R using something similar to the basic syntax below. 

**`exp = read.table('exp.txt', header = T, row.names = 1, sep = '\t')`**

### Optional
**Survival data** (event and time to event) for the set of patients is required if the effect of immune infiltration in patients prognosis is to be determined. The sample names in the expression matrix and the survival data need to be the same. The package will automatically consider only the patients that are common to both the expression matrix and survival data. After installing and loading *ImSig*, type `head(example_cli)` in the R console to view an example survival file.

The survival metadata can be imported into R using something similar to the basic syntax below.

**`cli = read.table('cli.txt', header = T, row.names = 1, sep = '\t')`**

# Usage

#### Install package

**`install.packages("imsig")`**    
**`library("imsig")`**

#### To use the latest developmental version

```
if( !require(devtools) ) install.packages("devtools")
devtools::install_github( "ajitjohnson/imsig", INSTALL_opts = "--no-multiarch")

# Load the package
library("imsig")

```

#### Get an idea of the imported data

**`gene_stat (exp = exp, r = 0.7)`**

This function returns a table of basic stats on the data. I would recommend to first take a look at the number of genes that overlap between *ImSig* and the imported expression data. Like all deconvolution methods, *ImSig* works best when all genes are available. An overlap of at least 75% is recommended. Secondly, I would look at the number of feature selected genes. Feature selection essentially removes genes below the user-defined correlation threshold (*r*). A simple interpretation of this field would be that, when a large number of signature genes are lost due to feature selection, it may indicate the absence of the cell type. This is because, in contrary to other signatures and deconvolution methods *ImSig* was designed to be co-expressed in tissues and so when they do not co-express, it may indicate the absence of that cell type. Finally, I would look at the median correlation values. Again poor median correlation values may indicate the absence of the cell type. For example, if you run a blood dataset, you will find that the macrophage signature genes will show a very low median correlation value since macrophages are not present in the blood.

This function can be run with different correlation thresholds (*r*) to determine the optimal value to be used for subsequent analysis. I would suggest picking a value that retains a large number of genes after feature selection, with a reasonably high median correlation value across cell types.

#### Run *ImSig* algorithm

**`imsig (exp = exp, r = 0.7)`**

This function returns a table of the relative abundance of immune cells across samples and ordered based on the relative abundance of T cells. If you would like to plot the results you could export the results or use the inbuilt function,

**`plot_abundance (exp = exp, r = 0.7)`**

Note, that the patients are ordered in the same way as in the table generated above. 

It is also possible to generate a network graph of the *ImSig* genes. This can complement the `gene_stat` function, to visually determine which cell types are likely to be present in the user's dataset. The network graph can be generated by running the following command.

**`plot_network (exp = exp, r = 0.7)`**

The network can be fine tuned by adding additional arguments. check `?plot_network` for all available arguments.

#### Survival analysis

If survival information is available, the Cox proportional hazard ratio can be calculated. The function first calculates the relative abundance of immune cells and then divides the patients into two groups ('high' and 'low' immune content) based on their median abundance value. A survival analysis is then carried out between these two groups.

**`imsig_survival (exp = exp, cli = cli, r = 0.7, time = 'time', status= 'status')`**

Here, `cli` is the survival dataframe, `time` is the columnname of the time to event data and `status` is the columnname of the event data within the survival dataframe (`cli`).

Similarly, the hazard ratio's can also be plotted using the following function.

**`plot_survival (exp = exp, cli = cli, r = 0.7, time = 'time', status= 'status')`**

# Issues
If there are any issues please report it at https://github.com/ajitjohnson/imsig/issues
