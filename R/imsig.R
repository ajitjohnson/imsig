#' @title Estimate the relative abundance of tissue-infiltrating immune subpopulations abundances using gene expression data
#' @description Estimates the relative abundance of immune cells across patients/samples.
#' @param exp Dataframe of transcriptomic data (natural scale) containing genes as rows and samples as columns. Note: Gene names should be set as row names and duplicates are not allowed. Missing values are not allowed within the expression matrix. Check example- head(example_data): \code{\link{example_data}}.
#' @param r Use a value between 0 and 1. Default is 0.6. This is a user defined correlation cut-off to perform feature selection (\code{\link{feature_select}}). Feature selection aids to enrich the prediction of relative abundance of immune cells by filtering off poorly correlated ImSig genes. To get an idea of what cut-off to use check the results of (\code{\link{gene_stat}}) and choose a cut-off that displays high median correlation and maintains a high proportion of genes after feature selection.
#' @param sort Sort the samples based on abundance of a particular cell type. `Set sort = FALSE` if you wish not to apply sorting. By default the function sorts by abundance of T cells. The cell type of interest for sorting can be controlled by the `sort_by` parameter.
#' @param sort_by Can be used to sort the samples by predicted abundance of a particular cell type. All other cell types follow this sorting. By default it is sorted by `T cells`
#' @return Relative abundance of immune cells across samples. Returns a dataframe.
#' @examples
#' cell_abundance = imsig (exp = example_data, r = 0.7, sort=TRUE, sort_by='T cells')
#' head(cell_abundance)
#' @seealso \code{\link{feature_select}}, \code{\link{example_data}}
#' @export

imsig <- function(exp, r = 0.6, sort=TRUE, sort_by='T cells'){
  exp <- pp_exp(exp)
  sig <- pp_sig(exp)
  fg <- feature_select(exp, r)
  sig_fg <- sig[ which(as.character(sig$gene) %in% fg), ]
  avg_data <- exp[row.names(exp) %in% fg,]
  cc <- data.frame(matrix(nrow = ncol(exp)))
  cc <- cc[,-1]
  for (i in levels(sig_fg$cell)){
    s <- sig_fg[sig_fg$cell %in% i,]
    e <- avg_data[as.character(s$gene),]
    e_avg <- data.frame(colMeans(e, na.rm = TRUE))
    colnames(e_avg) <- i
    cc <- cbind(cc, e_avg)
  }
  if (sort == TRUE){cc <- cc[sort.list(unlist(data.frame(cc[,sort_by]))),]}
  return(cc)
}

