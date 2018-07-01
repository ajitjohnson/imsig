#' @title Correlation matrix
#' @description Creates a correlation matrix of ImSig signature genes.
#' @param exp Dataframe of transcriptomic data (natural scale) containing genes as rows and samples as columns. Note: Gene names should be set as row names and duplicates are not allowed. Missing values are not allowed within the expression matrix. Check example- head(example_data): \code{\link{example_data}}.
#' @param r Use a value between 0 and 1. Default is 0.6. This is a user defined correlation cut-off to perform feature selection (\code{\link{feature_select}}). Feature selection aids to enrich the prediction of relative abundance of immune cells by filtering off poorly correlated ImSig genes. To get an idea of what cut-off to use check the results of (\code{\link{gene_stat}} and choose a cut-off that displays high median correlation and maintains a high proportion of genes after feature selection.
#' @return Gene-gene correlation matrix of ImSig genes.
#' @importFrom HiClimR fastCor

corr_matrix <- function(exp, r){
  exp <- pp_exp(exp)
  sig <- pp_sig(exp)
  cor_data <- fastCor(t(exp))
  cor_data[cor_data <= r]= 0
  diag(cor_data) = 0
  return(cor_data)
}
