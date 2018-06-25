#' @title Generation of a correlation matrix
#' @description Creates a correlation matrix using ImSig signature genes.
#' @param exp Expression matrix with genes as rows and samples as columns [data.frame]
#' @param r Correlation cut-off to use. Correlations below the defined cut-off will be replaced with zero before returning the final correlation matrix.
#' @return Gene-gene correlation matrix of the ImSig genes.
#' @importFrom HiClimR fastCor
#' @export

corr_matrix <- function(exp, r){
  sig_su <- sig[which(as.character(sig$gene) %in% row.names(exp)), ]
  data <- exp[as.character(sig_su$gene),]
  cor_data <- fastCor(t(data))
  cor_data[cor_data <= r]= 0
  diag(cor_data) = 0
  return(cor_data)
}
