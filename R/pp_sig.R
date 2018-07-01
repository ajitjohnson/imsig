#' @title Pre-processing ImSig file
#' @description Subsets ImSig genes based on the genes that are common to the users dataset and ImSig
#' @param exp Dataframe of transcriptomic data (natural scale) containing genes as rows and samples as columns. Note: Gene names should be set as row names and duplicates are not allowed. Missing values are not allowed within the expression matrix. Check example- head(example_data): \code{\link{example_data}}.
#' @return ImSig dataframe

pp_sig <- function(exp){
  g <- Reduce(intersect, list(as.character(row.names(exp), sig$gene)))
  sig <- sig[sig$gene %in% as.character(g),]
  return(sig)
}
