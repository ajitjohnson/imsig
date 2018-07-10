#' @title Feature selection of signature genes
#' @description ImSig genes were designed to be co-expressed in tissue transcriptomic data. However, depending on the dataset some of the genes may not co-express with the dominant module. In order to remove such deviant genes, a feature selection can be carried out based on correlation. This function removes genes that exhibit a poor correlation (less than the defined r value) with the dominant ImSig module. This step of feature selection is recommended to enrich the prediction of relative abundance of immune cells.
#' @param exp Dataframe of transcriptomic data (natural scale) containing genes as rows and samples as columns. Note: Gene names should be set as row names and duplicates are not allowed. Missing values are not allowed within the expression matrix. Check example- head(example_data): \code{\link{example_data}}.
#' @param r Use a value between 0 and 1. Default is 0.6. This is a user defined correlation cut-off to perform feature selection. To get an idea of what cut-off to use check the results of (\code{\link{gene_stat}}) and choose a cut-off that displays high median correlation and maintains a high proportion of genes after feature selection.
#' @return Returns a list of 'feature selected' genes based on the set r value.
#' @examples
#' feature_select (exp = example_data, r = 0.7)
#' @export

feature_select <- function(exp, r = 0.6){
  exp <- pp_exp(exp)
  sig <- pp_sig(exp)
  fg_all <- list()
  for (i in levels(sig$cell)){
    sig_subb <- sig[sig$cell %in% i,]
    exp_subb <- exp[as.character(sig_subb$gene),]
    cor_data <- corr_matrix (exp_subb, r)
    fg_subb <- apply(cor_data,1,function(x)!all(x==0))
    fg_subb <- row.names(cor_data[fg_subb,fg_subb])
    fg_all <- c(fg_all, fg_subb)
  }
  fg <- do.call(c, fg_all)
  return(fg)
}
