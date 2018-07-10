#' @title Plot relative abundance of immune cells
#' @description Barplots of relative abundance of immune cells across samples.The order of the samples are the same as that of \code{\link{imsig}}.
#' @param exp Dataframe of transcriptomic data (natural scale) containing genes as rows and samples as columns. Note: Gene names should be set as row names and duplicates are not allowed. Missing values are not allowed within the expression matrix. Check example- head(example_data): \code{\link{example_data}}.
#' @param r Use a value between 0 and 1. Default is 0.6. This is a user defined correlation cut-off to perform feature selection (\code{\link{feature_select}}). Feature selection aids to enrich the prediction of relative abundance of immune cells by filtering off poorly correlated ImSig genes. To get an idea of what cut-off to use check the results of (\code{\link{gene_stat}}) and choose a cut-off that displays high median correlation and maintains a high proportion of genes after feature selection.
#' @return ggplot
#' @import ggplot2
#' @import gridExtra
#' @examples
#' plot_abundance (exp = example_data, r = 0.7)
#' @seealso \code{\link{feature_select}}, \code{\link{example_data}}
#' @export

plot_abundance <- function(exp, r = 0.6){
  cell <- imsig(exp, r)
  cell$samples <- row.names(cell)
  cell$samples <- factor(cell$samples, levels = cell$samples)
  plots = lapply(1:(ncol(cell)-1), function(x) ggplot(cell, aes(x = cell$samples, y = cell[,x]))
                 + geom_bar(stat = "identity") + theme_classic() +
                   theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank())+
                   ggtitle(colnames(cell)[x]))
  do.call(grid.arrange,  plots)
}
