#' @title Network graph of ImSig genes
#' @description A Network visualization displays undirected graph structures and highlights the relationships between entities. The nodes are ImSig genes and the edges represent the correlation between them. The nodes are coloured based on cell type. Try using a correlation cut-off of '0' to get a complete picture.
#' @param exp Dataframe of transcriptomic data (natural scale) containing genes as rows and samples as columns. Note: Gene names should be set as row names and duplicates are not allowed. Missing values are not allowed within the expression matrix. Check example- head(example_data): \code{\link{example_data}}.
#' @param r Use a value between 0 and 1. Default is 0.6. This is a user defined correlation cut-off to perform feature selection (\code{\link{feature_select}}). Feature selection aids to enrich the prediction of relative abundance of immune cells by filtering off poorly correlated ImSig genes. To get an idea of what cut-off to use check the results of (\code{\link{gene_stat}}) and choose a cut-off that displays high median correlation and maintains a high proportion of genes after feature selection.
#' @param pt.cex expansion factor(s) for the points.
#' @param cex character expansion factor relative to current par("cex"). Used for text, and provides the default for pt.cex.
#' @param inset inset distance(s) from the margins as a fraction of the plot region when legend is placed by keyword.
#' @param x.intersp character interspacing factor for horizontal (x) spacing.
#' @param vertex.size Node size of network graph
#' @param vertex.label Add gene names to the network graph. Default set to NA.
#' @param layout Layout algorithm to be used for building network. Default set to force-directed layout algorithm by Fruchterman and Reingold. Read documentation of 'igraph' for other available algorithms.
#' @return Network graph
#' @seealso \code{\link{feature_select}}
#' @import RColorBrewer
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph layout_with_fr
#' @import graphics
#' @examples
#' plot_network (exp = example_data, r = 0.7)
#' @export

plot_network <- function(exp, r=0.6, pt.cex = 2, cex = 1, inset = 0, x.intersp=2, vertex.size = 3, vertex.label= NA, layout = layout_with_fr){
  exp <- pp_exp(exp)
  sig <- pp_sig(exp)
  fg <- feature_select(exp, r)
  cor_data <- corr_matrix(exp, r)
  cor_data <- cor_data[fg,fg]
  network <- graph_from_adjacency_matrix(cor_data, weighted=T, mode="undirected", diag=F)
  # Colour the nodes
  gene <- sig$gene
  sig_sub <- subset(sig, gene %in% fg)
  sig_sub <- sig_sub[match(fg, sig_sub$gene),]
  coul <- brewer.pal(nlevels(as.factor(sig_sub$cell)), "Set3")
  # Map the color to cell types
  my_color <- coul[sig_sub$cell]
  par(mfrow=c(1,2))
  plot(network,
       vertex.size = vertex.size,
       vertex.label= vertex.label,
       layout = layout,
       vertex.color=my_color)
  plot.new()
  legend('center', x.intersp=x.intersp, inset = inset, legend=levels(sig_sub$cell), col = coul, pch=20 , pt.cex = pt.cex, cex = cex, text.col="black", box.lty=0)
}

