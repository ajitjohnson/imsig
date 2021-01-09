#' @title Survival analysis based on relative abundance of immune infiltration estimated by ImSig
#' @description Patients are split into two groups based on their immune cell abundance (median aundance value) and a regular survival analyis is carried out.
#' @param exp Dataframe of transcriptomic data (natural scale) containing genes as rows and samples as columns. Note: Gene names should be set as row names and duplicates are not allowed. Missing values are not allowed within the expression matrix. Check example- head(example_data): \code{\link{example_data}}.
#' @param cli Clinical metadata containting the event data (dead or alive) and time to event data. Samples names should be in rownames and same as that in the expression file. Check head() of \code{\link{example_cli}} for an example clinical data.
#' @param time Column name of time-to-event parameter.
#' @param status Column name of event (dead or alive) parameter.
#' @param r Use a value between 0 and 1. Default is 0.6. This is a user defined correlation cut-off to perform feature selection (\code{\link{feature_select}}). Feature selection aids to enrich the prediction of relative abundance of immune cells by filtering off poorly correlated ImSig genes. To get an idea of what cut-off to use check the results of (\code{\link{gene_stat}}) and choose a cut-off that displays high median correlation and maintains a high proportion of genes after feature selection.
#' @return Hazard Ratio
#' @examples \donttest{
#' survival = imsig_survival (exp = example_data, cli = example_cli)
#' head(survival)
#' }
#' @import survival
#' @import stats
#' @seealso \code{\link{feature_select}}, \code{\link{example_data}}, \code{\link{example_cli}}
#' @export

imsig_survival <- function(exp, cli, time = 'time', status= 'status', r = 0.6){
  sig <- sig
  cell <- imsig(exp,r)
  cell_cli <- merge(cell, cli, by = 'row.names')
  row.names(cell_cli) <- cell_cli[,1]
  cell_cli <- cell_cli[,-1]
  HR <- data.frame()
  for (i in 1: length(levels(sig$cell))){
    cell_ordered <- cell_cli[unlist(sort.list(cell_cli[,i])),]
    cell_ordered$group <- ifelse(cell_ordered[,i] <= median(cell_cli[,i]), "low", "high")
    cox <- survival::coxph(Surv(time, status) ~ group, cell_ordered)
    x <- summary(cox)[8]
    HR.OS <- data.frame(log2(x[[1]][1]))
    HR.UP <- data.frame(log2(x[[1]][4]))
    HR.LOW <- data.frame(log2(x[[1]][3]))
    HR.PVAL <- data.frame(summary(cox)[[9]][[3]][1])
    Hazard_Ratio <- cbind(HR.OS, HR.UP, HR.LOW, HR.PVAL)
    HR <- rbind(HR, Hazard_Ratio)
  }
  row.names(HR) <- levels(sig$cell)
  colnames(HR) <- c("Log2-Hazard ratio of low expression", "95% CI upper limit", "95% CI lower limit", "P value")
  return(HR)
}
