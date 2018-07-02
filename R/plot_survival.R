#' @title Forest plot of survial analysis by ImSig
#' @description Patients are split into two groups based on their immune cell abundance (median aundance value) and a regular survival analyis is carried out. Raw values can be obtained from \code{\link{imsig_survival}}.
#' @param exp Dataframe of transcriptomic data (natural scale) containing genes as rows and samples as columns. Note: Gene names should be set as row names and duplicates are not allowed. Missing values are not allowed within the expression matrix. Check example- head(example_data): \code{\link{example_data}}.
#' @param cli Clinical metadata containting the event data (dead or alive) and time to event data. Samples names should be in rownames and same as that in the expression file. Check head() of \code{\link{example_cli}} for an example clinical data.
#' @param time Column name of time-to-event parameter.
#' @param status Column name of event (dead or alive) parameter.
#' @param r Use a value between 0 and 1. Default is 0.6. This is a user defined correlation cut-off to perform feature selection (\code{\link{feature_select}}). Feature selection aids to enrich the prediction of relative abundance of immune cells by filtering off poorly correlated ImSig genes. To get an idea of what cut-off to use check the results of (\code{\link{gene_stat}} and choose a cut-off that displays high median correlation and maintains a high proportion of genes after feature selection.
#' @return Forest plot
#' @import ggplot2
#' @examples
#' plot_survival (exp = example_data, r = 0.7, cli = example_cli, time = 'time', status= 'status')
#' @seealso \code{\link{feature_select}}, \code{\link{example_data}}, \code{\link{example_cli}}
#' @export

plot_survival <- function(exp, cli, time = 'time', status= 'status', r = 0.6){
  surv_data <- imsig_survival (exp, cli, time, status, r)
  surv_data$p <- ifelse(surv_data$`P value` <= 0.05, "P < 0.05", "non-significant")
  ggplot(data=surv_data, aes(x=row.names(surv_data), y=surv_data$`Log2-Hazard ratio of low expression`, ymin=surv_data$`95% CI lower limit`, ymax=surv_data$`95% CI upper limit`, colour = as.factor(surv_data$p))) +
    geom_pointrange() +
    geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("Cell types") + ylab("Hazard Ratio (95% CI)") +
    theme_bw()  +
    labs(shape = "", linetype = "") +
    labs(colour = "Hazard Ratio P value")+
    theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill=NA),
          aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12))
}
