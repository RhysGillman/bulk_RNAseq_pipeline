#' Generate basic visualisations for RNAseq data
#'
#'
#' @param data Input dataframe
#' @param labels Column with point labels (eg. gene_ID)
#' @param pcol Column with pvalues
#' @param l2fccol Column with l2fc values
#' @param rankcol Column with values to rank top N genes
#' @param invert_rank Use if rankcol implies lower values=better
#' @param topn Top N points to label
#' @param p_thresholds Thresholds of p-values to distinguish colours 
#' @param title Plot Title
#' @return ggplot Plots
#' @export


generate_volcano_plot <- function(data=NULL,
                                   labels="gene_ID",
                                   pcol=NULL,
                                   l2fccol=NULL,
                                   rank_col=NULL,
                                   invert_rank=FALSE,
                                   topn=20,
                                   p_thresholds=c(0.05,0.01),
                                   title=NULL){
  
  results <- data[c(labels,pcol,l2fccol,rank_col)]
  colnames(results) <- c("gene_ID","p","lf2c","stat")
  
  p_thresholds <- sort(p_thresholds, decreasing = T)
  
  results <- results %>%
    # colour groups
    mutate(group=ifelse(p > p_thresholds[1], paste0("p > ",p_thresholds[1]),
                        ifelse(p > p_thresholds[2], paste0("p < ",p_thresholds[1]),
                               paste0("p < ",p_thresholds[2])
                        ))
    ) %>%
    # plot gene labels
    { if (invert_rank) arrange(., stat) else arrange(., desc(stat)) } %>%
    mutate(label=ifelse(row_number() <= topn, gene_ID, NA))
  
  p <- ggplot(results,aes(x = lf2c, y = -log10(p), colour = group, label=label)) +
    geom_point() +
    scale_colour_manual(values = c("black","limegreen","blue"), breaks = c(paste0("p > ",p_thresholds[1]),
                                                                           paste0("p < ",p_thresholds[1]),
                                                                           paste0("p < ",p_thresholds[2])
                                                                           )) +
    guides(colour=guide_legend(title = "p-value Cutoff")) +
    geom_text_repel(size = 3, max.overlaps = Inf, show.legend = F) +
    ylab("-log10(Adjusted p Value)") +
    xlab("log2(Fold Change in Gene Expression)") +
    ggtitle(title) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)

}

generate_PCA_plot <- function(counts, metadata, title) {
  
  # data transformations
  md <- prep(t(counts), scale = "none", centre = FALSE)
  pc <- pca(md, method="svd", center=FALSE, nPcs=ncol(counts))
  #var_3 <- R2cum(pc)[3] # accumulated variance
  pc_1 <- round(pc@R2[1]*100, 2)
  pc_2 <- round(pc@R2[2]*100, 2)
  
  df <- as.data.frame(scores(pc)) %>%
    rownames_to_column("file") %>%
    mutate(sample_ID=gsub("_summarised_counts.txt","",file)) %>%
    left_join(metadata %>% dplyr::select(sample_ID,condition) %>% unique(), by="sample_ID")
  
  ggplot(df, aes(x=PC1, y=PC2, colour=condition, label=sample_ID)) +
    geom_point(size=3) +
    xlab(paste0("PC1 - ", round(pc_1, 2), "%")) +
    ylab(paste0("PC2 - ", round(pc_2, 2), "%")) +
    geom_text_repel(show.legend = F) +
    guides(colour=guide_legend(title = "Condition")) +
    ggtitle(title) +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))
  
}
