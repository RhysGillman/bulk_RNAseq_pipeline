
run_DOSE <- function(results, out, name){
  source("scripts/random_functions.R")
  
  # Prepare ranked gene list with correct column names
  ranked_genes <- results %>%
    as.data.frame() %>%
    dplyr::filter(!is.na(LogFC), !is.na(p_intersect), !is.na(entrez_ID), symbol != "") %>%
    dplyr::mutate(rank_metric = -log10(p_intersect) * sign(LogFC)) %>%
    dplyr::arrange(desc(rank_metric)) %>%
    dplyr::select(entrez_ID,rank_metric) %>%
    unique() %>%
    deframe()
  
  enriched_DO <- gseDO(
    geneList=ranked_genes,
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE,
    seed = FALSE,
    by = "fgsea"
  )
  
  saveRDS(enriched_DO,paste0(out,name,"_DO_results_object.rds"))
  
  png(paste0(out,name,"_DO_dot_plot.png"), 
      width = 10, height = 15, units = "in", res = 300)
  
  dotplot(enriched_DO)
  
  dev.off()
  
  
  
  
}