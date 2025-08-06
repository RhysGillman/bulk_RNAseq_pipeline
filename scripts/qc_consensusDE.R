#qc_consensusDE.R


# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(tximport, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(consensusDE, quietly = T))
suppressPackageStartupMessages (library(org.Hs.eg.db, quietly = T))
suppressPackageStartupMessages (library(ggrepel, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(patchwork, quietly = T))
suppressPackageStartupMessages (library(pcaMethods, quietly = T))
suppressPackageStartupMessages (library(EDASeq, quietly = T))

# Handling input arguments
option_list = list(
  make_option(c("-m", "--metadata"), type="character", default=NULL, 
              help="Path to metadata file", metavar ="MetaData"),
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Path to consensusDE directory", metavar ="InputPath"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Path to output directory for results files.", metavar ="OutputPath"),
  make_option(c("-d", "--DEmethods"), type="character", default="deseq,edger,voom", 
              help="Comma separated list of DE methods to run", metavar ="OutputPath"),
  make_option(c("-t", "--threads"), type="numeric", default=4, 
              help="Number of threads", metavar ="Threads"),
  make_option(c("-s", "--species"), type="character", default="human", 
              help="Species being studies (human/homo sapiens or mouse/mus musculus", metavar ="Species")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
metadata_path <- opt$metadata
in_path <- opt$input
out_path <- opt$output
DE_methods <- unlist(str_split(opt$DEmethods, ","))
threads <- opt$threads
species <- opt$species

# Functions

source("scripts/visualisations.R")
source("scripts/random_functions.R")

# paths need trailing /
in_path <- trail_slash(in_path)
out_path <- trail_slash(out_path)


metadata <- fread(metadata_path)

CDE_raw <- readRDS(paste0(in_path, "CDE_se_raw.rds"))
#CDE_run <- readRDS(paste0(in_path, "CDE_run.rds"))


# check results 

result_paths <- list.files(paste0(in_path,"results/"), pattern = "*combined_results_recal.tsv", full.names = T)

for(comparison in result_paths){
  
  merged_results <- fread(comparison) %>% as.data.frame()
  
  comparison_name <- gsub("_combined_results_recal.tsv","",basename(comparison))
  
  # Volcano plot using p_intersect
  p1 <- generate_volcano_plot(data = merged_results,
                         labels = "symbol",
                         pcol = "p_intersect",
                         l2fccol = "LogFC",
                         rank_col = "rank_sum",
                         invert_rank = TRUE,
                         topn = 20,
                         p_thresholds = c(0.05,0.01),
                         title = "p Intersect")
  
  p_methods <- foreach(method=DE_methods) %do% {
    
    generate_volcano_plot(data = merged_results,
                          labels = "symbol",
                          pcol = paste0(method,"_adj_p"),
                          l2fccol = "LogFC",
                          rank_col = paste0(method,"_rank"),
                          invert_rank = TRUE,
                          topn = 20,
                          p_thresholds = c(0.05,0.01),
                          title = method)
    
  }
  names(p_methods) <- DE_methods
  
  (p1 / (wrap_plots(p_methods, nrow = 1) + plot_layout(axis_titles = "collect"))) +
    plot_layout(guides = "collect",axis_titles = "collect")
  
  ggsave(paste0(out_path,paste0("plots/",comparison_name,"_all_volcano_plots.png")), height = 30, width=20, units="cm")
}

# PCAs

se_qc <- newSeqExpressionSet(assays(CDE_raw)$counts,
                             phenoData = data.frame(colData(CDE_raw)),
                             row.names = colnames(assays(CDE_raw)$counts))

raw_logCPM <- consensusDE:::check_normalise(se_in = se_qc)

generate_PCA_plot(raw_logCPM, metadata, "Raw")

ggsave(paste0(out_path,"plots/raw_data_pca.png"))

#se_qc_after <- newSeqExpressionSet(assays(CDE_run)$counts,
#                             phenoData = data.frame(colData(CDE_run)),
#                             row.names = colnames(assays(CDE_run)$counts))

#final_logCPM <- consensusDE:::check_normalise(se_in = se_qc_after)

#generate_PCA_plot(final_logCPM, metadata, "Normalised - Corrected")
#ggsave(paste0(out_path,"plots/corrected_normalised_data_pca.png"))


norm_logCPM <- fread(paste0(in_path, "/results/normalised_log_cpm.tsv")) %>%
  column_to_rownames("V1") %>%
  as.matrix()
 
generate_PCA_plot(norm_logCPM, metadata, "Normalised - Corrected")

ggsave(paste0(out_path,"plots/corrected_normalised_data_pca.png"))
