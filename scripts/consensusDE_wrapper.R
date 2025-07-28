#consensusDE_wrapper.R


# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(tximport, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(consensusDE, quietly = T))
suppressPackageStartupMessages (library(org.Hs.eg.db, quietly = T))
suppressPackageStartupMessages (library(ggfortify, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))

# Handling input arguments
option_list = list(
  make_option(c("-m", "--metadata"), type="character", default=NULL, 
              help="Path to metadata file", metavar ="MetaData"),
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Path to directory containing gene-level counts in .txt files", metavar ="InputPath"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Path to output directory for results files.", metavar ="OutputPath"),
  make_option(c("-d", "--DEmethods"), type="character", default="deseq,edger,voom", 
              help="Comma separated list of DE methods to run", metavar ="OutputPath"),
  make_option(c("-a", "--alpha"), type="numeric", default=0.05, 
              help="Alpha (p threshold) for determining DEGs", metavar ="Threads"),
  make_option(c("-g", "--gtf"), type="character", default=NULL, 
              help="Path to gtf", metavar ="GTFPath"),
  make_option(c("-p", "--paired"), type="logical", default=TRUE, 
              help="Is library paired? (Default = TRUE)", metavar ="Paired"),
  make_option(c("-s", "--strandedness"), type="character", default="unstranded", 
              help="Strandedness of library ('unstranded','forward-stranded', 'reverse-stranded'", metavar ="Strandedness"),
  make_option(c("-t", "--threads"), type="numeric", default=4, 
              help="Number of threads", metavar ="Threads")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

metadata_path <- opt$metadata
in_path <- opt$input
out_path <- opt$output
DE_methods <- unlist(str_split(opt$DEmethods, ","))
gtf_path <- opt$gtf
is_paired <- opt$paired
strandedness <- opt$strandedness
threads <- opt$threads


source("scripts/random_functions.R")

# paths need trailing /
in_path <- trail_slash(in_path)
out_path <- trail_slash(out_path)

message(paste0("Results, Plots and Logs will be stored in ", out_path ))

# Create Dirs

dir.create(paste0(out_path,"results"), recursive = T, showWarnings = F)
dir.create(paste0(out_path,"plots"), recursive = T, showWarnings = F)
dir.create(paste0(out_path,"logs"), recursive = T, showWarnings = F)

metadata <- fread(metadata_path)

# work out if data is paired

paired_data <- metadata %>%
  dplyr::select(subject_ID, condition) %>%
  unique()
# check if the nrows of above is equal to n conditions * n subjects
paired_data=ifelse(nrow(paired_data)==(length(unique(paired_data$condition)))*length(unique(paired_data$subject_ID)),
                   "paired","unpaired")

sample_table <- data.frame(file=list.files(in_path)) %>%
  mutate(sample_ID=str_extract(basename(file), "(^.*)_summarised_counts.txt$", group = 1)) %>%
  left_join(metadata, by = "sample_ID") %>%
  dplyr::select(file, group=condition, batch) %>%
  unique()

message("Building consensDE Summarised Object...")

CDE_summarised_object <- buildSummarized(sample_table = sample_table,
                                  htseq_dir = in_path,
                                  read_format = ifelse(is_paired,"paired","single"),
                                  output_log=paste0(out_path,"logs/"),
                                  gtf=gtf_path,
                                  filter=TRUE,
                                  strand_mode=ifelse(
                                    strandedness=="unstranded",0,
                                    ifelse(strandedness=="forward-stranded",1,
                                           ifelse(strandedness=="reverse-stranded",2,NA))
                                  ),
                                  n_cores=threads,
                                  force_build = TRUE)

saveRDS(CDE_summarised_object,paste0(out_path,"CDE_se_raw.rds"))




if(length(unique(sample_table$batch))>1){
  
  message("Detected multiple batches in metadata, running ComBat Seq...")
  
  # batch correction
  # Extract counts from summarized experiment
  counts <- assay(CDE_summarised_object)
  
  # Apply ComBat_seq correction
  batch_corrected <- ComBat_seq(counts = counts, 
                                batch = sample_table$batch, 
                                group = sample_table$group)
  
  # Create a new SummarizedExperiment with the corrected counts
  CDE_summarised_object <- SummarizedExperiment(
    assays = list(counts = batch_corrected),
    colData = colData(CDE_summarised_object),
    rowData = rowData(CDE_summarised_object),
    metadata = metadata(CDE_summarised_object)
  )
  
  saveRDS(CDE_summarised_object,paste0(out_path,"CDE_se_batch_corrected.rds"))
  
  message("Running consensusDE without RUV correction...")
  CDE_run <- multi_de_pairs(summarized = CDE_summarised_object,
                            paired = paired_data,
                            plot_dir = paste0(out_path,"plots/"),
                            output_voom = paste0(out_path, "results/"),
                            output_edger = paste0(out_path, "results/"),
                            output_deseq = paste0(out_path, "results/"),
                            output_combined = paste0(out_path,"results/"),
                            ruv_correct = F)
  
}else{

  message("Running consensusDE with RUV correction...")
  CDE_run <- multi_de_pairs(summarized = CDE_summarised_object,
                            paired = paired_data,
                            plot_dir = paste0(out_path,"plots/"),
                            output_voom = paste0(out_path, "results/"),
                            output_edger = paste0(out_path, "results/"),
                            output_deseq = paste0(out_path, "results/"),
                            output_combined = paste0(out_path,"results/"),
                            ensembl_annotate=org.Hs.eg.db,
                            ruv_correct =T)
  
}


# modify p_intersect

result_paths <- list.files(paste0(out_path,"results/"), pattern = "*combined_results.tsv", full.names = T)

for(comparison in result_paths){
  
  results <- fread(comparison) %>% as.data.frame()
  recal_results <- recalculate_p_intersect(results,DE_methods=DE_methods)
  recal_results_path <- gsub(".tsv","_recal.tsv",comparison)
  write_tsv(recal_results,recal_results_path)
  
}


saveRDS(CDE_run,paste0(out_path,"CDE_run.rds"))


# Summarise Results


result_paths <- list.files(paste0(out_path,"results/"), pattern = "*_recal.tsv", full.names = T)

merged_results <- foreach(comp=result_paths, .combine = "bind_rows") %do% {
  
  fread(comparison) %>% mutate(comparison=gsub("_combined_results_recal.tsv","",basename(comp)))
  
}
  
DEG_alpha=0.05

summary_DEG <- merged_results %>%
  dplyr::select(comparison,ID,p_intersect,deseq_adj_p,edger_adj_p,voom_adj_p) %>%
  pivot_longer(p_intersect:voom_adj_p,names_to="method", values_to = "p_adj") %>%
  filter(p_adj<DEG_alpha) %>%
  group_by(comparison, method) %>%
  summarise(n_degs=n(),.groups="drop") %>%
  pivot_wider(names_from = "method", values_from = "n_degs") %>%
  dplyr::select(comparison,p_intersect,deseq_adj_p,edger_adj_p,voom_adj_p)

# Mark methods that are excluded from p intersect
colnames <- c("comparison","p_intersect","deseq","edger","voom")
excluded <- which(!colnames %in% c("comparison","p_intersect", DE_methods))
colnames[excluded] <- paste0(colnames[excluded]," (Excluded)")

colnames(summary_DEG) <- colnames

write_tsv(summary_DEG, paste0(out_path,"DEG_summary.tsv"))

# Venn diagram?