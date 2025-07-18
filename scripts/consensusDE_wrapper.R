#consensusDE_wrapper.R


# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(tximport, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(consensusDE, quietly = T))
suppressPackageStartupMessages (library(org.Hs.eg.db, quietly = T))
suppressPackageStartupMessages (library(ggfortify, quietly = T))

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
gtf_path <- opt$gtf
is_paired <- opt$paired
strandedness <- opt$strandedness
threads <- opt$threads
DE_methods <- unlist(str_split(opt$DEmethods, ","))

source("scripts/random_functions.R")

# paths need trailing /
in_path <- trail_slash(in_path)
out_path <- trail_slash(out_path)

message(paste0("ConsensusDE is being run using these methods: (", paste0(DE_methods, collapse = ";"),")"))
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
  
  message("Running consensusDE without RUV correction...")
  CDE_run <- multi_de_pairs(summarized = CDE_summarised_object,
                            paired = paired_data,
                            plot_dir = paste0(out_path,"plots/"),
                            output_voom = paste0(out_path,"results/"),
                            output_edger = paste0(out_path,"results/"),
                            output_deseq = paste0(out_path,"results/"),
                            output_combined = paste0(out_path,"results/"),
                            ruv_correct = F)
  
}else{

  message("Running consensusDE with RUV correction...")
  CDE_run <- multi_de_pairs(summarized = CDE_summarised_object,
                            paired = paired_data,
                            plot_dir = paste0(out_path,"plots/"),
                            output_voom     = if ("voom"  %in% DE_methods) paste0(out_path, "results/") else NULL,
                            output_edger    = if ("edger" %in% DE_methods) paste0(out_path, "results/") else NULL,
                            output_deseq    = if ("deseq" %in% DE_methods) paste0(out_path, "results/") else NULL,
                            output_combined = paste0(out_path,"results/"),
                            ensembl_annotate=org.Hs.eg.db,
                            ruv_correct =T)
  
}

saveRDS(CDE_run,paste0(out_path,"CDE_run.rds"))

