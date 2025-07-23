#enrichment_analyses.R


# Load Required Packages
message("Loading required packages...")
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
suppressPackageStartupMessages (library(doParallel, quietly = T))

# Handling input arguments
option_list = list(
  make_option(c("-m", "--metadata"), type="character", default=NULL, 
              help="Path to metadata file", metavar ="MetaData"),
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Path to consensusDE directory", metavar ="InputPath"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Path to output directory for results files.", metavar ="OutputPath"),
  make_option(c("-t", "--threads"), type="numeric", default=4, 
              help="Number of threads", metavar ="Threads"),
  make_option(c("-g", "--gene-sets"), type="character", default=".*\\.all", 
              help="Regex defining gene sets to be studied using gsea", metavar ="GeneSets"),
  make_option(c("-G", "--gsea-path"), type="character", default=NULL, 
              help="Path to GSEA installation", metavar ="GSEAPath"),
  make_option(c("-r", "--ref-gmt"), type="character", default=NULL, 
              help="Path to directory of GMT gene set files", metavar ="GMT"),
  make_option(c("-n", "--gsea-nperm"), type="numeric", default=1000, 
              help="GSEA nperm", metavar ="nPerm"),
  make_option(c("-z", "--gsea-set-max"), type="numeric", default=500, 
              help="GSEA set_max", metavar ="SetMax"),
  make_option(c("-y", "--gsea-set-min"), type="numeric", default=15, 
              help="GSEA set_min", metavar ="SetMax"),
  make_option(c("-l", "--log"), type="character", default=15, 
              help="Path to log file for gsea parallel workers", metavar ="Log")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
metadata_path <- opt$metadata
in_path <- opt$input
out_path <- opt$output
threads <- opt$threads
gene_sets <- opt$`gene-sets`
gsea_path <- opt$`gsea-path`
gmt_dir <- opt$`ref-gmt`
gsea_nperm <- opt$`gsea-nperm`
gsea_set_max <- opt$`gsea-set-max`
gsea_set_min <- opt$`gsea-set-min`
gsea_log <- opt$log

# Sourcing Scripts
source("scripts/random_functions.R")
source("scripts/gsea_functions.R")

# paths need trailing /
in_path <- trail_slash(in_path)
out_path <- trail_slash(out_path)
gsea_path <- trail_slash(gsea_path)

# Run GSEA

result_paths <- list.files(paste0(in_path,"results"), pattern = "*combined_results_recal.tsv", full.names = T)

message(paste0("Detected the following results tables: \n", paste(result_paths, collapse = "\n")))

for(comparison in result_paths){
  
  message(paste0("Beginning GSEA for comparison: ", gsub("_combined_results_recal.tsv","",basename(comparison))))
  
  merged_results <- fread(comparison) %>% as.data.frame()

  run_gsea(results=merged_results,
           name=gsub("_combined_results_recal.tsv","",basename(comparison)),
           gene_sets=gene_sets,
           output_dir=out_path, 
           gsea_path=gsea_path,
           gmt_dir=gmt_dir,
           nperm=gsea_nperm,
           set_max=gsea_set_max,
           set_min=gsea_set_min,
           threads=threads,
           threads_per_job=6,
           log=gsea_log
           )
}

# Get GSEA Results

all_results <- read_gsea_results(out_path)


# Generate Plots

dir.create(paste0(out_path,"plots"))

for(comp in unique(all_results$comparison)){
  
  for(db in all_results %>% filter(comparison==comp) %>% pull(database) %>% unique()){
    
    create_gsea_dotplot(data = all_results,
                        comparison_name = comp,
                        database_name = db,
                        output_file = paste0(out_path,"plots/", comp,"-", db,".png"))
    
  }
}


create_gsea_network(data = ,
                    comparison_name = ,
                    database_name = ,
                    gmt_dir = ,
                    output_file = ,
                    gsea_base_dir = )

