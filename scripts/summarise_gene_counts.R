#summarise_gene_counts.R


# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(tximport, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))

# Handling input arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Path to directory containing RSEM isoforms.results files", metavar ="Input Path"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Path to output directory. Defaults (NULL) to same as input", metavar ="Output Path"),
  make_option(c("-t", "--tx2gene"), type="character", default=NULL, 
              help="Path to tx2gene map. By default (NULL) will be generated from input files", metavar ="Output Path")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

in_path <- opt$input
out_path <- opt$output
tx2gene_path <- opt$tx2gene


# remove trailing '/' from paths
in_path <- gsub("/$","",in_path)



if(!dir.exists(in_path)){
  stop(paste0("IMPORT FAILED: ",in_path," does not exist"))
}

# get full file paths

input_files <- list.files(in_path, pattern = "\\.isoforms\\.results", full.names = TRUE)

# add sample names
sample_names <- str_extract(basename(input_files), "(^.*)\\.isoforms\\.results", group = 1)
names(input_files) <- sample_names

if(length(input_files)==0){
  stop(paste0("IMPORT FAILED: no '.isoforms.results files in '", in_path,"'"))
}

if(is.null(out_path)){
  out_path=in_path
}

out_path <- gsub("/$","",out_path)

#tx2gene map

if(!is.null(tx2gene_path)){
  tx2gene <- read_tsv(tx2gene, col_names = c("gene_ID","transcript_ID")) %>%
    dplyr::select(transcript_ID,gene_ID) %>%
    unique
}else{
  present <- all(c("gene_id","transcript_id") %in% colnames(fread(input_files[1], nrows = 1)))
  
  if(present){
    message("Making tx2gene from input files")
    tx2gene <- fread(input_files[1], select = c("transcript_id","gene_id")) %>%
      unique()
  }else{
    stop("ERROR: Error, no tx2gene map provided and can't find 'gene_id' and 'transcript_id' in input files")
  }
}

txi_genes <- tximport(input_files,
                      type = "rsem",
                      txIn = T,
                      txOut = F,
                      tx2gene = tx2gene,
                      countsFromAbundance="lengthScaledTPM"
)

counts <- round(txi_genes$counts)

for(i in colnames(counts)){
  message(paste0("Creating File: ",out_path,"/",i,"_summarised_counts.txt"))
  
  write_tsv(counts[,i] %>% as.data.frame() %>% rownames_to_column("gene_ID"), col_names = F,
            paste0(out_path,"/",i,"_summarised_counts.txt"))
}

