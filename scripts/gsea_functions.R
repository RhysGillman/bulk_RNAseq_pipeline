#gsea_functions.R

#' @param example description
#' @return TRUE
#' @export

run_gsea <- function(results=NULL,
                     name="unnamed_comparison",
                     gene_sets=".*\\.all", #match hierarchical terms from https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=H
                     output_dir=NULL, 
                     gsea_path=NULL,
                     gmt_dir=NULL,
                     nperm=1000,
                     set_max=500,
                     set_min=15,
                     threads=6,
                     threads_per_job=6,
                     log="logs/gsea.log") {
  
  # prepare multithreading
  
  if(threads >= threads_per_job){
    jobs <- threads %/% threads_per_job
  }else{
    jobs=1
    threads_per_job=threads
  }
  
  cl <- makeCluster(jobs, outfile=log)
  clusterExport(cl, varlist = "threads_per_job",
                envir   = environment())
  # tell every worker to see only threads_per_job cores
  clusterEvalQ(cl, Sys.setenv(JAVA_TOOL_OPTIONS = paste0("-XX:ActiveProcessorCount=", threads_per_job)))
  registerDoParallel(cl)
  
  
  source("scripts/random_functions.R")
  
  output_dir <- trail_slash(output_dir)
  gmt_dir <- trail_slash(gmt_dir)
  
  
  dir.create(paste0(output_dir,name))
  
  # Prepare ranked gene list with correct column names
  ranked_genes <- results %>%
    as.data.frame() %>%
    dplyr::filter(!is.na(LogFC), !is.na(p_intersect), !is.na(symbol), symbol != "") %>%
    dplyr::mutate(rank_metric = -log10(p_intersect) * sign(LogFC)) %>%
    dplyr::arrange(desc(rank_metric)) %>%
    dplyr::select(symbol, rank_metric)
  
  message(paste0(nrow(ranked_genes),"/",nrow(results)," results rows included in analysis\n"))
  
  rnk_file <- paste0(output_dir,name,".rnk")
  
  # Write ranked file
  write_tsv(ranked_genes, rnk_file, col_names = F, quote = "none")
  
  # Run GSEA for each gene set collection
  
  # get GMT files based on regex provided
  gmt_files <- foreach(gene_set=gene_sets, .combine = "c") %do% {
    tmp_files <- list.files(gmt_dir, pattern = paste0(gene_set,".*symbols.gmt"), full.names = T)
    # abort script if no matches
    if(length(tmp_files)<1){
      rlang::abort(paste0("\n\nERROR! Your gene set selection: (", gene_set, ") didn't match any files\n",
                   "If you think the gene set regex is correct, check the ref_gmt_dir"
                   ))
    }else{
      tmp_files
    }
    
  }
  
  message(paste0("Running GSEA in ", jobs, " workers with ", threads_per_job, "threads each\n"))
  
  foreach(gmt=gmt_files, .packages = c("stringr","rlang"),
          .export   = c("gsea_path","rnk_file","name","nperm","set_max","set_min","output_dir")
          ) %dopar% {
    # get a gene set identifier
    gmt_name <- str_extract(basename(gmt), "(^.*)\\.v[0-9]{4}.*", group = 1)
  
    gsea_cmd <- paste(
      paste0(gsea_path,"gsea-cli.sh"), "GSEAPreranked",
      "-gmx", gmt,
      "-norm meandiv",
      "-nperm",nperm,
      "-rnk", rnk_file,
      "-scoring_scheme weighted",
      "-rpt_label", paste0(gmt_name),
      "-create_svgs false",
      "-make_sets true",
      "-plot_top_x 40",
      "-rnd_seed timestamp",
      "-set_max", set_max,
      "-set_min", set_min,
      "-zip_report false",
      "-out", paste0(output_dir,name)
    )
    
    message(paste0("Running GSEA for", name, "with", gmt_name, "\n"))
    status <- system(gsea_cmd)
    
    if (status != 0) {
      rlang::abort(paste("GSEA failed for", gmt_name, "exit code", status))
    }
  }
  stopCluster(cl)
  return(TRUE)
}


read_gsea_results <- function(gsea_dir) {
  # First, get all comparison directories
  comparison_dirs <- list.dirs(gsea_dir, full.names = TRUE, recursive = FALSE) %>%
    # keep only the correct directories
    grep(pattern="/.*-.*$", value = T)
    
  
  all_results <- list()
  
  # Process each comparison directory
  for (comp_dir in comparison_dirs) {
    comparison_name <- basename(comp_dir)
    
    # Find all GSEA result directories within this comparison
    gsea_result_dirs <- list.dirs(comp_dir, full.names = TRUE, recursive = FALSE)
    gsea_result_dirs <- gsea_result_dirs[grepl("GseaPreranked", gsea_result_dirs)]
    
    for (result_dir in gsea_result_dirs) {
      # Extract database info from directory name
      #dir_name <- basename(result_dir)
      
      # Parse the directory name (format: database.comparison.GseaPreranked.timestamp)
      #parts <- strsplit(dir_name, "\\.")[[1]]
      #if (length(parts) >= 3) {
      db <- str_extract(basename(result_dir), "(.*)\\.GseaPreranked", group = 1)
      
      # Find report files
      report_files <- list.files(result_dir, 
                                 pattern = "gsea_report_.*\\.tsv$", 
                                 full.names = TRUE)
      
      for (report_file in report_files) {
        if (!file.exists(report_file) || file.size(report_file) == 0) next
        
        # Determine if this is positive or negative enrichment
        is_pos <- grepl("pos", basename(report_file))
        
        # Read the report
        result <- tryCatch({
          read.delim(report_file, check.names = FALSE, stringsAsFactors = FALSE)
        }, error = function(e) {
          message("Error reading file ", report_file, ": ", e$message)
          return(NULL)
        })
        
        if (!is.null(result) && nrow(result) > 0) {
          result$database <- db
          result$comparison <- comparison_name
          result$direction <- ifelse(is_pos, "positive", "negative")
          
          # Remove columns that are entirely NAs
          result <- dplyr::select(result, where(function(x) !all(is.na(x))))
          
          # Create unique key
          key <- paste(db, comparison_name, ifelse(is_pos, "pos", "neg"), sep = "_")
          all_results[[key]] <- result
        }
        
        
        
        
      }
    }
  }
  
  if (length(all_results) > 0) {
    # Define numeric columns
    numeric_columns <- c("NES", "NOM p-val", "FDR q-val", "FWER p-val", "ES", "SIZE", "RANK AT MAX")
    
    # Clean up data
    for (i in seq_along(all_results)) {
      # Convert numeric columns
      for (col in numeric_columns) {
        if (col %in% names(all_results[[i]])) {
          all_results[[i]][[col]] <- as.numeric(as.character(all_results[[i]][[col]]))
        }
      }
      
      # Remove rows with NA NES values
      if ("NES" %in% names(all_results[[i]])) {
        all_results[[i]] <- all_results[[i]][!is.na(all_results[[i]]$NES), ]
      }
      
    }
    
    # Combine all results
    combined_results <- dplyr::bind_rows(all_results)
    return(combined_results)
  } else {
    message("No GSEA results found in ", gsea_dir)
    return(NULL)
  }
}

create_gsea_dotplot <- function(data, comparison_name, database_name, output_file) {
  # Filter data for specific comparison and database
  plot_data <- data %>%
    filter(comparison == comparison_name, database == database_name) %>%
    filter(abs(NES) > 1, `FDR q-val` < 0.25) %>%
    arrange(`FDR q-val`) %>%
    slice_head(n = 30)  # Top 30 pathways
  
  if (nrow(plot_data) == 0) {
    message("No significant pathways for ", comparison_name, " - ", database_name)
    return(NULL)
  }
  
  # Prepare data for plotting
  plot_data <- plot_data %>%
    mutate(
      # Truncate long pathway names
      NAME_display = str_trunc(NAME, 50, "right"),
      # Calculate -log10(FDR) for sizing
      neg_log10_fdr = -log10(pmax(`FDR q-val`, 1e-10)),
      neg_log10_fdr = pmin(neg_log10_fdr, 10),
      # Ensure unique names
      NAME_display = make.unique(NAME_display)
    ) %>%
    arrange(NES)
  
  # Convert to factor for proper ordering
  plot_data$NAME_display <- factor(plot_data$NAME_display, 
                                   levels = plot_data$NAME_display)
  
  # Create plot
  p <- ggplot(plot_data, aes(x = NES, y = NAME_display, 
                             size = neg_log10_fdr, 
                             color = `FDR q-val`)) +
    geom_point() +
    scale_color_gradientn(
      colors = colorRampPalette(c("darkred", "red", "orange", "yellow"))(100),
      trans = "log10",
      name = "FDR q-val",
      guide = guide_colorbar(reverse = TRUE)
    ) +
    scale_size_continuous(range = c(2, 8), name = "-log10(FDR)") +
    labs(
      title = paste("GSEA Results:", comparison_name),
      subtitle = paste("Database:", database_name),
      x = "Normalized Enrichment Score (NES)",
      y = "Pathway"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "right",
      panel.grid.major.y = element_line(size = 0.2),
      panel.grid.minor = element_blank()
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, linewidth=2)
  
  # Save plot
  ggsave(output_file, p, width = 10, height = 8, dpi = 300)
  message("Saved dotplot: ", output_file)
  
  return(p)
}

create_gsea_network <- function(data, comparison_name, database_name, 
                                gmt_dir, output_file, gsea_base_dir = NULL) {
  # Filter data
  plot_data <- data %>%
    filter(comparison == comparison_name, database == database_name) %>%
    filter(abs(NES) > 1, `FDR q-val` < 0.25) %>%
    slice_max(abs(NES), n = 30)  # Top 20 pathways for network
  
  if (nrow(plot_data) < 3) {
    message("Not enough pathways for network plot: ", comparison_name, " - ", database_name)
    return(NULL)
  }
  
  pathway_names <- plot_data$NAME
  
  # Read gene sets from GMT file
  gmt_file <- file.path(gmt_dir, paste0(database_name, ".all.v2024.1.Mm.symbols.gmt"))
  
  # If GMT file not found in the standard location, try to find it in GSEA results
  if (!file.exists(gmt_file) && !is.null(gsea_base_dir)) {
    # Look for the GMT file in the comparison directory
    comp_dir <- file.path(gsea_base_dir, comparison_name)
    gsea_dirs <- list.dirs(comp_dir, full.names = TRUE, recursive = FALSE)
    gsea_dirs <- gsea_dirs[grepl(database_name, gsea_dirs)]
    
    if (length(gsea_dirs) > 0) {
      # Look for GMT files in the GSEA output directory
      gmt_files_found <- list.files(gsea_dirs[1], pattern = "\\.gmt$", 
                                    full.names = TRUE, recursive = TRUE)
      if (length(gmt_files_found) > 0) {
        gmt_file <- gmt_files_found[1]
      }
    }
  }
  
  if (!file.exists(gmt_file)) {
    message("GMT file not found: ", gmt_file)
    # Try alternate naming conventions
    alt_gmt_file <- file.path(gmt_dir, paste0(database_name, ".all.v2023.2.Mm.symbols.gmt"))
    if (file.exists(alt_gmt_file)) {
      gmt_file <- alt_gmt_file
    } else {
      return(NULL)
    }
  }
  
  # Parse GMT file
  gene_sets <- list()
  gmt_lines <- readLines(gmt_file)
  
  for (line in gmt_lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      pathway_name <- parts[1]
      if (pathway_name %in% pathway_names) {
        gene_sets[[pathway_name]] <- parts[-c(1,2)]
      }
    }
  }
  
  if (length(gene_sets) < 2) {
    message("Not enough gene sets found for network")
    return(NULL)
  }
  
  # Calculate Jaccard similarity between pathways
  pathway_pairs <- combn(names(gene_sets), 2, simplify = FALSE)
  edges <- data.frame()
  
  for (pair in pathway_pairs) {
    genes1 <- gene_sets[[pair[1]]]
    genes2 <- gene_sets[[pair[2]]]
    
    intersection <- length(intersect(genes1, genes2))
    union_size <- length(union(genes1, genes2))
    
    if (union_size > 0) {
      jaccard <- intersection / union_size
      
      if (jaccard > 0.1) {  # Threshold for edge creation
        edges <- rbind(edges, data.frame(
          from = pair[1],
          to = pair[2],
          weight = jaccard
        ))
      }
    }
  }
  
  if (nrow(edges) == 0) {
    message("No significant connections found for network")
    return(NULL)
  }
  
  # Create nodes data frame
  nodes <- plot_data %>%
    select(NAME, NES, `FDR q-val`) %>%
    rename(name = NAME) %>%
    mutate(
      size = 5 + 3 * abs(NES),
      color = NES
    )
  
  
  # Create graph
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  
  # Create plot
  set.seed(123)  # For reproducibility
  
  # Create layout first
  layout <- create_layout(g, layout = "fr")
  
  # Create the plot with fixed aesthetics mapping
  p <- ggraph(layout) +
    geom_edge_link(aes(width = weight), alpha = 0.3, color = "grey50") +
    geom_node_point(aes(size = size, color = color)) +
    geom_node_text(aes(label = name), 
                   repel = TRUE, size = 3, max.overlaps = 20) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                          midpoint = 0, name = "NES") +
    scale_size_continuous(range = c(4, 12), name = "Effect size") +
    scale_edge_width_continuous(range = c(0.5, 3), name = "Similarity") +
    labs(
      title = paste("GSEA Network:", comparison_name),
      subtitle = paste("Database:", database_name)
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # Save plot
  ggsave(output_file, p, width = 12, height = 10, dpi = 300)
  message("Saved network plot: ", output_file)
  
  return(p)
}
