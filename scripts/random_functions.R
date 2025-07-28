trail_slash <- function(path){
  if(str_extract(path,"(.{1})$",group = 1)!="/"){
    return(paste0(path,"/"))
  } else {
    return(path)
  }
}

recalculate_p_intersect <- function(results=NULL,DE_methods=DE_methods){
  
  p <- results %>%
    # get all adj_values
    dplyr::select(ends_with("adj_p")) %>%
    # get only the DEmethods
    dplyr::select(starts_with(DE_methods))
  
  # extract max p values from selected methods
  p_intersect <- do.call(pmax, c(p, na.rm = T))
  
  # replace p_intersect
  return_results <- results %>%
    dplyr::select(-p_intersect) %>%
    cbind(p_intersect)
  
  return(return_results)
  
}

add_gene_col <- function(data, gene_col, from, to, species, multi_map="all"){
  
  from_ensembl <- toupper(from) %in% c("ENSEMBL","ENS")
  from_symbol <- toupper(from) %in% c("SYMBOL","SYM","HGNC","NAME")
  from_entrez <- toupper(from) %in% c("ENTREZ")
  to_ensembl <- toupper(to) %in% c("ENSEMBL","ENS")
  to_symbol <- toupper(to) %in% c("SYMBOL","SYM","HGNC","NAME")
  to_entrez <- toupper(to) %in% c("ENTREZ")
  
  if(from_ensembl){
    if(to_entrez){
      map <- as.data.frame(org.Hs.egENSEMBL)
      colnames(map) <- c("entrez_ID","ensembl_ID")
      out_data <- data %>%
        left_join(map, by=c(setNames("ensembl_ID",gene_col)))
    }
  }
  if(multi_map=="all"){
    return(out_data)
  }
  
  
  
  
  
}
