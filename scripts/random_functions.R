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