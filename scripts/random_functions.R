trail_slash <- function(path){
  if(str_extract(path,"(.{1})$",group = 1)!="/"){
    return(paste0(path,"/"))
  } else {
    return(path)
  }
}