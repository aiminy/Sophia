#' Title
#'
#' @param input_dir
#'
#' @return
#' @export
#'
#' @examples
Function2ReadCEL_or_cel<-function(input_dir){
  affybatch <- read.celfiles(list.celfiles(setwd(input_dir)))
  eset <- rma(affybatch)
  re<-list(aff_data=affybatch,ed_norma=eset)
  return(re)
}
