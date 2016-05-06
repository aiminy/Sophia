#' Title
#'
#' @param data.ER.PR.sample.info
#' @param data.set1.normalized
#' @param sample_type
#'
#' @return
#' @export
#'
#' @examples
MapSample2CELfile<-function(data.ER.PR.sample.info,data.set1.normalized,sample_type){

  cat(dim(data.ER.PR.sample.info),dim(data.set1.normalized),sample_type,"\n")

  num.sample1=length(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==1),][,1])
  cel_file_index_4_sample<-array()

  for(i in 1:num.sample1) {
    tmp<-grep(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==sample_type),][i,1],colnames(data.set1.normalized))
    cel_file_index_4_sample<-c(cel_file_index_4_sample,tmp)
  }
  cel_file_index_4_sample<-cel_file_index_4_sample[-1]

  colnames(data.set1.normalized)[cel_file_index_4_sample]

  re<-data.set1.normalized[,cel_file_index_4_sample]

  return(re)

}
