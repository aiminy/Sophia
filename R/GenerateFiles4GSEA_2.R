#Use all data to perform GSEA analysis
#' Title
#'
#' @param ed.set123.frma.cancer
#' @param out_dir
#' @param out_file_name
#' @param subtypeA
#' @param subtypeB
#'
#' @return
#' @export
#'
#' @examples
GenerateFiles4GSEA<-function(ed.set123.frma.cancer,out_dir,out_file_name,subtypeA,subtypeB){

  index.st4<-which(colnames(ed.set123.frma.cancer) %in% cel.file.name.key.set.symbol.st134[which(cel.file.name.key.set.symbol.st134[,5]==subtypeA),2])
  index.st1<-which(colnames(ed.set123.frma.cancer) %in% cel.file.name.key.set.symbol.st134[which(cel.file.name.key.set.symbol.st134[,5]==subtypeB),2])

  ed.set123.frma.cancer.st41<-cbind(ed.set123.frma.cancer[,index.st4],ed.set123.frma.cancer[,index.st1])
  cat(subtypeA,"\t",length(index.st4),subtypeB,"\t",length(index.st1),"\t",dim(ed.set123.frma.cancer.st41),"\n")

  ed.set123.frma.cancer.st41.reformat<-cbind(rownames(ed.set123.frma.cancer.st41),rep("NA",dim(ed.set123.frma.cancer.st41)[1]),ed.set123.frma.cancer.st41)
  colnames(ed.set123.frma.cancer.st41.reformat)[c(1,2)]=c("NAME","Description")

  write.table(ed.set123.frma.cancer.st41.reformat,file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = TRUE)

}
