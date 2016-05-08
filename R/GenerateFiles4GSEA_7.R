#Use all data to perform GSEA analysis
#' Title
#'
#' @param ed.set123.frma.cancer
#' @param cel.file.name.key.set.symbol.st134
#' @param cutoff
#' @param out_dir
#' @param out_file_name
#' @param subtypeA
#' @param subtypeB
#'
#' @return
#' @export
#'
#' @examples
GenerateFiles4GSEA<-function(ed.set123.frma.cancer,cel.file.name.key.set.symbol.st134,cutoff,out_dir,out_file_name,subtypeA,subtypeB){

  index.st4<-which(colnames(ed.set123.frma.cancer) %in% cel.file.name.key.set.symbol.st134[which(cel.file.name.key.set.symbol.st134[,5]==subtypeA),2])
  index.st1<-which(colnames(ed.set123.frma.cancer) %in% cel.file.name.key.set.symbol.st134[which(cel.file.name.key.set.symbol.st134[,5]==subtypeB),2])

  ed.set123.frma.cancer.st41<-cbind(ed.set123.frma.cancer[,index.st4],ed.set123.frma.cancer[,index.st1])
  cat(subtypeA,"\t",length(index.st4),subtypeB,"\t",length(index.st1),"\t",dim(ed.set123.frma.cancer.st41),"\n")

  sink(paste0(out_dir,cutoff,"_",subtypeA,"vs",subtypeB,"_phenotype.cls"))
  cat(dim(ed.set123.frma.cancer.st41)[2],2,1,"\n")
  cat("#",subtypeA,subtypeB,"\n")
  cat(c(rep(0,length(index.st4)),rep(1,length(index.st1))),"\n")
  sink()

  ed.set123.frma.cancer.st41.reformat<-cbind(rownames(ed.set123.frma.cancer.st41),rep("NA",dim(ed.set123.frma.cancer.st41)[1]),ed.set123.frma.cancer.st41)
  colnames(ed.set123.frma.cancer.st41.reformat)[c(1,2)]=c("NAME","Description")

  cat("#1.2","\n",file=paste0(out_dir,cutoff,"_",out_file_name), sep="",append=TRUE)
  cat(dim(ed.set123.frma.cancer.st41)[1]," ",dim(ed.set123.frma.cancer.st41)[2], "\n",file=paste0(out_dir,cutoff,"_",out_file_name), sep="",append=TRUE)
  write.table(ed.set123.frma.cancer.st41.reformat,file = paste0(out_dir,cutoff,"_",out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = TRUE,append=TRUE)

}
