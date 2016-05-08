#' Generate PCA plot only
#'
#' @param reorder.column.ed.set123.frma.cancer
#' @param GeneSym.all
#' @param data.set123.raw.with.set.label.and.type.cancer
#' @param output_file
#' @param output_dir
#'
#' @return
#' @export
#'
#' @examples
GeneratePCA<-function(reorder.column.ed.set123.frma.cancer,GeneSym.all,data.set123.raw.with.set.label.and.type.cancer,output_file,output_dir){

  ndata<-reorder.column.ed.set123.frma.cancer

  geneSymbol=GeneSym.all
  tempdata.byGSym =data.frame(ndata, Symbol = geneSymbol)

  colnames(tempdata.byGSym)

  tempdata.byGSym.2<-tempdata.byGSym[-which(is.na(tempdata.byGSym$Symbol)),]

  data.byGSym = ddply(tempdata.byGSym.2, c("Symbol"),function(h)
    summary = apply(h[,1:(dim(tempdata.byGSym.2)[2]-1)],2,max)
  )

  rownames(data.byGSym)=data.byGSym$Symbol
  colnames(data.byGSym)

  data.byGSym.2<-data.byGSym[,-1]
  temp = apply(data.byGSym.2,2,as.numeric)
  rownames(temp) = rownames(data.byGSym.2)
  colnames(temp) = colnames(ndata)

  data.byGSym.2 = temp
  colnames(data.byGSym.2)

  Draw_PCA(data.byGSym.2,output_file,output_dir,data.set123.raw.with.set.label.and.type.cancer[,2])

}
