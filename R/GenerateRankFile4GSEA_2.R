#' Title
#'
#' @param TopTableSt41.gene.frma
#' @param out_dir
#' @param out_file_name
#' @param up_down_top
#'
#' @return
#' @export
#'
#' @examples
GenerateRankFile4GSEA<-function(TopTableSt41.gene.frma,out_dir,out_file_name,up_down_top){

  FC.sign=sign(TopTableSt41.gene.frma[,1])
  p.inve=1/TopTableSt41.gene.frma[,4]
  Gene.Rank.Score=FC.sign*p.inve

  TopTableSt41.gene.frma.GSEA<-cbind(TopTableSt41.gene.frma,Gene.Rank.Score)

  TopTableSt41.gene.frma.GSEA.sorted.by.score<-TopTableSt41.gene.frma.GSEA[order(TopTableSt41.gene.frma.GSEA[,7],decreasing = TRUE),]

  Re.GSEA<-cbind(rownames(TopTableSt41.gene.frma.GSEA.sorted.by.score),TopTableSt41.gene.frma.GSEA.sorted.by.score[,7])
  colnames(Re.GSEA)=c("GeneName","Score")

  Re.GSEA.2<-rbind(head(Re.GSEA,up_down_top),tail(Re.GSEA,up_down_top))

  write.table(Re.GSEA.2, file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = FALSE)

  return(TopTableSt41.gene.frma.GSEA.sorted.by.score)
}
