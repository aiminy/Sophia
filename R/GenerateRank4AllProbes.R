#head(RE41)

#' Title
#'
#' @param TopTableSt41.all.probes.frma
#' @param up_down_top
#'
#' @return
#' @export
#'
#' @examples
GenerateRank4AllProbes<-function(TopTableSt41.all.probes.frma,up_down_top){

  FC.sign=sign(TopTableSt41.all.probes.frma[,1])
  p.inve=1/TopTableSt41.all.probes.frma[,4]
  Gene.Rank.Score=FC.sign*p.inve

  TopTableSt41.gene.frma.GSEA<-cbind(TopTableSt41.all.probes.frma,Gene.Rank.Score)

  TopTableSt41.gene.frma.GSEA.sorted.by.score<-TopTableSt41.gene.frma.GSEA[order(TopTableSt41.gene.frma.GSEA[,7],decreasing = TRUE),]

  Re.GSEA<-cbind(rownames(TopTableSt41.gene.frma.GSEA.sorted.by.score),TopTableSt41.gene.frma.GSEA.sorted.by.score[,7])
  colnames(Re.GSEA)=c("Probe","Score")

  Re.GSEA.2<-rbind(head(Re.GSEA,up_down_top),tail(Re.GSEA,up_down_top))

  #write.table(Re.GSEA.2, file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = FALSE)

  return(Re.GSEA.2)
}
