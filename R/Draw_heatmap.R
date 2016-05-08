#' draw heatmap only
#'
#' @param Data
#' @param output_heatmap
#' @param out_dir
#' @param g_level
#'
#' @return
#' @export
#'
#' @examples
Draw_heatmap = function(Data,output_heatmap, out_dir, g_level = NULL){
  hmcol<-rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
  if(is.null(g_level)){
    type_level = 1:ncol(Data)
    col_level = "black"
  }else{
    type_level = 16
    TEMP = factor(g_level)
    uniq_label =  levels(TEMP)
    levels(TEMP) = hmcol[ceiling(seq(length.out=length(levels(TEMP)),from=1,to=256))]
    col_level = as.character(TEMP)
    uniq_col = levels(TEMP)
  }

  Data.hc = hclust(dist(Data), method="average")
  rowInd <- order.dendrogram(as.dendrogram(Data.hc))
  pdf(file=paste0(out_dir,"/",output_heatmap), width=10, height=12)
  heatmap.2(Data[rowInd,], col=hmcol, Rowv = F, dendrogram = "column", scale="row",labRow = NA,key=TRUE, keysize=0.55, symkey=FALSE,density.info="none", trace="none",cexCol=1,margins=c(8,12))
  dev.off()
}
