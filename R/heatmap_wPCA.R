#' Draw heamap and PCA based on gene expression data
#'
#' @param Data
#' @param output_heatmap
#' @param output_pca
#' @param out_dir
#' @param g_level
#'
#' @return
#' @export
#'
#' @examples
heatmap_wPCA = function(Data, output_heatmap, output_pca, out_dir, g_level = NULL){
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

  pdf(file=paste0(out_dir,"/",output_pca), width=12, height=12)
  Data.pca = prcomp(t(Data))
  with(data.frame(Data.pca$x), scatter3D(PC1, PC2, PC3, colvar = NULL, type="h",
                                         ticktype = "detailed", bty="b2", cex=1,
                                         xlab="PC 1",	ylab="PC 2",zlab="PC 3", theta = 40, phi = 40, pch=type_level,
                                         col=col_level,
                                         main = "Principal component analysis")
  )
  legend("topright", legend = uniq_label, pch=type_level,
         col = uniq_col,
         cex=1, inset=c(0.02))
  with(data.frame(Data.pca$x), text3D(x=PC1, y=PC2,
                                      z=PC3, colnames(Data), col = "black", add=TRUE, colkey = FALSE, cex=0.5)
  )
  dev.off()
}
