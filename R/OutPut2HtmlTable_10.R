OutPut2HtmlTable<-function(TopTableSt41.gene.frma,out_dir,out_file_name,out_title){

  hgnc.gene.symbol.ENTREZID<-select(Homo.sapiens, keys=rownames(TopTableSt41.gene.frma),columns=c("SYMBOL","ENTREZID"), keytype="SYMBOL")

  TopTableSt41.gene.frma.2<-data.frame(rownames(TopTableSt41.gene.frma),TopTableSt41.gene.frma)
  colnames(TopTableSt41.gene.frma.2)[1]="SYMBOL"

  TopTableSt41.gene.frma.4<-merge(TopTableSt41.gene.frma.2,hgnc.gene.symbol.ENTREZID,by="SYMBOL",sort=FALSE)

  htmlpage(list(TopTableSt41.gene.frma.4$ENTREZID),filename=paste0(out_dir,out_file_name),
           title=out_title,
           othernames=TopTableSt41.gene.frma.4[,-8],
           table.head=c("ENTREZID",colnames(TopTableSt41.gene.frma.3[,-8])),
           table.center=TRUE, digits=6)
}
