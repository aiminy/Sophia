OutPut2HtmlTableAllProbes<-function(TopTableSt41.all.probes.frma,out_dir,out_file_name,out_title){


  GeneSym.all <- as.data.frame(getSYMBOL(rownames(TopTableSt41.all.probes.frma), "hgu133plus2.db"))

  #probes.GeneSym.all<-cbind(rownames(TopTableSt41.all.probes.frma),GeneSym.all)

  #rownames(probes.GeneSym.all)
  #hgnc.gene.symbol.ENTREZID<-select(Homo.sapiens, keys=rownames(TopTableSt41.gene.frma),columns=c("SYMBOL","ENTREZID"), keytype="SYMBOL")


  TopTableSt41.all.probes.frma.2<-merge(GeneSym.all,TopTableSt41.all.probes.frma,by=0,sort=FALSE)
  rownames(TopTableSt41.all.probes.frma.2)=TopTableSt41.all.probes.frma.2[,1]

  colnames(TopTableSt41.all.probes.frma.2)[1]="Probe_ID"
  colnames(TopTableSt41.all.probes.frma.2)[2]="SYMBOL"

  TopTableSt41.all.probes.frma.3<-data.frame(TopTableSt41.all.probes.frma.2,Index=seq(1,dim(TopTableSt41.all.probes.frma.2)[1]))

  #colnames(TopTableSt41.gene.frma.2)[1]="SYMBOL"
  ll <- getEG(rownames(TopTableSt41.all.probes.frma.3),"hgu133plus2.db")

  #TopTableSt41.gene.frma.4<-merge(TopTableSt41.gene.frma.2,hgnc.gene.symbol.ENTREZID,by="SYMBOL",sort=FALSE)


  htmlpage(list(ll),filename=paste0(out_dir,out_file_name),
           title=out_title,
           othernames=TopTableSt41.all.probes.frma.3,
           table.head=c("Locus_ID",colnames(TopTableSt41.all.probes.frma.3)),
           table.center=TRUE, digits=6)

  return(TopTableSt41.all.probes.frma.3)
}
