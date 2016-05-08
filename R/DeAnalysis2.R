
#' perform differential gene expression analysis
#'
#' @param Re.cutoff.50.old
#'
#' @return
#' @export
#'
#' @examples
DeAnalysis2<-function(Re.cutoff.50.old){

  cel.file.sample.infor.no.8=Re.cutoff.50.old[[1]]
  ed.set123.frma.cancer=Re.cutoff.50.old[[2]]

  f.st134.frma <- factor(cel.file.sample.infor.no.8$subtype)
  design.st134.frma <- model.matrix(~0+f.st134.frma)
  colnames(design.st134.frma) <- levels(f.st134.frma)

  reorder.index.all.probes<-match(cel.file.sample.infor.no.8[,6],colnames(ed.set123.frma.cancer))
  ed.set123.frma.cancer.reorder<-ed.set123.frma.cancer[,reorder.index.all.probes]

  check.match<-cbind(colnames(ed.set123.frma.cancer.reorder),design.st134.frma,cel.file.sample.infor.no.8)

  fit.st134.all.probes.frma <- lmFit(ed.set123.frma.cancer.reorder, design.st134.frma)
  #cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134.frma)

  #print(f.st134.frma)
  #cat(length(levels(f.st134.frma)))

  #print(combn(levels(f.st134.frma),2))

  Contrast.com<-combn(levels(f.st134.frma),2)

  print(dim(Contrast.com))

  x=paste0(Contrast.com[1,],"-",Contrast.com[2,])

  cont.matrix.st134.frma <- makeContrasts(contrasts=x,levels=design.st134.frma)

  fit2.st134.all.probes.frma  <- contrasts.fit(fit.st134.all.probes.frma, cont.matrix.st134.frma)
  fit2.st134.all.probes.frma  <- eBayes(fit2.st134.all.probes.frma)

  re<-list()

  re[[1]]<-topTable(fit2.st134.all.probes.frma,coef=1,n=dim(ed.set123.frma.cancer.reorder)[1])
  names(re)[1]=colnames(cont.matrix.st134.frma)[1]

  num.contrast=dim(cont.matrix.st134.frma)[2]
  for(i in 2:num.contrast){

    re[[i]]<-topTable(fit2.st134.all.probes.frma,coef=i,n=dim(ed.set123.frma.cancer.reorder)[1])
    names(re)[i]=colnames(cont.matrix.st134.frma)[i]
  }

  re[[num.contrast+1]]=cont.matrix.st134.frma
  names(re)[num.contrast+1]="ct"

  #TopTableSt34.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=1,n=dim(ed.set123.frma.cancer.reorder)[1])
  #TopTableSt41.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=2,n=dim(ed.set123.frma.cancer.reorder)[1])
  #TopTableSt13.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=3,n=dim(ed.set123.frma.cancer.reorder)[1])
  # TopTableSt12.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=4,n=dim(ed.set123.frma.cancer.reorder)[1])
  # TopTableSt23.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=5,n=dim(ed.set123.frma.cancer.reorder)[1])
  # TopTableSt24.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=6,n=dim(ed.set123.frma.cancer.reorder)[1])
  #
  #
  #re<-list(ct1=TopTableSt34.all.probes.frma,ct2=TopTableSt41.all.probes.frma,ct3=TopTableSt13.all.probes.frma,ct=cont.matrix.st134.frma)

  return(re)
}







