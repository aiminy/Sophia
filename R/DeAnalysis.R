##Use all probes for DE
#' Title
#'
#' @param Re.cutoff.50.old
#'
#' @return
#' @export
#'
#' @examples
DeAnalysis<-function(Re.cutoff.50.old){

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
  cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",st13="st1-st3",st12="st1-st2",st23="st2-st3",st24="st2-st4",levels=design.st134.frma)

  fit2.st134.all.probes.frma  <- contrasts.fit(fit.st134.all.probes.frma, cont.matrix.st134.frma)
  fit2.st134.all.probes.frma  <- eBayes(fit2.st134.all.probes.frma)

  TopTableSt34.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=1,n=dim(ed.set123.frma.cancer.reorder)[1])
  TopTableSt41.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=2,n=dim(ed.set123.frma.cancer.reorder)[1])
  TopTableSt13.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=3,n=dim(ed.set123.frma.cancer.reorder)[1])
  TopTableSt12.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=4,n=dim(ed.set123.frma.cancer.reorder)[1])
  TopTableSt23.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=5,n=dim(ed.set123.frma.cancer.reorder)[1])
  TopTableSt24.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=6,n=dim(ed.set123.frma.cancer.reorder)[1])


  re<-list(ct1=TopTableSt34.all.probes.frma,ct2=TopTableSt41.all.probes.frma,ct3=TopTableSt13.all.probes.frma,
           ct4=TopTableSt12.all.probes.frma,ct5=TopTableSt23.all.probes.frma,ct6=TopTableSt24.all.probes.frma,
           ct=cont.matrix.st134.frma)

  return(re)

}
