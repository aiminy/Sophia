#Classify cancer patients to 4 subtypes

#' Title
#'
#' @param cel.file.cancer.set1.set2.set3
#' @param data.ER.PR.sample.info
#' @param cut_off
#'
#' @return
#' @export
#'
#' @examples
ClassifyCancerSamplesIntoSubType<-function(cel.file.cancer.set1.set2.set3,data.ER.PR.sample.info,cut_off){
  subtype.1<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,cut_off]==1),1])
  subtype.2<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,cut_off]==2),1])
  subtype.3<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,cut_off]==3),1])
  subtype.4<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,cut_off]==4),1])

  cancer.data.set123.st1<-cel.file.cancer.set1.set2.set3[which(cel.file.cancer.set1.set2.set3[,4] %in% subtype.1),]
  cancer.data.set123.st2<-cel.file.cancer.set1.set2.set3[which(cel.file.cancer.set1.set2.set3[,4] %in% subtype.2),]
  cancer.data.set123.st3<-cel.file.cancer.set1.set2.set3[which(cel.file.cancer.set1.set2.set3[,4] %in% subtype.3),]
  cancer.data.set123.st4<-cel.file.cancer.set1.set2.set3[which(cel.file.cancer.set1.set2.set3[,4] %in% subtype.4),]

  cancer.data.set123.st1.subtype<-cbind(cancer.data.set123.st1,rep("st1",dim(cancer.data.set123.st1)[1]))
  colnames(cancer.data.set123.st1.subtype)[5]="subtype"
  cancer.data.set123.st2.subtype<-cbind(cancer.data.set123.st2,rep("st2",dim(cancer.data.set123.st2)[1]))
  colnames(cancer.data.set123.st2.subtype)[5]="subtype"
  cancer.data.set123.st3.subtype<-cbind(cancer.data.set123.st3,rep("st3",dim(cancer.data.set123.st3)[1]))
  colnames(cancer.data.set123.st3.subtype)[5]="subtype"
  cancer.data.set123.st4.subtype<-cbind(cancer.data.set123.st4,rep("st4",dim(cancer.data.set123.st4)[1]))
  colnames(cancer.data.set123.st4.subtype)[5]="subtype"

  #Use subtype 1,2,3,4
  cel.file.name.key.set.symbol.subtype<-rbind(cancer.data.set123.st1.subtype,cancer.data.set123.st2.subtype,cancer.data.set123.st3.subtype,cancer.data.set123.st4.subtype)

  return(cel.file.name.key.set.symbol.subtype)

}
