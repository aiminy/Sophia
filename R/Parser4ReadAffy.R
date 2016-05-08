#' Title
#'
#' @return
#' @export
#'
#' @examples
Parser4ReadAffy <- function() {
  # Read the CEL files in set1,2 and set3(cancer)
  setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM")
  data.set1 <- ReadAffy()
  ed.raw.set1 <- exprs(data.set1)
  samp.set1 <- sampleNames(data.set1)
  probes.set1 <- featureNames(data.set1)

  setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set2_corrected")
  data.set2 <- ReadAffy()
  ed.raw.set2 <- exprs(data.set2)
  samp.set2 <- sampleNames(data.set2)
  probes.set2 <- featureNames(data.set2)

  setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set3_HGSCancer_only")
  data.set3 <- ReadAffy()
  ed.raw.set3 <- exprs(data.set3)
  samp.set3 <- sampleNames(data.set3)
  probes.set3 <- featureNames(data.set3)

  ll(dim=T)

  data.set1.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/Microarray_Set1_FTE_HGSC_LCM/CASE_DESCRIPTION_SET1.xlsx")
  data.set2.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Set2_header_missing/GSE28044_George.xlsx")
  data.ER.PR.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx")

  original.cel.file.all<-rbind(cbind(unique(samp.set1),rep(1,length(unique(samp.set1)))),
                               cbind(unique(samp.set2),rep(2,length(unique(samp.set2)))),
                               cbind(unique(samp.set3),rep(3,length(unique(samp.set3)))))

  original.cel.file.names.1<-gsub(" - ","_",original.cel.file.all[,1])
  original.cel.file.names.2<-gsub(".CEL","",original.cel.file.names.1)
  original.cel.file.names.3<-gsub(".cel","",original.cel.file.names.2)
  original.cel.file.names.key<-cbind(original.cel.file.all,gsub("-","_",original.cel.file.names.3))
  colnames(original.cel.file.names.key)=c("cel_file_name","set","key")

  #Extract all cel files belonging to the cancer samples in dat set2
  data.set2.cancer.GSM.cel.mapping<-data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)]
  data.set2.cancer.GSM.cel.mapping.2<-data.set2.cancer.GSM.cel.mapping[,c(1,8)]
  data.set2.cancer.GSM.part<-as.character(data.set2.cancer.GSM.cel.mapping.2[which(as.character(data.set2.cancer.GSM.cel.mapping.2[,1])!=""),1])
  data.set2.cancer.GSM.part.cel.name<-as.character(data.set2.cancer.GSM.cel.mapping.2[which(as.character(data.set2.cancer.GSM.cel.mapping.2[,1])!=""),2])



  data.set2.cancer.GSM.cel.mapping.combine<-cbind(data.set2.cancer.GSM.part,ReformatCelName(data.set2.cancer.GSM.part.cel.name))
  GSM.based.cel.index<-which(original.cel.file.names.key[,3] %in% data.set2.cancer.GSM.cel.mapping.combine[,1])
  cel.key.based.index<-which(original.cel.file.names.key[,3] %in% data.set2.cancer.GSM.cel.mapping.combine[,3])
  com.GSM.cel.index<-intersect(GSM.based.cel.index,cel.key.based.index)
  diff.GSM.cel.index<-setdiff(GSM.based.cel.index,cel.key.based.index)
  diff.cel.GSM.index<-setdiff(cel.key.based.index,GSM.based.cel.index)
  original.cel.file.names.key.rm.dup<-original.cel.file.names.key[-diff.cel.GSM.index,]

  original.cel.file.names.set2.cancer.index<-which(original.cel.file.names.key.rm.dup[,3] %in% data.set2.cancer.GSM.cel.mapping.combine[,1])
  original.cel.file.names.set2.cancer<-original.cel.file.names.key.rm.dup[original.cel.file.names.set2.cancer.index,]
  data.set2.cancer.cel.part.name<-as.character(data.set2.cancer.GSM.cel.mapping.2[which(as.character(data.set2.cancer.GSM.cel.mapping.2[,1])==""),2])
  data.set2.cancer.cel.part.name.combine<-cbind(data.set2.cancer.cel.part.name,ReformatCelName(data.set2.cancer.cel.part.name))
  original.cel.file.names.set2.cancer2.index<-which(original.cel.file.names.key.rm.dup[,3] %in% data.set2.cancer.cel.part.name.combine[,3])
  original.cel.file.names.set2.cancer2<-original.cel.file.names.key.rm.dup[original.cel.file.names.set2.cancer2.index,]
  original.cel.file.names.set2.cancer.all<-rbind(original.cel.file.names.set2.cancer,original.cel.file.names.set2.cancer2)
  colnames(original.cel.file.names.set2.cancer.all)=c("cel_file_name","set","key")
  data.set2.cancer.GSM.cel.key<-rbind(data.set2.cancer.GSM.cel.mapping.combine[,c(1,4)],data.set2.cancer.cel.part.name.combine[,c(3,4)])
  colnames(data.set2.cancer.GSM.cel.key)=c("key","symbol")
  original.cel.file.names.set2.cancer.all.2<-merge(original.cel.file.names.set2.cancer.all,data.set2.cancer.GSM.cel.key,by="key",sort=FALSE)

  #Extract all cel files belonging to the cancer samples in dat set1
  cel.file.cancer.data.set1<-as.character(data.set1.sample.info[which(data.set1.sample.info[,4]=="cancer"),2])
  data.set1.cancer.celname.combine<-cbind(cel.file.cancer.data.set1,ReformatCelName(cel.file.cancer.data.set1))
  data.set1.cancer.celname.combine.2<-data.set1.cancer.celname.combine[,3:4]
  colnames(data.set1.cancer.celname.combine.2)=c("key","symbol")
  original.cel.file.names.set1.cancer<-merge(original.cel.file.names.key,data.set1.cancer.celname.combine.2,by="key",sort=FALSE)

  #Generate the symbol for all cel files belonging to the cancer samples in dat set3
  original.cel.file.names.key.data.set3<-original.cel.file.names.key[which(original.cel.file.names.key[,2]==3),]
  original.cel.file.names.key.data.set3.2<-cbind(original.cel.file.names.key.data.set3,sapply(strsplit(original.cel.file.names.key.data.set3[,3],"_"),"[[",4))
  colnames(original.cel.file.names.key.data.set3.2)[4]="symbol"
  original.cel.file.names.key.data.set3.3<-original.cel.file.names.key.data.set3.2[,c(3,1,2,4)]

  #Combine the cancer samples in date set1,2,3
  cel.file.cancer.set1.set2.set3<-rbind(original.cel.file.names.set1.cancer,original.cel.file.names.set2.cancer.all.2,original.cel.file.names.key.data.set3.3)

  #Classify cancer patients to 4 subtypes

  subtype.1<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==1),1])
  subtype.2<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==2),1])
  subtype.3<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==3),1])
  subtype.4<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==4),1])

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

  #Use subtype 1,3,4
  cel.file.name.key.set.symbol.st134<-cel.file.name.key.set.symbol.subtype[which(cel.file.name.key.set.symbol.subtype[,5]!="st2"),]
}
