ReformatCelName<-function(data.set2.cancer.GSM.cel.mapping.3){
  data.set2.cancer.GSM.cel.mapping.3.2<-gsub(".cel","",data.set2.cancer.GSM.cel.mapping.3)
  data.set2.cancer.GSM.cel.mapping.3.2.2<-gsub(".CEL","",as.character(data.set2.cancer.GSM.cel.mapping.3.2))
  data.set2.cancer.GSM.cel.mapping.3.2.2.2<-gsub("-","_",as.character(data.set2.cancer.GSM.cel.mapping.3.2.2))
  data.set2.cancer.GSM.cel.mapping.4<-cbind(data.set2.cancer.GSM.cel.mapping.3,data.set2.cancer.GSM.cel.mapping.3.2.2.2,sapply(strsplit(data.set2.cancer.GSM.cel.mapping.3.2.2.2,"_"),"[[",1))
  return(data.set2.cancer.GSM.cel.mapping.4)
}
