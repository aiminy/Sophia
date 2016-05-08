MapingCelFile4ReadFram<-function(samp.frma,data.set123,ed.set123.frma,cel.file.name.key.set.symbol.st134){

  samp.frma.1<-gsub(" - ","_",samp.frma)
  samp.frma.2<-gsub(".CEL","",samp.frma.1)
  samp.frma.3<-gsub(".cel","",samp.frma.2)
  samp.frma.4<-gsub("-","_",samp.frma.3)
  samp.frma.5<-cbind(samp.frma,samp.frma.4)

  samp.frma.cancer.index<-which(samp.frma.5[,2] %in% cel.file.name.key.set.symbol.st134[,1])
  samp.frma.6<-samp.frma.5[samp.frma.cancer.index,]

  data.set123.frma.cancer<-data.set123.frma[,samp.frma.cancer.index]
  ed.set123.frma.cancer<-ed.set123.frma[,samp.frma.cancer.index]
  #dim(ed.set123.frma.cancer)

  samp.frma.7<-cbind(samp.frma.6,colnames(ed.set123.frma.cancer))
  samp.frma.8<-samp.frma.7[,2:3]
  colnames(samp.frma.8)<-c("key","filename_frma")

  cel.file.sample.infor.no.8<-merge(cel.file.name.key.set.symbol.st134,samp.frma.8,by="key",sort=FALSE)

  Re<-list(fileName=cel.file.sample.infor.no.8,fileData=ed.set123.frma.cancer)

  return(Re)
}
