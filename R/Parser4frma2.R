#Use frma method to do normalization
#' Title
#'
#' @return
#' @export
#'
#' @examples
#' Parser4frma2("/media/H_driver/2015/Sophia/Cel_file_frma/","heatmap_allsample_frma.pdf","PCA_allsample_frma.pdf","/media/H_driver/2015/Sophia/")
#'
Parser4frma2<- function(dir_cel_file_frma,cel.file.sample.infor.no.2,heatmap_output_file,PCA_output_file,output_dir_re) {
  setwd(dir_cel_file_frma)
  data.set123.raw <- ReadAffy()
  data.set123.frma<-frma(data.set123.raw)
  ed.set123.frma <- exprs(data.set123.frma)
  samp.frma <- sampleNames(data.set123.frma)
  probes.frma <- featureNames(data.set123.frma)

  cel.file.sample.infor.no.3<-gsub("X","",cel.file.sample.infor.no.2[,1])
  cel.file.sample.infor.no.4<-gsub(".CEL","",cel.file.sample.infor.no.3)
  cel.file.sample.infor.no.5<-gsub(".cel","",cel.file.sample.infor.no.4)
  cel.file.sample.infor.no.6<-gsub("\\.","_",cel.file.sample.infor.no.5)
  cel.file.sample.infor.no.7<-cbind(cel.file.sample.infor.no.2,cel.file.sample.infor.no.6)
  colnames(cel.file.sample.infor.no.7)=c("filename","subtype","match_key")

  samp.frma.1<-gsub(" - ","_",samp.frma)
  samp.frma.2<-gsub(".CEL","",samp.frma.1)
  samp.frma.3<-gsub(".cel","",samp.frma.2)
  samp.frma.4<-gsub("-","_",samp.frma.3)
  samp.frma.5<-cbind(samp.frma,samp.frma.4)
  samp.frma.cancer.index<-which(samp.frma.5[,2] %in% cel.file.sample.infor.no.7[,3])
  samp.frma.6<-samp.frma.5[samp.frma.cancer.index,]

  setdiff(cel.file.sample.infor.no.3[,3],samp.frma.2[,2])
  samp.frma.1[grep(1443,samp.frma.1[,2]),]

  data.set123.frma.cancer<-data.set123.frma[,samp.frma.cancer.index]
  ed.set123.frma.cancer<-ed.set123.frma[,samp.frma.cancer.index]

  dim(ed.set123.frma.cancer)
  samp.frma.7<-cbind(samp.frma.6,colnames(ed.set123.frma.cancer))
  samp.frma.8<-samp.frma.7[,2:3]
  colnames(samp.frma.8)<-c("match_key","filename_frma")
  cel.file.sample.infor.no.8<-merge(cel.file.sample.infor.no.7,samp.frma.8,by="match_key")

  head(samp.frma.7)
  dim(ed.set123.frma.cancer)

  GeneSym.all <- getSYMBOL(rownames(ed.set123.frma.cancer), "hgu133plus2.db")
  ndata<-ed.set123.frma.cancer
  geneSymbol=GeneSym.all
  tempdata.byGSym = data.frame(ndata, Symbol = geneSymbol)

  tempdata.byGSym.2<-tempdata.byGSym[-which(is.na(tempdata.byGSym$Symbol)),]

  rownames(tempdata.byGSym.2) = NULL
  data.byGSym = ddply(tempdata.byGSym.2, c("Symbol"),function(h)
    summary = apply(h,2,max)
  )

  data.byGSym.index = ddply(tempdata.byGSym.2, c("Symbol"),function(h)
    summary = apply(h,2,which.max)
  )
  #class(data.byGSym)
  #which(is.na(data.byGSym$Symbol))
  rownames(data.byGSym)=data.byGSym$Symbol
  data.byGSym.2<-data.byGSym[,-61]
  temp = apply(data.byGSym.2,2,as.numeric)
  rownames(temp) = rownames(data.byGSym.2)
  data.byGSym.2 = temp
  SampleType<-cel.file.sample.infor.no.8$subtype
  heatmap_wPCA(data.byGSym.2,heatmap_output_file,PCA_output_file,output_dir_re,SampleType)
}
