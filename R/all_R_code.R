biocLite("SpikeInSubset")
biocLite("hgu95av2cdf")
library(SpikeInSubset)
library(limma)
data(spikein95)
rma.eset <- rma(spikein95)
rma.e <- exprs(rma.eset)

head(rma.e)

design <- model.matrix(~factor(rma.eset$population))
head(rma.eset)

biocLite("hgu95av2.db")
library(limma)
targets <- readTargets("targets.txt")
data(RG)
RG <- backgroundCorrect(RG, method="normexp")
MA <- normalizeWithinArrays(RG)
targets <- data.frame(Cy3=I(rep("Pool",6)),Cy5=I(c("WT","WT","WT","KO","KO","KO")))
design <- modelMatrix(targets, ref="Pool")
arrayw <- arrayWeightsSimple(MA, design)
fit <- lmFit(MA, design, weights=arrayw)
fit2 <- contrasts.fit(fit, contrasts=c(-1,1))
fit2 <- eBayes(fit2)

setwd("/media/H_driver/Aimin_project/estrogen")
dir()
targets <- readTargets("estrogen.txt", sep="")
targets
library(affy)
library(hgu95av2cdf)
abatch <- ReadAffy(filenames=targets$filename)
eset <- rma(abatch)
expr(eset)
f <- paste(targets$estrogen,targets$time.h,sep="")

f <- factor(f)
f
design <- model.matrix(~0+f)
colnames(design) <- levels(f)
design
fit <- lmFit(eset, design)
names(fit)

cont.matrix <- makeContrasts(E10="present10-absent10",E48="present48-absent48",Time="absent48-absent10",levels=design)
cont.matrix

fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)

colnames(fit2)
topTable(fit2,coef=1)
topTable(fit2,coef=1,adjust="fdr")
topTable(fit2,coef=2)
topTable(fit2,coef=2,adjust="fdr")

results <- decideTests(fit2)
summary(results)
vennDiagram(results)

library(hgu95av2cdf)
library(hgu95av2.db)

geneIDs <- ls(hgu95av2cdf)

geneSymbols <- as.character(unlist(lapply(mget(geneIDs,env=hgu95av2SYMBOL),
                                          function (symbol) { return(paste(symbol,collapse="; ")) } )))
geneNames <- as.character(unlist(lapply(mget(geneIDs,env=hgu95av2GENENAME),
                                        function (name) { return(paste(name,collapse="; ")) } )))
unigene <- as.character(unlist(lapply(mget(geneIDs,env=hgu95av2UNIGENE),
                                      function (unigeneID) { return(paste(unigeneID,collapse="; ")) } )))

geneNames <- substring(geneNames,1,40)
unigene <- gsub("Hs\\.","",unigene)

genelist <- data.frame(GeneID=geneIDs,GeneSymbol=geneSymbols,GeneName=geneNames,
                       UniGeneHsID=paste("<a href=http://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Hs&CID=",
                                         unigene,">",unigene,"</a>",sep=""))

unigeneTopTableEst10 <- topTable(fit2,coef=1,n=20,genelist=genelist)
unigeneTopTableEst48 <- topTable(fit2,coef=2,n=20,genelist=genelist)

library(xtable)
xtableUnigeneEst10 <- xtable(unigeneTopTableEst10,display=c("d","s","s","s","s","g","g","g","e","g","g"))
xtableUnigeneEst48 <- xtable(unigeneTopTableEst48,display=c("d","s","s","s","s","g","g","g","e","g","g"))

cat(file="estrogenUniGeneE10.html","<html>\n<body>")
print.xtable(xtableUnigeneEst10,type="html",file="estrogenUniGeneE10.html",append=TRUE)
cat(file="estrogenUniGeneE10.html","</body>\n</html>",append=TRUE)

cat(file="estrogenUniGeneE48.html","<html>\n<body>")
print.xtable(xtableUnigeneEst48,type="html",file="estrogenUniGeneE48.html",append=TRUE)
cat(file="estrogenUniGeneE48.html","</body>\n</html>",append=TRUE)
#load the affy library
library(gdata)
library(affy)
library(affyio)
library(annotate)

biocLite("mogene20sttranscriptcluster.db")
library(mogene20stprobeset.db)
library(mogene20sttranscriptcluster.db)

biocLite("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)

biocLite("pd.hg.u133.plus.2")
library(pd.hg.u133.plus.2)

# Read in the CEL files in the directory, then normalize the data
setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM")
data.set1 <- ReadAffy()
ed.raw.set1 <- exprs(data.set1)
samp.set1 <- sampleNames(data.set1)
probes.set1 <- featureNames(data.set1)
ed.normalized.set1<- rma(data.set1)

setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set2_corrected")
data.set2 <- ReadAffy()
ed.raw.set2 <- exprs(data.set2)
samp.set2 <- sampleNames(data.set2)
probes.set2 <- featureNames(data.set2)
ed.normalized.set2<- rma(data.set2)

setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set3_HGSCancer_only")
data.set3 <- ReadAffy()
ed.raw.set3 <- exprs(data.set3)
samp.set3 <- sampleNames(data.set3)
probes.set3 <- featureNames(data.set3)
ed.normalized.set3<- rma(data.set3)

data.set1.normalized<-ed.normalized.set1
data.set2.normalized<-ed.normalized.set2
data.set3.normalized<-ed.normalized.set3

ll(dim=T)

colnames(data.set1.normalized)
colnames(data.set2.normalized)
colnames(data.set3.normalized)

length(colnames(data.set1.normalized))
length(colnames(data.set2.normalized))
length(colnames(data.set3.normalized))

unique(colnames(data.set1.normalized))
length(unique(colnames(data.set1.normalized)))
length(unique(colnames(data.set2.normalized)))
length(unique(colnames(data.set3.normalized)))

data.set1.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/Microarray_Set1_FTE_HGSC_LCM/CASE_DESCRIPTION_SET1.xlsx")
data.set2.sample.info<-read.xls("/home/aiminyan/Downloads/GSE28044_George.xlsx")
data.ER.PR.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx")
head(data.ER.PR.sample.info)

head(data.set1.sample.info)
data.set1.sample.info[,2]
data.set1.normalize
data.set2.normalized

grep(2367,unique(data.set1.sample.info[,2]))

unique(data.set2.sample.info[,2])
head(data.ER.PR.sample.info)

dim(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==1),])
dim(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==2),])
dim(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==3),])
dim(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==4),])



#ER-PR-
ER.PR.sample1.from.set1<-MapSample2CELfile(data.ER.PR.sample.info,data.set1.normalized,1)
ER.PR.sample1.from.set2<-MapSample2CELfile(data.ER.PR.sample.info,data.set2.normalized,1)
ER.PR.sample1.from.set3<-MapSample2CELfile(data.ER.PR.sample.info,data.set3.normalized,1)

#ER-PR+
ER.PR.sample2.from.set1<-MapSample2CELfile(data.ER.PR.sample.info,data.set1.normalized,2)
ER.PR.sample2.from.set2<-MapSample2CELfile(data.ER.PR.sample.info,data.set2.normalized,2)
ER.PR.sample2.from.set3<-MapSample2CELfile(data.ER.PR.sample.info,data.set3.normalized,2)

#ER+PR-
ER.PR.sample3.from.set1<-MapSample2CELfile(data.ER.PR.sample.info,data.set1.normalized,3)
ER.PR.sample3.from.set2<-MapSample2CELfile(data.ER.PR.sample.info,data.set2.normalized,3)
ER.PR.sample3.from.set3<-MapSample2CELfile(data.ER.PR.sample.info,data.set3.normalized,3)

#ER+PR+
ER.PR.sample4.from.set1<-MapSample2CELfile(data.ER.PR.sample.info,data.set1.normalized,4)
ER.PR.sample4.from.set2<-MapSample2CELfile(data.ER.PR.sample.info,data.set2.normalized,4)
ER.PR.sample4.from.set3<-MapSample2CELfile(data.ER.PR.sample.info,data.set3.normalized,4)

head(ER.PR.sample1.from.set1)
head(ER.PR.sample1.from.set2)
head(ER.PR.sample1.from.set3)

head(ER.PR.sample3.from.set1)
head(ER.PR.sample3.from.set2)
head(ER.PR.sample3.from.set3)

head(data.set1.sample.info)
head(data.set2.sample.info)
head(data.ER.PR.sample.info)

data.set1.sample.info[,2]

data.set2.sample.info[,8:9:10]

cel.file.all<-c(unique(colnames(data.set1.normalized)),unique(colnames(data.set2.normalized)),unique(colnames(data.set3.normalized)))

length(cel.file.all)
cel.file.all

cel.file.all.2<-rbind(cbind(unique(colnames(data.set1.normalized)),rep(1,length(unique(colnames(data.set1.normalized))))),
cbind(unique(colnames(data.set2.normalized)),rep(2,length(unique(colnames(data.set2.normalized))))),
cbind(unique(colnames(data.set3.normalized)),rep(3,length(unique(colnames(data.set3.normalized))))))

cel.file.cancer.data.set2<-as.character(data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)][,8])
cel.file.cancer.data.set2.2<-as.character(data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)][,1])
cel.file.cancer.data.set2.2.2<-cel.file.cancer.data.set2.2[-which(cel.file.cancer.data.set2.2=="")]

data.set2.cancer.GSM.cel.mapping<-data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)]
data.set2.cancer.GSM.cel.mapping.2<-data.set2.cancer.GSM.cel.mapping[,c(1,8)]
data.set2.cancer.GSM.cel.mapping.3<-data.set2.cancer.GSM.cel.mapping.2[which(as.character(data.set2.cancer.GSM.cel.mapping.2[,1])!=""),]

data.set2.cancer.GSM.cel.mapping.3.2<-gsub(".cel","",as.character(data.set2.cancer.GSM.cel.mapping.3[,2]))
data.set2.cancer.GSM.cel.mapping.3.2.2<-gsub(".CEL","",as.character(data.set2.cancer.GSM.cel.mapping.3.2))
data.set2.cancer.GSM.cel.mapping.3.2.2.2<-gsub("-","_",as.character(data.set2.cancer.GSM.cel.mapping.3.2.2))
data.set2.cancer.GSM.cel.mapping.4<-cbind(data.set2.cancer.GSM.cel.mapping.3,sapply(strsplit(data.set2.cancer.GSM.cel.mapping.3.2.2.2,"_"),"[[",1))

cel.file.cancer.data.set1<-as.character(data.set1.sample.info[which(data.set1.sample.info[,4]=="cancer"),2])

cel.file.cancer.set1.set2<-rbind(cbind(cel.file.cancer.data.set1,rep(1,length(cel.file.cancer.data.set1))),
cbind(cel.file.cancer.data.set2,rep(2,length(cel.file.cancer.data.set2))),
cbind(cel.file.cancer.data.set2.2.2,rep(2,length(cel.file.cancer.data.set2.2.2))))

cel.file.cancer.set1.set2.set3<-rbind(cel.file.cancer.set1.set2,cel.file.all.2[which(cel.file.all.2[,2]==3),])

cel.file.cancer.set1.set2.set3.name<-gsub(".CEL","",cel.file.cancer.set1.set2.set3[,1])
cel.file.cancer.set1.set2.set3.name.1<-gsub(".cel","",cel.file.cancer.set1.set2.set3.name)
cel.file.cancer.set1.set2.set3.name.2<-gsub("-","_",cel.file.cancer.set1.set2.set3.name.1)

data.set123.normalized<-cbind(data.set1.normalized,data.set2.normalized,data.set3.normalized)
dim(data.set123.normalized)
dim(data.set123.normalized[,which(colnames(data.set123.normalized) %in% cel.file.cancer.set1.set2.set3[,1])])
colnames(data.set123.normalized)[grep(54,colnames(data.set123.normalized))]

original.cel.file.names<-colnames(data.set123.normalized)
original.cel.file.names.1<-gsub("X","",original.cel.file.names)
original.cel.file.names.2<-gsub(".CEL","",original.cel.file.names.1)
original.cel.file.names.3<-gsub(".cel","",original.cel.file.names.2)
original.cel.file.names.4<-gsub("\\.","_",original.cel.file.names.3)

index.4.cancer.sample<-which(original.cel.file.names.4 %in% cel.file.cancer.set1.set2.set3.name.2)
original.cel.file.name.with.mapped.names<-cbind(original.cel.file.names[index.4.cancer.sample],original.cel.file.names.4[index.4.cancer.sample])
#index.4.no.cancer.sample.1<-cbind(original.cel.file.names[-index.4.cancer.sample],original.cel.file.names.4[-index.4.cancer.sample])
setdiff(cel.file.cancer.set1.set2.set3.name.2,original.cel.file.names.4[index.4.cancer.sample])
which(original.cel.file.names.4[-index.4.cancer.sample] %in% cel.file.cancer.set1.set2.set3.name.2)

cancer.data.set123.normalized<-data.set123.normalized[,index.4.cancer.sample]
dim(cancer.data.set123.normalized)

cancer.cel.file.name.72<-colnames(cancer.data.set123.normalized)
cancer.cel.file.name.72.1<-gsub("X","",cancer.cel.file.name.72)
cancer.cel.file.name.72.2<-gsub(".CEL","",cancer.cel.file.name.72.1)
cancer.cel.file.name.72.3<-gsub(".cel","",cancer.cel.file.name.72.2)
cancer.cel.file.name.72.4<-gsub("\\.","_",cancer.cel.file.name.72.3)
cancer.cel.file.name.72.5<-cancer.cel.file.name.72.4

data.set2.cancer.GSM.cel.mapping.44<-as.character(data.set2.cancer.GSM.cel.mapping.4[which(data.set2.cancer.GSM.cel.mapping.4[,1] %in% cancer.cel.file.name.72.4),3])
cancer.cel.file.name.72.5[which(cancer.cel.file.name.72.5 %in% data.set2.cancer.GSM.cel.mapping.4[,1])]<-data.set2.cancer.GSM.cel.mapping.44
cancer.cel.file.name.72.6<-cbind(cancer.cel.file.name.72,cancer.cel.file.name.72.5)
cancer.cel.file.name.72.7<-c(sapply(strsplit(cancer.cel.file.name.72.6[1:29,2],"_"),"[[",1),
sapply(strsplit(cancer.cel.file.name.72.6[30:72,2],"_"),"[[",4))
cancer.cel.file.name.72.8<-cbind(cancer.cel.file.name.72.6,cancer.cel.file.name.72.7)

#Classify cancer patients to 4 subtypes

subtype.1<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==1),1])
subtype.2<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==2),1])
subtype.3<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==3),1])
subtype.4<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==4),1])

cancer.data.set123.normalized.st1<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.1)]
cancer.data.set123.normalized.st2<-as.data.frame(cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.2)])
colnames(cancer.data.set123.normalized.st2)<-colnames(cancer.data.set123.normalized)[which(cancer.cel.file.name.72.8[,3] %in% subtype.2)]
cancer.data.set123.normalized.st3<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.3)]
cancer.data.set123.normalized.st4<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.4)]

dim(cancer.data.set123.normalized.st1)
dim(cancer.data.set123.normalized.st2)
dim(cancer.data.set123.normalized.st3)
dim(cancer.data.set123.normalized.st4)

head(cancer.data.set123.normalized.st1)
head(cancer.data.set123.normalized.st2)
head(cancer.data.set123.normalized.st3)
head(cancer.data.set123.normalized.st4)

colnames(cancer.data.set123.normalized.st1)
colnames(cancer.data.set123.normalized.st2)
colnames(cancer.data.set123.normalized.st3)
colnames(cancer.data.set123.normalized.st4)

#n.st<-length(colnames(cancer.data.set123.normalized.st1))

# #Use subtype1,2,3,4
# cel.file.sample.infor<-as.data.frame(rbind(
# cbind(colnames(cancer.data.set123.normalized.st1),rep("st1",length(colnames(cancer.data.set123.normalized.st1)))),
# cbind(colnames(cancer.data.set123.normalized.st2),rep("st2",length(colnames(cancer.data.set123.normalized.st2)))),
# cbind(colnames(cancer.data.set123.normalized.st3),rep("st3",length(colnames(cancer.data.set123.normalized.st3)))),
# cbind(colnames(cancer.data.set123.normalized.st4),rep("st4",length(colnames(cancer.data.set123.normalized.st4))))
# ))
# colnames(cel.file.sample.infor)=c("filename","subtype")

#Use subtype 1,3,4
cel.file.sample.infor.no.2<-as.data.frame(rbind(
  cbind(colnames(cancer.data.set123.normalized.st1),rep("st1",length(colnames(cancer.data.set123.normalized.st1)))),
  cbind(colnames(cancer.data.set123.normalized.st3),rep("st3",length(colnames(cancer.data.set123.normalized.st3)))),
  cbind(colnames(cancer.data.set123.normalized.st4),rep("st4",length(colnames(cancer.data.set123.normalized.st4))))
))

colnames(cel.file.sample.infor.no.2)=c("filename","subtype")
f.st134 <- factor(cel.file.sample.infor.no.2$subtype)
design.st134 <- model.matrix(~0+f.st134)
colnames(design.st134) <- levels(f.st134)

cancer.data.st134<-cbind(cancer.data.set123.normalized.st1,
cancer.data.set123.normalized.st3,
cancer.data.set123.normalized.st4)
head(cancer.data.st134)

fit.st134 <- lmFit(cancer.data.st134, design.st134)

cont.matrix.st134 <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134)
fit2.st134  <- contrasts.fit(fit.st134, cont.matrix.st134)
fit2.st134  <- eBayes(fit2.st134)

str(fit2.st134)

TopTableSt34.all<-topTable(fit2.st134,coef=1,n=54674)

TopTableSt34.100<-topTable(fit2.st134,coef=1,adjust="fdr",n=100)

genenames <- as.character(rownames(TopTableSt34.54675))

length(genenames)

genenames

annotation(ed.normalized.set1)
annotation(ed.normalized.set2)
annotation(ed.normalized.set3)
library("hgu133plus2.db")

map <- getpAnnMa("ENTREZID", "hgu133plus2", load=TRUE, type=c("env", "db"))
class(map)
ll <- getEG(genenames,"hgu133plus2.db")
GeneSym <- getSYMBOL(genenames, "hgu133plus2.db")

Probe.gene.sym<-(cbind(genenames,GeneSym))

Probe.gene.sym

length(which(is.na(Probe.gene.sym[,2])))


tab <- data.frame(GeneSym, TopTableSt34.100)
tab <- data.frame(rownames(tab), tab)

colnames(tab)[1] <- c("Probe ID")
ll <- list(ll)
htmlpage(ll, filename="/media/H_driver/2015/Sophia/St34-4.html", title="HTML report",
         othernames=tab, table.head=c("Locus ID",colnames(tab)), table.center=TRUE, digits=6)

colnames(fit2.st134)
topTable(fit2.st134,coef=1)
topTable(fit2.st134,coef=1,adjust="fdr")

topTable(fit2.st134,coef=2)
topTable(fit2.st134,coef=2,adjust="fdr")

results.st134 <- decideTests(fit2.st134,p.value=0.99)

summary(results.st134)
vennDiagram(results.st134)

fit2.st134$genes

unigeneTopTableSt34 <- topTable(fit2.st134,coef=1,n=20,genelist=genelist)
unigeneTopTableSt41 <- topTable(fit2.st134,coef=2,n=20,genelist=genelist)

library(xtable)
xtableUnigeneSt34 <- xtable(unigeneTopTableSt34,display=c("s","s","s","s","g","g","g","e","e","g","g"))
xtableUnigeneSt41 <- xtable(unigeneTopTableSt41,display=c("s","s","s","s","g","g","g","e","e","g","g"))

cat(file="/media/H_driver/2015/Sophia/St34.html","<html>\n<body>")
print.xtable(xtableUnigeneSt34,type="html",file="/media/H_driver/2015/Sophia/St34.html",append=TRUE)
cat(file="/media/H_driver/2015/Sophia/St34.html","</body>\n</html>",append=TRUE)

cat(file="/media/H_driver/2015/Sophia/St41.html","<html>\n<body>")
print.xtable(xtableUnigeneSt41,type="html",file="/media/H_driver/2015/Sophia/St41.html",append=TRUE)
cat(file="/media/H_driver/2015/Sophia/St41.html","</body>\n</html>",append=TRUE)

#fit#Differential gene expression analysis for st3 vs st4
#probe<-rownames(cancer.data.set123.normalized.st1)
#geneSymbol<-getSYMBOL(probe, "hugene10sttranscriptcluster.db")
#Differential gene expression analysis for st2 vs st1
#Differential gene expression analysis for st4 vs st1
save.image(file="/media/H_driver/Aimin_project/Sophia/Data_set_3.RData")
#load the affy library
biocLite("mogene20sttranscriptcluster.db")
library(gdata)
library(affy)
library(affyio)
library(annotate)
library(mogene20stprobeset.db)
library(mogene20sttranscriptcluster.db)

biocLite("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)

# Read in the CEL files in the directory, then normalize the data
setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM")
data.set1 <- ReadAffy()
ed.raw.set1 <- exprs(data.set1)
samp.set1 <- sampleNames(data.set1)
probes.set1 <- featureNames(data.set1)
ed.normalized.set1<- rma(data.set1)

setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set2_corrected")
data.set2 <- ReadAffy()
ed.raw.set2 <- exprs(data.set2)
samp.set2 <- sampleNames(data.set2)
probes.set2 <- featureNames(data.set2)
ed.normalized.set2<- rma(data.set2)

setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set3_HGSCancer_only")
data.set3 <- ReadAffy()
ed.raw.set3 <- exprs(data.set3)
samp.set3 <- sampleNames(data.set3)
probes.set3 <- featureNames(data.set3)
ed.normalized.set3<- rma(data.set3)

data.set1.normalized<-ed.normalized.set1
data.set2.normalized<-ed.normalized.set2
data.set3.normalized<-ed.normalized.set3

ll(dim=T)

colnames(data.set1.normalized)
colnames(data.set2.normalized)
colnames(data.set3.normalized)

length(colnames(data.set1.normalized))
length(colnames(data.set2.normalized))
length(colnames(data.set3.normalized))

unique(colnames(data.set1.normalized))
length(unique(colnames(data.set1.normalized)))
length(unique(colnames(data.set2.normalized)))
length(unique(colnames(data.set3.normalized)))

data.set1.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/Microarray_Set1_FTE_HGSC_LCM/CASE_DESCRIPTION_SET1.xlsx")
data.set2.sample.info<-read.xls("/home/aiminyan/Downloads/GSE28044_George.xlsx")
data.ER.PR.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx")
head(data.ER.PR.sample.info)

head(data.set1.sample.info)
data.set1.sample.info[,2]
data.set1.normalize
data.set2.normalized

grep(2367,unique(data.set1.sample.info[,2]))

unique(data.set2.sample.info[,2])
head(data.ER.PR.sample.info)

dim(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==1),])
dim(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==2),])
dim(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==3),])
dim(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==4),])

MapSample2CELfile<-function(data.ER.PR.sample.info,data.set1.normalized,sample_type){

cat(dim(data.ER.PR.sample.info),dim(data.set1.normalized),sample_type,"\n")

num.sample1=length(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==1),][,1])
cel_file_index_4_sample<-array()

for(i in 1:num.sample1) {
tmp<-grep(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==sample_type),][i,1],colnames(data.set1.normalized))
cel_file_index_4_sample<-c(cel_file_index_4_sample,tmp)
}
cel_file_index_4_sample<-cel_file_index_4_sample[-1]

colnames(data.set1.normalized)[cel_file_index_4_sample]

re<-data.set1.normalized[,cel_file_index_4_sample]

return(re)

}

#ER-PR-
ER.PR.sample1.from.set1<-MapSample2CELfile(data.ER.PR.sample.info,data.set1.normalized,1)
ER.PR.sample1.from.set2<-MapSample2CELfile(data.ER.PR.sample.info,data.set2.normalized,1)
ER.PR.sample1.from.set3<-MapSample2CELfile(data.ER.PR.sample.info,data.set3.normalized,1)

#ER-PR+
ER.PR.sample2.from.set1<-MapSample2CELfile(data.ER.PR.sample.info,data.set1.normalized,2)
ER.PR.sample2.from.set2<-MapSample2CELfile(data.ER.PR.sample.info,data.set2.normalized,2)
ER.PR.sample2.from.set3<-MapSample2CELfile(data.ER.PR.sample.info,data.set3.normalized,2)

#ER+PR-
ER.PR.sample3.from.set1<-MapSample2CELfile(data.ER.PR.sample.info,data.set1.normalized,3)
ER.PR.sample3.from.set2<-MapSample2CELfile(data.ER.PR.sample.info,data.set2.normalized,3)
ER.PR.sample3.from.set3<-MapSample2CELfile(data.ER.PR.sample.info,data.set3.normalized,3)

#ER+PR+
ER.PR.sample4.from.set1<-MapSample2CELfile(data.ER.PR.sample.info,data.set1.normalized,4)
ER.PR.sample4.from.set2<-MapSample2CELfile(data.ER.PR.sample.info,data.set2.normalized,4)
ER.PR.sample4.from.set3<-MapSample2CELfile(data.ER.PR.sample.info,data.set3.normalized,4)

head(ER.PR.sample1.from.set1)
head(ER.PR.sample1.from.set2)
head(ER.PR.sample1.from.set3)

head(ER.PR.sample3.from.set1)
head(ER.PR.sample3.from.set2)
head(ER.PR.sample3.from.set3)

head(data.set1.sample.info)
head(data.set2.sample.info)
head(data.ER.PR.sample.info)

data.set1.sample.info[,2]

data.set2.sample.info[,8:9:10]

cel.file.all<-c(unique(colnames(data.set1.normalized)),unique(colnames(data.set2.normalized)),unique(colnames(data.set3.normalized)))

length(cel.file.all)
cel.file.all

cel.file.all.2<-rbind(cbind(unique(colnames(data.set1.normalized)),rep(1,length(unique(colnames(data.set1.normalized))))),
cbind(unique(colnames(data.set2.normalized)),rep(2,length(unique(colnames(data.set2.normalized))))),
cbind(unique(colnames(data.set3.normalized)),rep(3,length(unique(colnames(data.set3.normalized))))))

cel.file.cancer.data.set2<-as.character(data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)][,8])
cel.file.cancer.data.set2.2<-as.character(data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)][,1])
cel.file.cancer.data.set2.2.2<-cel.file.cancer.data.set2.2[-which(cel.file.cancer.data.set2.2=="")]

data.set2.cancer.GSM.cel.mapping<-data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)]
data.set2.cancer.GSM.cel.mapping.2<-data.set2.cancer.GSM.cel.mapping[,c(1,8)]
data.set2.cancer.GSM.cel.mapping.3<-data.set2.cancer.GSM.cel.mapping.2[which(as.character(data.set2.cancer.GSM.cel.mapping.2[,1])!=""),]

data.set2.cancer.GSM.cel.mapping.3.2<-gsub(".cel","",as.character(data.set2.cancer.GSM.cel.mapping.3[,2]))
data.set2.cancer.GSM.cel.mapping.3.2.2<-gsub(".CEL","",as.character(data.set2.cancer.GSM.cel.mapping.3.2))
data.set2.cancer.GSM.cel.mapping.3.2.2.2<-gsub("-","_",as.character(data.set2.cancer.GSM.cel.mapping.3.2.2))
data.set2.cancer.GSM.cel.mapping.4<-cbind(data.set2.cancer.GSM.cel.mapping.3,sapply(strsplit(data.set2.cancer.GSM.cel.mapping.3.2.2.2,"_"),"[[",1))

cel.file.cancer.data.set1<-as.character(data.set1.sample.info[which(data.set1.sample.info[,4]=="cancer"),2])

cel.file.cancer.set1.set2<-rbind(cbind(cel.file.cancer.data.set1,rep(1,length(cel.file.cancer.data.set1))),
cbind(cel.file.cancer.data.set2,rep(2,length(cel.file.cancer.data.set2))),
cbind(cel.file.cancer.data.set2.2.2,rep(2,length(cel.file.cancer.data.set2.2.2))))

cel.file.cancer.set1.set2.set3<-rbind(cel.file.cancer.set1.set2,cel.file.all.2[which(cel.file.all.2[,2]==3),])

cel.file.cancer.set1.set2.set3.name<-gsub(".CEL","",cel.file.cancer.set1.set2.set3[,1])
cel.file.cancer.set1.set2.set3.name.1<-gsub(".cel","",cel.file.cancer.set1.set2.set3.name)
cel.file.cancer.set1.set2.set3.name.2<-gsub("-","_",cel.file.cancer.set1.set2.set3.name.1)

data.set123.normalized<-cbind(data.set1.normalized,data.set2.normalized,data.set3.normalized)
dim(data.set123.normalized)
dim(data.set123.normalized[,which(colnames(data.set123.normalized) %in% cel.file.cancer.set1.set2.set3[,1])])
colnames(data.set123.normalized)[grep(54,colnames(data.set123.normalized))]

original.cel.file.names<-colnames(data.set123.normalized)
original.cel.file.names.1<-gsub("X","",original.cel.file.names)
original.cel.file.names.2<-gsub(".CEL","",original.cel.file.names.1)
original.cel.file.names.3<-gsub(".cel","",original.cel.file.names.2)
original.cel.file.names.4<-gsub("\\.","_",original.cel.file.names.3)

index.4.cancer.sample<-which(original.cel.file.names.4 %in% cel.file.cancer.set1.set2.set3.name.2)
original.cel.file.name.with.mapped.names<-cbind(original.cel.file.names[index.4.cancer.sample],original.cel.file.names.4[index.4.cancer.sample])
#index.4.no.cancer.sample.1<-cbind(original.cel.file.names[-index.4.cancer.sample],original.cel.file.names.4[-index.4.cancer.sample])
setdiff(cel.file.cancer.set1.set2.set3.name.2,original.cel.file.names.4[index.4.cancer.sample])
which(original.cel.file.names.4[-index.4.cancer.sample] %in% cel.file.cancer.set1.set2.set3.name.2)

cancer.data.set123.normalized<-data.set123.normalized[,index.4.cancer.sample]
dim(cancer.data.set123.normalized)

cancer.cel.file.name.72<-colnames(cancer.data.set123.normalized)
cancer.cel.file.name.72.1<-gsub("X","",cancer.cel.file.name.72)
cancer.cel.file.name.72.2<-gsub(".CEL","",cancer.cel.file.name.72.1)
cancer.cel.file.name.72.3<-gsub(".cel","",cancer.cel.file.name.72.2)
cancer.cel.file.name.72.4<-gsub("\\.","_",cancer.cel.file.name.72.3)
cancer.cel.file.name.72.5<-cancer.cel.file.name.72.4

data.set2.cancer.GSM.cel.mapping.44<-as.character(data.set2.cancer.GSM.cel.mapping.4[which(data.set2.cancer.GSM.cel.mapping.4[,1] %in% cancer.cel.file.name.72.4),3])
cancer.cel.file.name.72.5[which(cancer.cel.file.name.72.5 %in% data.set2.cancer.GSM.cel.mapping.4[,1])]<-data.set2.cancer.GSM.cel.mapping.44
cancer.cel.file.name.72.6<-cbind(cancer.cel.file.name.72,cancer.cel.file.name.72.5)
cancer.cel.file.name.72.7<-c(sapply(strsplit(cancer.cel.file.name.72.6[1:29,2],"_"),"[[",1),
sapply(strsplit(cancer.cel.file.name.72.6[30:72,2],"_"),"[[",4))
cancer.cel.file.name.72.8<-cbind(cancer.cel.file.name.72.6,cancer.cel.file.name.72.7)

#Classify cancer patients to 4 subtypes

subtype.1<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==1),1])
subtype.2<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==2),1])
subtype.3<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==3),1])
subtype.4<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==4),1])

cancer.data.set123.normalized.st1<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.1)]
cancer.data.set123.normalized.st2<-as.data.frame(cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.2)])
colnames(cancer.data.set123.normalized.st2)<-colnames(cancer.data.set123.normalized)[which(cancer.cel.file.name.72.8[,3] %in% subtype.2)]
cancer.data.set123.normalized.st3<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.3)]
cancer.data.set123.normalized.st4<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.4)]

dim(cancer.data.set123.normalized.st1)
dim(cancer.data.set123.normalized.st2)
dim(cancer.data.set123.normalized.st3)
dim(cancer.data.set123.normalized.st4)

head(cancer.data.set123.normalized.st1)
head(cancer.data.set123.normalized.st2)
head(cancer.data.set123.normalized.st3)
head(cancer.data.set123.normalized.st4)

colnames(cancer.data.set123.normalized.st1)
colnames(cancer.data.set123.normalized.st2)
colnames(cancer.data.set123.normalized.st3)
colnames(cancer.data.set123.normalized.st4)

#n.st<-length(colnames(cancer.data.set123.normalized.st1))

cel.file.sample.infor<-as.data.frame(rbind(
cbind(colnames(cancer.data.set123.normalized.st1),rep("st1",length(colnames(cancer.data.set123.normalized.st1)))),
cbind(colnames(cancer.data.set123.normalized.st2),rep("st2",length(colnames(cancer.data.set123.normalized.st2)))),
cbind(colnames(cancer.data.set123.normalized.st3),rep("st3",length(colnames(cancer.data.set123.normalized.st3)))),
cbind(colnames(cancer.data.set123.normalized.st4),rep("st4",length(colnames(cancer.data.set123.normalized.st4))))
))
colnames(cel.file.sample.infor)=c("filename","subtype")

f <- factor(cel.file.sample.infor$subtype)
f
design <- model.matrix(~0+f)
#Differential gene expression analysis for st3 vs st4
probe<-rownames(cancer.data.set123.normalized.st1)
geneSymbol<-getSYMBOL(probe, "hugene10sttranscriptcluster.db")

#Differential gene expression analysis for st2 vs st1


#Differential gene expression analysis for st4 vs st1

save.image(file="/media/H_driver/Aimin_project/Sophia/Data_set_3.RData")
#load the affy library
.libPaths()
library(gdata)
library(affy)
library(affyio)
library(annotate)
library(plyr)

biocLite("mogene20sttranscriptcluster.db")
library(mogene20stprobeset.db)
library(mogene20sttranscriptcluster.db)

biocLite("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)

biocLite("pd.hg.u133.plus.2")
library(pd.hg.u133.plus.2)
biocLite("ggplot2")
library(ggplot2)
biocLite("reshape")
biocLite("plot3D")
biocLite("gplots")
biocLite("ggdendro")
biocLite("RColorBrewer")
biocLite("frma")

library("reshape")
library("plot3D")
library("gplots")
library(ggdendro)
library(RColorBrewer)
library(RColorBrewer)
library(frma)
biocLite("hgu133plus2frmavecs")
#R CMD INSTALL -l R/x86_64-pc-linux-gnu-library/3.2/ /media/H_driver/Aimin_project/hgu133plus2frmavecs_1.5.0.tar.gz
library(hgu133plus2frmavecs)
biocLite("Homo.sapiens")
library("Homo.sapiens")

# Read in the CEL files in the directory, then normalize the data
read.celfiles(list.celfiles("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM"))
affybatch <- read.celfiles(list.celfiles("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM"))
eset <- rma(affybatch)



Re.set1<-Function2ReadCEL_or_cel("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM")


setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM")
CEL.filenames.data.set1=dir(pattern="CEL","/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM")
data.set1 <- ReadAffy(filenames = CEL.filenames.data.set1)
ed.raw.set1 <- exprs(data.set1)
samp.set1 <- sampleNames(data.set1)
probes.set1 <- featureNames(data.set1)
ed.normalized.set1<- rma(data.set1)

setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set2_corrected")
data.set2 <- ReadAffy()
ed.raw.set2 <- exprs(data.set2)
samp.set2 <- sampleNames(data.set2)
probes.set2 <- featureNames(data.set2)
ed.normalized.set2<- rma(data.set2)

setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set3_HGSCancer_only")
data.set3 <- ReadAffy()
ed.raw.set3 <- exprs(data.set3)
samp.set3 <- sampleNames(data.set3)
probes.set3 <- featureNames(data.set3)
ed.normalized.set3<- rma(data.set3)

data.set1.normalized<-ed.normalized.set1
data.set2.normalized<-ed.normalized.set2
data.set3.normalized<-ed.normalized.set3

ll(dim=T)

data.set1.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/Microarray_Set1_FTE_HGSC_LCM/CASE_DESCRIPTION_SET1.xlsx")
data.set2.sample.info<-read.xls("/home/aiminyan/Downloads/GSE28044_George.xlsx")
data.ER.PR.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx")


cel.file.all<-c(unique(colnames(data.set1.normalized)),unique(colnames(data.set2.normalized)),unique(colnames(data.set3.normalized)))

length(cel.file.all)
cel.file.all

cel.file.all.2<-rbind(cbind(unique(colnames(data.set1.normalized)),rep(1,length(unique(colnames(data.set1.normalized))))),
cbind(unique(colnames(data.set2.normalized)),rep(2,length(unique(colnames(data.set2.normalized))))),
cbind(unique(colnames(data.set3.normalized)),rep(3,length(unique(colnames(data.set3.normalized))))))

cel.file.cancer.data.set2<-as.character(data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)][,8])
cel.file.cancer.data.set2.2<-as.character(data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)][,1])
cel.file.cancer.data.set2.2.2<-cel.file.cancer.data.set2.2[-which(cel.file.cancer.data.set2.2=="")]

data.set2.cancer.GSM.cel.mapping<-data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)]
data.set2.cancer.GSM.cel.mapping.2<-data.set2.cancer.GSM.cel.mapping[,c(1,8)]
data.set2.cancer.GSM.cel.mapping.3<-data.set2.cancer.GSM.cel.mapping.2[which(as.character(data.set2.cancer.GSM.cel.mapping.2[,1])!=""),]

data.set2.cancer.GSM.cel.mapping.3.2<-gsub(".cel","",as.character(data.set2.cancer.GSM.cel.mapping.3[,2]))
data.set2.cancer.GSM.cel.mapping.3.2.2<-gsub(".CEL","",as.character(data.set2.cancer.GSM.cel.mapping.3.2))
data.set2.cancer.GSM.cel.mapping.3.2.2.2<-gsub("-","_",as.character(data.set2.cancer.GSM.cel.mapping.3.2.2))
data.set2.cancer.GSM.cel.mapping.4<-cbind(data.set2.cancer.GSM.cel.mapping.3,sapply(strsplit(data.set2.cancer.GSM.cel.mapping.3.2.2.2,"_"),"[[",1))

cel.file.cancer.data.set1<-as.character(data.set1.sample.info[which(data.set1.sample.info[,4]=="cancer"),2])

cel.file.cancer.set1.set2<-rbind(cbind(cel.file.cancer.data.set1,rep(1,length(cel.file.cancer.data.set1))),
cbind(cel.file.cancer.data.set2,rep(2,length(cel.file.cancer.data.set2))),
cbind(cel.file.cancer.data.set2.2.2,rep(2,length(cel.file.cancer.data.set2.2.2))))

cel.file.cancer.set1.set2.set3<-rbind(cel.file.cancer.set1.set2,cel.file.all.2[which(cel.file.all.2[,2]==3),])

cel.file.cancer.set1.set2.set3.name<-gsub(".CEL","",cel.file.cancer.set1.set2.set3[,1])
cel.file.cancer.set1.set2.set3.name.1<-gsub(".cel","",cel.file.cancer.set1.set2.set3.name)
cel.file.cancer.set1.set2.set3.name.2<-gsub("-","_",cel.file.cancer.set1.set2.set3.name.1)

data.set123.normalized<-cbind(data.set1.normalized,data.set2.normalized,data.set3.normalized)

dim(data.set123.normalized)
dim(data.set123.normalized[,which(colnames(data.set123.normalized) %in% cel.file.cancer.set1.set2.set3[,1])])
colnames(data.set123.normalized)[grep(54,colnames(data.set123.normalized))]

original.cel.file.names<-colnames(data.set123.normalized)
original.cel.file.names.1<-gsub("X","",original.cel.file.names)
original.cel.file.names.2<-gsub(".CEL","",original.cel.file.names.1)
original.cel.file.names.3<-gsub(".cel","",original.cel.file.names.2)
original.cel.file.names.4<-gsub("\\.","_",original.cel.file.names.3)

index.4.cancer.sample<-which(original.cel.file.names.4 %in% cel.file.cancer.set1.set2.set3.name.2)
original.cel.file.name.with.mapped.names<-cbind(original.cel.file.names[index.4.cancer.sample],original.cel.file.names.4[index.4.cancer.sample])
#index.4.no.cancer.sample.1<-cbind(original.cel.file.names[-index.4.cancer.sample],original.cel.file.names.4[-index.4.cancer.sample])
setdiff(cel.file.cancer.set1.set2.set3.name.2,original.cel.file.names.4[index.4.cancer.sample])
which(original.cel.file.names.4[-index.4.cancer.sample] %in% cel.file.cancer.set1.set2.set3.name.2)

cancer.data.set123.normalized<-data.set123.normalized[,index.4.cancer.sample]

cancer.cel.file.name.72<-colnames(cancer.data.set123.normalized)
cancer.cel.file.name.72.1<-gsub("X","",cancer.cel.file.name.72)
cancer.cel.file.name.72.2<-gsub(".CEL","",cancer.cel.file.name.72.1)
cancer.cel.file.name.72.3<-gsub(".cel","",cancer.cel.file.name.72.2)
cancer.cel.file.name.72.4<-gsub("\\.","_",cancer.cel.file.name.72.3)
cancer.cel.file.name.72.5<-cancer.cel.file.name.72.4

data.set2.cancer.GSM.cel.mapping.44<-as.character(data.set2.cancer.GSM.cel.mapping.4[which(data.set2.cancer.GSM.cel.mapping.4[,1] %in% cancer.cel.file.name.72.4),3])
cancer.cel.file.name.72.5[which(cancer.cel.file.name.72.5 %in% data.set2.cancer.GSM.cel.mapping.4[,1])]<-data.set2.cancer.GSM.cel.mapping.44
cancer.cel.file.name.72.6<-cbind(cancer.cel.file.name.72,cancer.cel.file.name.72.5)
cancer.cel.file.name.72.7<-c(sapply(strsplit(cancer.cel.file.name.72.6[1:29,2],"_"),"[[",1),
sapply(strsplit(cancer.cel.file.name.72.6[30:72,2],"_"),"[[",4))
cancer.cel.file.name.72.8<-cbind(cancer.cel.file.name.72.6,cancer.cel.file.name.72.7)

#Classify cancer patients to 4 subtypes

subtype.1<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==1),1])
subtype.2<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==2),1])
subtype.3<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==3),1])
subtype.4<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==4),1])

cancer.data.set123.normalized.st1<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.1)]
cancer.data.set123.normalized.st2<-as.data.frame(cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.2)])
colnames(cancer.data.set123.normalized.st2)<-colnames(cancer.data.set123.normalized)[which(cancer.cel.file.name.72.8[,3] %in% subtype.2)]
cancer.data.set123.normalized.st3<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.3)]
cancer.data.set123.normalized.st4<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.4)]

#n.st<-length(colnames(cancer.data.set123.normalized.st1))

# #Use subtype1,2,3,4
# cel.file.sample.infor<-as.data.frame(rbind(
# cbind(colnames(cancer.data.set123.normalized.st1),rep("st1",length(colnames(cancer.data.set123.normalized.st1)))),
# cbind(colnames(cancer.data.set123.normalized.st2),rep("st2",length(colnames(cancer.data.set123.normalized.st2)))),
# cbind(colnames(cancer.data.set123.normalized.st3),rep("st3",length(colnames(cancer.data.set123.normalized.st3)))),
# cbind(colnames(cancer.data.set123.normalized.st4),rep("st4",length(colnames(cancer.data.set123.normalized.st4))))
# ))
# colnames(cel.file.sample.infor)=c("filename","subtype")

#Use subtype 1,3,4
cel.file.sample.infor.no.2<-as.data.frame(rbind(
  cbind(colnames(cancer.data.set123.normalized.st1),rep("st1",length(colnames(cancer.data.set123.normalized.st1)))),
  cbind(colnames(cancer.data.set123.normalized.st3),rep("st3",length(colnames(cancer.data.set123.normalized.st3)))),
  cbind(colnames(cancer.data.set123.normalized.st4),rep("st4",length(colnames(cancer.data.set123.normalized.st4))))
))

#Use frma method to do normalization
setwd("/media/H_driver/2015/Sophia/Cel_file_frma/")
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

# data.byGSym.index = ddply(tempdata.byGSym.2, c("Symbol"),function(h)
#   summary = apply(h,2,which.max)
# )

#class(data.byGSym)

#which(is.na(data.byGSym$Symbol))
rownames(data.byGSym)=data.byGSym$Symbol
data.byGSym$Symbol = NULL
data.byGSym.2<-data.byGSym
temp = apply(data.byGSym.2,2,as.numeric)
rownames(temp) = rownames(data.byGSym.2)
data.byGSym.2 = temp
SampleType<-cel.file.sample.infor.no.8$subtype
heatmap_wPCA(data.byGSym.2, "heatmap_allsample_frma.pdf","PCA_allsample_frma.pdf","/media/H_driver/2015/Sophia/", SampleType)

f.st134.frma <- factor(cel.file.sample.infor.no.8$subtype)
design.st134.frma <- model.matrix(~0+f.st134.frma)
colnames(design.st134.frma) <- levels(f.st134.frma)

colnames(data.byGSym.2)

fit.st134.gene.frma <- lmFit(data.byGSym.2, design.st134.frma)
cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134.frma)
fit2.st134.frma  <- contrasts.fit(fit.st134.gene.frma, cont.matrix.st134.frma)
fit2.st134.frma  <- eBayes(fit2.st134.frma)

TopTableSt34.gene.frma<-topTable(fit2.st134.frma,coef=1,n=60000)
TopTableSt41.gene.frma<-topTable(fit2.st134.frma,coef=2,n=60000)




OutPut2HtmlTable(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/","St41-5-frma.html","st4-vs-st1_HTML_report")
OutPut2HtmlTable(TopTableSt34.gene.frma,"/media/H_driver/2015/Sophia/","St34-5-frma.html","st3-vs-st4_HTML_report")

save.image(file="/media/H_driver/Aimin_project/Sophia/Data_set_3.RData")
#load the affy library
.libPaths()
library(gdata)
library(affy)
library(affyio)
library(annotate)
library(plyr)

biocLite("mogene20sttranscriptcluster.db")
library(mogene20stprobeset.db)
library(mogene20sttranscriptcluster.db)

biocLite("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)

biocLite("pd.hg.u133.plus.2")
library(pd.hg.u133.plus.2)
biocLite("ggplot2")
library(ggplot2)
biocLite("reshape")
biocLite("plot3D")
biocLite("gplots")
biocLite("ggdendro")
biocLite("RColorBrewer")
biocLite("frma")

library("reshape")
library("plot3D")
library("gplots")
library(ggdendro)
library(RColorBrewer)
library(RColorBrewer)
library(frma)
biocLite("hgu133plus2frmavecs")
#R CMD INSTALL -l R/x86_64-pc-linux-gnu-library/3.2/ /media/H_driver/Aimin_project/hgu133plus2frmavecs_1.5.0.tar.gz
library(hgu133plus2frmavecs)
biocLite("Homo.sapiens")
library("Homo.sapiens")

# Read in the CEL files in the directory, then normalize the data
read.celfiles(list.celfiles("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM"))
affybatch <- read.celfiles(list.celfiles("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM"))
eset <- rma(affybatch)



Re.set1<-Function2ReadCEL_or_cel("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM")


setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM")
CEL.filenames.data.set1=dir(pattern="CEL","/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM")
data.set1 <- ReadAffy(filenames = CEL.filenames.data.set1)
ed.raw.set1 <- exprs(data.set1)
samp.set1 <- sampleNames(data.set1)
probes.set1 <- featureNames(data.set1)
ed.normalized.set1<- rma(data.set1)

setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set2_corrected")
data.set2 <- ReadAffy()
ed.raw.set2 <- exprs(data.set2)
samp.set2 <- sampleNames(data.set2)
probes.set2 <- featureNames(data.set2)
ed.normalized.set2<- rma(data.set2)

setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set3_HGSCancer_only")
data.set3 <- ReadAffy()
ed.raw.set3 <- exprs(data.set3)
samp.set3 <- sampleNames(data.set3)
probes.set3 <- featureNames(data.set3)
ed.normalized.set3<- rma(data.set3)

data.set1.normalized<-ed.normalized.set1
data.set2.normalized<-ed.normalized.set2
data.set3.normalized<-ed.normalized.set3

ll(dim=T)

data.set1.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/Microarray_Set1_FTE_HGSC_LCM/CASE_DESCRIPTION_SET1.xlsx")
data.set2.sample.info<-read.xls("/home/aiminyan/Downloads/GSE28044_George.xlsx")
data.ER.PR.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx")


cel.file.all<-c(unique(colnames(data.set1.normalized)),unique(colnames(data.set2.normalized)),unique(colnames(data.set3.normalized)))

length(cel.file.all)
cel.file.all

cel.file.all.2<-rbind(cbind(unique(colnames(data.set1.normalized)),rep(1,length(unique(colnames(data.set1.normalized))))),
cbind(unique(colnames(data.set2.normalized)),rep(2,length(unique(colnames(data.set2.normalized))))),
cbind(unique(colnames(data.set3.normalized)),rep(3,length(unique(colnames(data.set3.normalized))))))

cel.file.cancer.data.set2<-as.character(data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)][,8])
cel.file.cancer.data.set2.2<-as.character(data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)][,1])
cel.file.cancer.data.set2.2.2<-cel.file.cancer.data.set2.2[-which(cel.file.cancer.data.set2.2=="")]

data.set2.cancer.GSM.cel.mapping<-data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)]
data.set2.cancer.GSM.cel.mapping.2<-data.set2.cancer.GSM.cel.mapping[,c(1,8)]
data.set2.cancer.GSM.cel.mapping.3<-data.set2.cancer.GSM.cel.mapping.2[which(as.character(data.set2.cancer.GSM.cel.mapping.2[,1])!=""),]

data.set2.cancer.GSM.cel.mapping.3.2<-gsub(".cel","",as.character(data.set2.cancer.GSM.cel.mapping.3[,2]))
data.set2.cancer.GSM.cel.mapping.3.2.2<-gsub(".CEL","",as.character(data.set2.cancer.GSM.cel.mapping.3.2))
data.set2.cancer.GSM.cel.mapping.3.2.2.2<-gsub("-","_",as.character(data.set2.cancer.GSM.cel.mapping.3.2.2))
data.set2.cancer.GSM.cel.mapping.4<-cbind(data.set2.cancer.GSM.cel.mapping.3,sapply(strsplit(data.set2.cancer.GSM.cel.mapping.3.2.2.2,"_"),"[[",1))

cel.file.cancer.data.set1<-as.character(data.set1.sample.info[which(data.set1.sample.info[,4]=="cancer"),2])

cel.file.cancer.set1.set2<-rbind(cbind(cel.file.cancer.data.set1,rep(1,length(cel.file.cancer.data.set1))),
cbind(cel.file.cancer.data.set2,rep(2,length(cel.file.cancer.data.set2))),
cbind(cel.file.cancer.data.set2.2.2,rep(2,length(cel.file.cancer.data.set2.2.2))))

cel.file.cancer.set1.set2.set3<-rbind(cel.file.cancer.set1.set2,cel.file.all.2[which(cel.file.all.2[,2]==3),])

cel.file.cancer.set1.set2.set3.name<-gsub(".CEL","",cel.file.cancer.set1.set2.set3[,1])
cel.file.cancer.set1.set2.set3.name.1<-gsub(".cel","",cel.file.cancer.set1.set2.set3.name)
cel.file.cancer.set1.set2.set3.name.2<-gsub("-","_",cel.file.cancer.set1.set2.set3.name.1)

data.set123.normalized<-cbind(data.set1.normalized,data.set2.normalized,data.set3.normalized)

dim(data.set123.normalized)
dim(data.set123.normalized[,which(colnames(data.set123.normalized) %in% cel.file.cancer.set1.set2.set3[,1])])
colnames(data.set123.normalized)[grep(54,colnames(data.set123.normalized))]

original.cel.file.names<-colnames(data.set123.normalized)
original.cel.file.names.1<-gsub("X","",original.cel.file.names)
original.cel.file.names.2<-gsub(".CEL","",original.cel.file.names.1)
original.cel.file.names.3<-gsub(".cel","",original.cel.file.names.2)
original.cel.file.names.4<-gsub("\\.","_",original.cel.file.names.3)

index.4.cancer.sample<-which(original.cel.file.names.4 %in% cel.file.cancer.set1.set2.set3.name.2)
original.cel.file.name.with.mapped.names<-cbind(original.cel.file.names[index.4.cancer.sample],original.cel.file.names.4[index.4.cancer.sample])
#index.4.no.cancer.sample.1<-cbind(original.cel.file.names[-index.4.cancer.sample],original.cel.file.names.4[-index.4.cancer.sample])
setdiff(cel.file.cancer.set1.set2.set3.name.2,original.cel.file.names.4[index.4.cancer.sample])
which(original.cel.file.names.4[-index.4.cancer.sample] %in% cel.file.cancer.set1.set2.set3.name.2)

cancer.data.set123.normalized<-data.set123.normalized[,index.4.cancer.sample]

cancer.cel.file.name.72<-colnames(cancer.data.set123.normalized)
cancer.cel.file.name.72.1<-gsub("X","",cancer.cel.file.name.72)
cancer.cel.file.name.72.2<-gsub(".CEL","",cancer.cel.file.name.72.1)
cancer.cel.file.name.72.3<-gsub(".cel","",cancer.cel.file.name.72.2)
cancer.cel.file.name.72.4<-gsub("\\.","_",cancer.cel.file.name.72.3)
cancer.cel.file.name.72.5<-cancer.cel.file.name.72.4

data.set2.cancer.GSM.cel.mapping.44<-as.character(data.set2.cancer.GSM.cel.mapping.4[which(data.set2.cancer.GSM.cel.mapping.4[,1] %in% cancer.cel.file.name.72.4),3])
cancer.cel.file.name.72.5[which(cancer.cel.file.name.72.5 %in% data.set2.cancer.GSM.cel.mapping.4[,1])]<-data.set2.cancer.GSM.cel.mapping.44
cancer.cel.file.name.72.6<-cbind(cancer.cel.file.name.72,cancer.cel.file.name.72.5)
cancer.cel.file.name.72.7<-c(sapply(strsplit(cancer.cel.file.name.72.6[1:29,2],"_"),"[[",1),
sapply(strsplit(cancer.cel.file.name.72.6[30:72,2],"_"),"[[",4))
cancer.cel.file.name.72.8<-cbind(cancer.cel.file.name.72.6,cancer.cel.file.name.72.7)

#Classify cancer patients to 4 subtypes

subtype.1<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==1),1])
subtype.2<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==2),1])
subtype.3<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==3),1])
subtype.4<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==4),1])

cancer.data.set123.normalized.st1<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.1)]
cancer.data.set123.normalized.st2<-as.data.frame(cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.2)])
colnames(cancer.data.set123.normalized.st2)<-colnames(cancer.data.set123.normalized)[which(cancer.cel.file.name.72.8[,3] %in% subtype.2)]
cancer.data.set123.normalized.st3<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.3)]
cancer.data.set123.normalized.st4<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.4)]

#n.st<-length(colnames(cancer.data.set123.normalized.st1))

# #Use subtype1,2,3,4
# cel.file.sample.infor<-as.data.frame(rbind(
# cbind(colnames(cancer.data.set123.normalized.st1),rep("st1",length(colnames(cancer.data.set123.normalized.st1)))),
# cbind(colnames(cancer.data.set123.normalized.st2),rep("st2",length(colnames(cancer.data.set123.normalized.st2)))),
# cbind(colnames(cancer.data.set123.normalized.st3),rep("st3",length(colnames(cancer.data.set123.normalized.st3)))),
# cbind(colnames(cancer.data.set123.normalized.st4),rep("st4",length(colnames(cancer.data.set123.normalized.st4))))
# ))
# colnames(cel.file.sample.infor)=c("filename","subtype")

#Use subtype 1,3,4
cel.file.sample.infor.no.2<-as.data.frame(rbind(
  cbind(colnames(cancer.data.set123.normalized.st1),rep("st1",length(colnames(cancer.data.set123.normalized.st1)))),
  cbind(colnames(cancer.data.set123.normalized.st3),rep("st3",length(colnames(cancer.data.set123.normalized.st3)))),
  cbind(colnames(cancer.data.set123.normalized.st4),rep("st4",length(colnames(cancer.data.set123.normalized.st4))))
))

#Use frma method to do normalization
setwd("/media/H_driver/2015/Sophia/Cel_file_frma/")
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
heatmap_wPCA(data.byGSym.2, "heatmap_allsample_frma.pdf","PCA_allsample_frma.pdf","/media/H_driver/2015/Sophia/", SampleType)

f.st134.frma <- factor(cel.file.sample.infor.no.8$subtype)
design.st134.frma <- model.matrix(~0+f.st134.frma)
colnames(design.st134.frma) <- levels(f.st134.frma)

fit.st134.gene.frma <- lmFit(data.byGSym.2, design.st134.frma)
cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134.frma)
fit2.st134.frma  <- contrasts.fit(fit.st134.gene.frma, cont.matrix.st134.frma)
fit2.st134.frma  <- eBayes(fit2.st134.frma)

TopTableSt34.gene.frma<-topTable(fit2.st134.frma,coef=1,n=60000)
TopTableSt41.gene.frma<-topTable(fit2.st134.frma,coef=2,n=60000)




OutPut2HtmlTable(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/","St41-5-frma.html","st4-vs-st1_HTML_report")
OutPut2HtmlTable(TopTableSt34.gene.frma,"/media/H_driver/2015/Sophia/","St34-5-frma.html","st3-vs-st4_HTML_report")

save.image(file="/media/H_driver/Aimin_project/Sophia/Data_set_3.RData")
#load the affy library
.libPaths()
library(gdata)
library(affy)
library(affyio)
library(annotate)
library(plyr)

biocLite("mogene20sttranscriptcluster.db")
library(mogene20stprobeset.db)
library(mogene20sttranscriptcluster.db)

biocLite("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)

biocLite("pd.hg.u133.plus.2")
library(pd.hg.u133.plus.2)
biocLite("ggplot2")
library(ggplot2)
biocLite("reshape")
biocLite("plot3D")
biocLite("gplots")
biocLite("ggdendro")
biocLite("RColorBrewer")
biocLite("frma")

library("reshape")
library("plot3D")
library("gplots")
library(ggdendro)
library(RColorBrewer)
library(RColorBrewer)
library(frma)
biocLite("hgu133plus2frmavecs")
#R CMD INSTALL -l R/x86_64-pc-linux-gnu-library/3.2/ /media/H_driver/Aimin_project/hgu133plus2frmavecs_1.5.0.tar.gz
library(hgu133plus2frmavecs)
biocLite("Homo.sapiens")
library("Homo.sapiens")

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

#Use frma method to do normalization
setwd("/media/H_driver/2015/Sophia/Cel_file_frma/")
data.set123.raw <- ReadAffy()
data.set123.frma<-frma(data.set123.raw)
ed.set123.frma <- exprs(data.set123.frma)
samp.frma <- sampleNames(data.set123.frma)
probes.frma <- featureNames(data.set123.frma)

samp.frma.1<-gsub(" - ","_",samp.frma)
samp.frma.2<-gsub(".CEL","",samp.frma.1)
samp.frma.3<-gsub(".cel","",samp.frma.2)
samp.frma.4<-gsub("-","_",samp.frma.3)
samp.frma.5<-cbind(samp.frma,samp.frma.4)
samp.frma.cancer.index<-which(samp.frma.5[,2] %in% cel.file.name.key.set.symbol.st134[,1])

samp.frma.6<-samp.frma.5[samp.frma.cancer.index,]

data.set123.frma.cancer<-data.set123.frma[,samp.frma.cancer.index]
ed.set123.frma.cancer<-ed.set123.frma[,samp.frma.cancer.index]

samp.frma.7<-cbind(samp.frma.6,colnames(ed.set123.frma.cancer))
samp.frma.8<-samp.frma.7[,2:3]
colnames(samp.frma.8)<-c("key","filename_frma")

cel.file.sample.infor.no.8<-merge(cel.file.name.key.set.symbol.st134,samp.frma.8,by="key",sort=FALSE)

GeneSym.all <- getSYMBOL(rownames(ed.set123.frma.cancer), "hgu133plus2.db")
ndata<-ed.set123.frma.cancer
geneSymbol=GeneSym.all
tempdata.byGSym = data.frame(ndata, Symbol = geneSymbol)
tempdata.byGSym.2<-tempdata.byGSym[-which(is.na(tempdata.byGSym$Symbol)),]
rownames(tempdata.byGSym.2) = NULL
data.byGSym = ddply(tempdata.byGSym.2, c("Symbol"),function(h)
  summary = apply(h,2,max)
)

rownames(data.byGSym)=data.byGSym$Symbol
data.byGSym.2<-data.byGSym[,-dim(data.byGSym)[2]]
temp = apply(data.byGSym.2,2,as.numeric)
rownames(temp) = rownames(data.byGSym.2)
data.byGSym.2 = temp
SampleType<-cel.file.sample.infor.no.8$subtype
heatmap_wPCA(data.byGSym.2, "heatmap_allsample_frma_2.pdf","PCA_allsample_frma_2.pdf","/media/H_driver/2015/Sophia/Results/", SampleType)

f.st134.frma <- factor(cel.file.sample.infor.no.8$subtype)
design.st134.frma <- model.matrix(~0+f.st134.frma)
colnames(design.st134.frma) <- levels(f.st134.frma)

fit.st134.gene.frma <- lmFit(data.byGSym.2, design.st134.frma)
cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134.frma)
fit2.st134.frma  <- contrasts.fit(fit.st134.gene.frma, cont.matrix.st134.frma)
fit2.st134.frma  <- eBayes(fit2.st134.frma)

TopTableSt34.gene.frma<-topTable(fit2.st134.frma,coef=1,n=dim(data.byGSym.2)[1])
TopTableSt41.gene.frma<-topTable(fit2.st134.frma,coef=2,n=dim(data.byGSym.2)[1])



OutPut2HtmlTable(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/Results/","St41-6-frma.html","st4-vs-st1_HTML_report")
OutPut2HtmlTable(TopTableSt34.gene.frma,"/media/H_driver/2015/Sophia/Results/","St34-6-frma.html","st3-vs-st4_HTML_report")

#Generate the rank file for GSEA



GenerateRankFile4GSEA(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/Results/","St41-4_GSEA_up_down_300.rnk",300)

#write.csv(TopTableSt41.gene.frma, file="/media/H_driver/2015/Sophia/Results/St41.csv",quote =FALSE)

save.image(file="/media/H_driver/Aimin_project/Sophia/Data_set_3_13_2016.RData")
savehistory(file="/media/H_driver/Aimin_project/Sophia/Data_set_3_13_2016.Rhistory")#load the affy library
.libPaths()
library(gdata)
library(affy)
library(affyio)
library(annotate)
library(plyr)

biocLite("mogene20sttranscriptcluster.db")
library(mogene20stprobeset.db)
library(mogene20sttranscriptcluster.db)

biocLite("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)

biocLite("pd.hg.u133.plus.2")
library(pd.hg.u133.plus.2)
biocLite("ggplot2")
library(ggplot2)
biocLite("reshape")
biocLite("plot3D")
biocLite("gplots")
biocLite("ggdendro")
biocLite("RColorBrewer")
biocLite("frma")

library("reshape")
library("plot3D")
library("gplots")
library(ggdendro)
library(RColorBrewer)
library(RColorBrewer)
library(frma)
biocLite("hgu133plus2frmavecs")
#R CMD INSTALL -l R/x86_64-pc-linux-gnu-library/3.2/ /media/H_driver/Aimin_project/hgu133plus2frmavecs_1.5.0.tar.gz
library(hgu133plus2frmavecs)
biocLite("Homo.sapiens")
library("Homo.sapiens")

# Read in the CEL files in the directory, then normalize the data
Function2ReadCEL_or_cel<-function(input_dir){
   affybatch <- read.celfiles(list.celfiles(setwd(input_dir)))
   eset <- rma(affybatch)
   re<-list(aff_data=affybatch,ed_norma=eset)
   return(re)
 }
Re.set1<-Function2ReadCEL_or_cel("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM")
Re.set2<-Function2ReadCEL_or_cel("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set2_corrected")

setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM")
data.set1 <- ReadAffy()
ed.raw.set1 <- exprs(data.set1)
samp.set1 <- sampleNames(data.set1)
probes.set1 <- featureNames(data.set1)
ed.normalized.set1<- rma(data.set1)

setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set2_corrected")
data.set2 <- ReadAffy()
ed.raw.set2 <- exprs(data.set2)
samp.set2 <- sampleNames(data.set2)
probes.set2 <- featureNames(data.set2)
ed.normalized.set2<- rma(data.set2)

setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set3_HGSCancer_only")
data.set3 <- ReadAffy()

ed.raw.set3 <- exprs(data.set3)
samp.set3 <- sampleNames(data.set3)
probes.set3 <- featureNames(data.set3)
ed.normalized.set3<- rma(data.set3)

data.set1.normalized<-ed.normalized.set1
data.set2.normalized<-ed.normalized.set2
data.set3.normalized<-ed.normalized.set3

ll(dim=T)

data.set1.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/Microarray_Set1_FTE_HGSC_LCM/CASE_DESCRIPTION_SET1.xlsx")
data.set2.sample.info<-read.xls("/home/aiminyan/Downloads/GSE28044_George.xlsx")
data.ER.PR.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx")


cel.file.all<-c(unique(colnames(data.set1.normalized)),unique(colnames(data.set2.normalized)),unique(colnames(data.set3.normalized)))

length(cel.file.all)
cel.file.all

cel.file.all.2<-rbind(cbind(unique(colnames(data.set1.normalized)),rep(1,length(unique(colnames(data.set1.normalized))))),
cbind(unique(colnames(data.set2.normalized)),rep(2,length(unique(colnames(data.set2.normalized))))),
cbind(unique(colnames(data.set3.normalized)),rep(3,length(unique(colnames(data.set3.normalized))))))

cel.file.cancer.data.set2<-as.character(data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)][,8])
cel.file.cancer.data.set2.2<-as.character(data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)][,1])
cel.file.cancer.data.set2.2.2<-cel.file.cancer.data.set2.2[-which(cel.file.cancer.data.set2.2=="")]

data.set2.cancer.GSM.cel.mapping<-data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)]
data.set2.cancer.GSM.cel.mapping.2<-data.set2.cancer.GSM.cel.mapping[,c(1,8)]
data.set2.cancer.GSM.cel.mapping.3<-data.set2.cancer.GSM.cel.mapping.2[which(as.character(data.set2.cancer.GSM.cel.mapping.2[,1])!=""),]

data.set2.cancer.GSM.cel.mapping.3.2<-gsub(".cel","",as.character(data.set2.cancer.GSM.cel.mapping.3[,2]))
data.set2.cancer.GSM.cel.mapping.3.2.2<-gsub(".CEL","",as.character(data.set2.cancer.GSM.cel.mapping.3.2))
data.set2.cancer.GSM.cel.mapping.3.2.2.2<-gsub("-","_",as.character(data.set2.cancer.GSM.cel.mapping.3.2.2))
data.set2.cancer.GSM.cel.mapping.4<-cbind(data.set2.cancer.GSM.cel.mapping.3,sapply(strsplit(data.set2.cancer.GSM.cel.mapping.3.2.2.2,"_"),"[[",1))

cel.file.cancer.data.set1<-as.character(data.set1.sample.info[which(data.set1.sample.info[,4]=="cancer"),2])

cel.file.cancer.set1.set2<-rbind(cbind(cel.file.cancer.data.set1,rep(1,length(cel.file.cancer.data.set1))),
cbind(cel.file.cancer.data.set2,rep(2,length(cel.file.cancer.data.set2))),
cbind(cel.file.cancer.data.set2.2.2,rep(2,length(cel.file.cancer.data.set2.2.2))))

cel.file.cancer.set1.set2.set3<-rbind(cel.file.cancer.set1.set2,cel.file.all.2[which(cel.file.all.2[,2]==3),])

cel.file.cancer.set1.set2.set3.name<-gsub(".CEL","",cel.file.cancer.set1.set2.set3[,1])
cel.file.cancer.set1.set2.set3.name.1<-gsub(".cel","",cel.file.cancer.set1.set2.set3.name)
cel.file.cancer.set1.set2.set3.name.2<-gsub("-","_",cel.file.cancer.set1.set2.set3.name.1)

data.set123.normalized<-cbind(data.set1.normalized,data.set2.normalized,data.set3.normalized)

dim(data.set123.normalized)
dim(data.set123.normalized[,which(colnames(data.set123.normalized) %in% cel.file.cancer.set1.set2.set3[,1])])
colnames(data.set123.normalized)[grep(54,colnames(data.set123.normalized))]

original.cel.file.names<-colnames(data.set123.normalized)
original.cel.file.names.1<-gsub("X","",original.cel.file.names)
original.cel.file.names.2<-gsub(".CEL","",original.cel.file.names.1)
original.cel.file.names.3<-gsub(".cel","",original.cel.file.names.2)
original.cel.file.names.4<-gsub("\\.","_",original.cel.file.names.3)

index.4.cancer.sample<-which(original.cel.file.names.4 %in% cel.file.cancer.set1.set2.set3.name.2)
original.cel.file.name.with.mapped.names<-cbind(original.cel.file.names[index.4.cancer.sample],original.cel.file.names.4[index.4.cancer.sample])
#index.4.no.cancer.sample.1<-cbind(original.cel.file.names[-index.4.cancer.sample],original.cel.file.names.4[-index.4.cancer.sample])
setdiff(cel.file.cancer.set1.set2.set3.name.2,original.cel.file.names.4[index.4.cancer.sample])
which(original.cel.file.names.4[-index.4.cancer.sample] %in% cel.file.cancer.set1.set2.set3.name.2)

cancer.data.set123.normalized<-data.set123.normalized[,index.4.cancer.sample]

cancer.cel.file.name.72<-colnames(cancer.data.set123.normalized)
cancer.cel.file.name.72.1<-gsub("X","",cancer.cel.file.name.72)
cancer.cel.file.name.72.2<-gsub(".CEL","",cancer.cel.file.name.72.1)
cancer.cel.file.name.72.3<-gsub(".cel","",cancer.cel.file.name.72.2)
cancer.cel.file.name.72.4<-gsub("\\.","_",cancer.cel.file.name.72.3)
cancer.cel.file.name.72.5<-cancer.cel.file.name.72.4

data.set2.cancer.GSM.cel.mapping.44<-as.character(data.set2.cancer.GSM.cel.mapping.4[which(data.set2.cancer.GSM.cel.mapping.4[,1] %in% cancer.cel.file.name.72.4),3])
cancer.cel.file.name.72.5[which(cancer.cel.file.name.72.5 %in% data.set2.cancer.GSM.cel.mapping.4[,1])]<-data.set2.cancer.GSM.cel.mapping.44
cancer.cel.file.name.72.6<-cbind(cancer.cel.file.name.72,cancer.cel.file.name.72.5)
cancer.cel.file.name.72.7<-c(sapply(strsplit(cancer.cel.file.name.72.6[1:29,2],"_"),"[[",1),
sapply(strsplit(cancer.cel.file.name.72.6[30:72,2],"_"),"[[",4))
cancer.cel.file.name.72.8<-cbind(cancer.cel.file.name.72.6,cancer.cel.file.name.72.7)

#Classify cancer patients to 4 subtypes

subtype.1<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==1),1])
subtype.2<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==2),1])
subtype.3<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==3),1])
subtype.4<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==4),1])

cancer.data.set123.normalized.st1<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.1)]
cancer.data.set123.normalized.st2<-as.data.frame(cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.2)])
colnames(cancer.data.set123.normalized.st2)<-colnames(cancer.data.set123.normalized)[which(cancer.cel.file.name.72.8[,3] %in% subtype.2)]
cancer.data.set123.normalized.st3<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.3)]
cancer.data.set123.normalized.st4<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.4)]

#n.st<-length(colnames(cancer.data.set123.normalized.st1))

# #Use subtype1,2,3,4
# cel.file.sample.infor<-as.data.frame(rbind(
# cbind(colnames(cancer.data.set123.normalized.st1),rep("st1",length(colnames(cancer.data.set123.normalized.st1)))),
# cbind(colnames(cancer.data.set123.normalized.st2),rep("st2",length(colnames(cancer.data.set123.normalized.st2)))),
# cbind(colnames(cancer.data.set123.normalized.st3),rep("st3",length(colnames(cancer.data.set123.normalized.st3)))),
# cbind(colnames(cancer.data.set123.normalized.st4),rep("st4",length(colnames(cancer.data.set123.normalized.st4))))
# ))
# colnames(cel.file.sample.infor)=c("filename","subtype")

#Use subtype 1,3,4
cel.file.sample.infor.no.2<-as.data.frame(rbind(
  cbind(colnames(cancer.data.set123.normalized.st1),rep("st1",length(colnames(cancer.data.set123.normalized.st1)))),
  cbind(colnames(cancer.data.set123.normalized.st3),rep("st3",length(colnames(cancer.data.set123.normalized.st3)))),
  cbind(colnames(cancer.data.set123.normalized.st4),rep("st4",length(colnames(cancer.data.set123.normalized.st4))))
))

#Use frma method to do normalization
setwd("/media/H_driver/2015/Sophia/Cel_file_frma/")
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
heatmap_wPCA(data.byGSym.2, "heatmap_allsample_frma.pdf","PCA_allsample_frma.pdf","/media/H_driver/2015/Sophia/", SampleType)

f.st134.frma <- factor(cel.file.sample.infor.no.8$subtype)
design.st134.frma <- model.matrix(~0+f.st134.frma)
colnames(design.st134.frma) <- levels(f.st134.frma)

fit.st134.gene.frma <- lmFit(data.byGSym.2, design.st134.frma)
cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134.frma)
fit2.st134.frma  <- contrasts.fit(fit.st134.gene.frma, cont.matrix.st134.frma)
fit2.st134.frma  <- eBayes(fit2.st134.frma)

TopTableSt34.gene.frma<-topTable(fit2.st134.frma,coef=1,n=60000)
TopTableSt41.gene.frma<-topTable(fit2.st134.frma,coef=2,n=60000)


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

OutPut2HtmlTable(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/","St41-5-frma.html","st4-vs-st1_HTML_report")
OutPut2HtmlTable(TopTableSt34.gene.frma,"/media/H_driver/2015/Sophia/","St34-5-frma.html","st3-vs-st4_HTML_report")

save.image(file="/media/H_driver/Aimin_project/Sophia/Data_set_3.RData")
#load the affy library
.libPaths()
library(gdata)
library(affy)
library(affyio)
library(annotate)
library(plyr)

#biocLite("mogene20sttranscriptcluster.db")
#library(mogene20stprobeset.db)
#library(mogene20sttranscriptcluster.db)

#biocLite("hugene10sttranscriptcluster.db")
#library(hugene10sttranscriptcluster.db)

biocLite("pd.hg.u133.plus.2")
library(pd.hg.u133.plus.2)
biocLite("ggplot2")
library(ggplot2)
biocLite("reshape")
biocLite("plot3D")
biocLite("gplots")
biocLite("ggdendro")
biocLite("RColorBrewer")
biocLite("frma")

library("reshape")
library("plot3D")
library("gplots")
library(ggdendro)
library(RColorBrewer)
library(RColorBrewer)
library(frma)
biocLite("hgu133plus2frmavecs")
#R CMD INSTALL -l R/x86_64-pc-linux-gnu-library/3.2/ /media/H_driver/Aimin_project/hgu133plus2frmavecs_1.5.0.tar.gz
library(hgu133plus2frmavecs)
biocLite("Homo.sapiens")
library("Homo.sapiens")

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

ReformatCelName<-function(data.set2.cancer.GSM.cel.mapping.3){
data.set2.cancer.GSM.cel.mapping.3.2<-gsub(".cel","",data.set2.cancer.GSM.cel.mapping.3)
data.set2.cancer.GSM.cel.mapping.3.2.2<-gsub(".CEL","",as.character(data.set2.cancer.GSM.cel.mapping.3.2))
data.set2.cancer.GSM.cel.mapping.3.2.2.2<-gsub("-","_",as.character(data.set2.cancer.GSM.cel.mapping.3.2.2))
data.set2.cancer.GSM.cel.mapping.4<-cbind(data.set2.cancer.GSM.cel.mapping.3,data.set2.cancer.GSM.cel.mapping.3.2.2.2,sapply(strsplit(data.set2.cancer.GSM.cel.mapping.3.2.2.2,"_"),"[[",1))
return(data.set2.cancer.GSM.cel.mapping.4)
}

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

#Use frma method to do normalization
setwd("/media/H_driver/2015/Sophia/Cel_file_frma/")
data.set123.raw <- ReadAffy()
data.set123.frma<-frma(data.set123.raw)
ed.set123.frma <- exprs(data.set123.frma)
samp.frma <- sampleNames(data.set123.frma)
probes.frma <- featureNames(data.set123.frma)

samp.frma.1<-gsub(" - ","_",samp.frma)
samp.frma.2<-gsub(".CEL","",samp.frma.1)
samp.frma.3<-gsub(".cel","",samp.frma.2)
samp.frma.4<-gsub("-","_",samp.frma.3)
samp.frma.5<-cbind(samp.frma,samp.frma.4)
samp.frma.cancer.index<-which(samp.frma.5[,2] %in% cel.file.name.key.set.symbol.st134[,1])

samp.frma.6<-samp.frma.5[samp.frma.cancer.index,]

data.set123.frma.cancer<-data.set123.frma[,samp.frma.cancer.index]
ed.set123.frma.cancer<-ed.set123.frma[,samp.frma.cancer.index]
dim(ed.set123.frma.cancer)

GenerateFiles4GSEA<-function(ed.set123.frma.cancer,out_dir,out_file_name,subtypeA,subtypeB){

index.st4<-which(colnames(ed.set123.frma.cancer) %in% cel.file.name.key.set.symbol.st134[which(cel.file.name.key.set.symbol.st134[,5]==subtypeA),2])
index.st1<-which(colnames(ed.set123.frma.cancer) %in% cel.file.name.key.set.symbol.st134[which(cel.file.name.key.set.symbol.st134[,5]==subtypeB),2])

ed.set123.frma.cancer.st41<-cbind(ed.set123.frma.cancer[,index.st4],ed.set123.frma.cancer[,index.st1])
cat(subtypeA,"\t",length(index.st4),subtypeB,"\t",length(index.st1),"\t",dim(ed.set123.frma.cancer.st41),"\n")

ed.set123.frma.cancer.st41.reformat<-cbind(rownames(ed.set123.frma.cancer.st41),rep("NA",dim(ed.set123.frma.cancer.st41)[1]),ed.set123.frma.cancer.st41)
colnames(ed.set123.frma.cancer.st41.reformat)[c(1,2)]=c("NAME","Description")

write.table(ed.set123.frma.cancer.st41.reformat,file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = TRUE)

}

GenerateFiles4GSEA(ed.set123.frma.cancer,"/media/H_driver/2015/Sophia/Results/","St41-4_GSEA.xls","st4","st1")
GenerateFiles4GSEA(ed.set123.frma.cancer,"/media/H_driver/2015/Sophia/Results/","St34-4_GSEA.xls","st3","st4")

samp.frma.7<-cbind(samp.frma.6,colnames(ed.set123.frma.cancer))
samp.frma.8<-samp.frma.7[,2:3]
colnames(samp.frma.8)<-c("key","filename_frma")

cel.file.sample.infor.no.8<-merge(cel.file.name.key.set.symbol.st134,samp.frma.8,by="key",sort=FALSE)

GeneSym.all <- getSYMBOL(rownames(ed.set123.frma.cancer), "hgu133plus2.db")
ndata<-ed.set123.frma.cancer
geneSymbol=GeneSym.all
tempdata.byGSym = data.frame(ndata, Symbol = geneSymbol)
tempdata.byGSym.2<-tempdata.byGSym[-which(is.na(tempdata.byGSym$Symbol)),]
rownames(tempdata.byGSym.2) = NULL
data.byGSym = ddply(tempdata.byGSym.2, c("Symbol"),function(h)
  summary = apply(h,2,max)
)

rownames(data.byGSym)=data.byGSym$Symbol
data.byGSym.2<-data.byGSym[,-dim(data.byGSym)[2]]
temp = apply(data.byGSym.2,2,as.numeric)
rownames(temp) = rownames(data.byGSym.2)
data.byGSym.2 = temp
SampleType<-cel.file.sample.infor.no.8$subtype
heatmap_wPCA(data.byGSym.2, "heatmap_allsample_frma_2.pdf","PCA_allsample_frma_2.pdf","/media/H_driver/2015/Sophia/Results/", SampleType)

f.st134.frma <- factor(cel.file.sample.infor.no.8$subtype)
design.st134.frma <- model.matrix(~0+f.st134.frma)
colnames(design.st134.frma) <- levels(f.st134.frma)

fit.st134.gene.frma <- lmFit(data.byGSym.2, design.st134.frma)
cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134.frma)
fit2.st134.frma  <- contrasts.fit(fit.st134.gene.frma, cont.matrix.st134.frma)
fit2.st134.frma  <- eBayes(fit2.st134.frma)

TopTableSt34.gene.frma<-topTable(fit2.st134.frma,coef=1,n=dim(data.byGSym.2)[1])
TopTableSt41.gene.frma<-topTable(fit2.st134.frma,coef=2,n=dim(data.byGSym.2)[1])

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

OutPut2HtmlTable(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/Results/","St41-6-frma.html","st4-vs-st1_HTML_report")
OutPut2HtmlTable(TopTableSt34.gene.frma,"/media/H_driver/2015/Sophia/Results/","St34-6-frma.html","st3-vs-st4_HTML_report")

#Generate the rank file for GSEA

GenerateRankFile4GSEA<-function(TopTableSt41.gene.frma,out_dir,out_file_name,up_down_top){

  FC.sign=sign(TopTableSt41.gene.frma[,1])
  p.inve=1/TopTableSt41.gene.frma[,4]
  Gene.Rank.Score=FC.sign*p.inve

  TopTableSt41.gene.frma.GSEA<-cbind(TopTableSt41.gene.frma,Gene.Rank.Score)

  TopTableSt41.gene.frma.GSEA.sorted.by.score<-TopTableSt41.gene.frma.GSEA[order(TopTableSt41.gene.frma.GSEA[,7],decreasing = TRUE),]

  Re.GSEA<-cbind(rownames(TopTableSt41.gene.frma.GSEA.sorted.by.score),TopTableSt41.gene.frma.GSEA.sorted.by.score[,7])
  colnames(Re.GSEA)=c("GeneName","Score")

  Re.GSEA.2<-rbind(head(Re.GSEA,up_down_top),tail(Re.GSEA,up_down_top))

  write.table(Re.GSEA.2, file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = FALSE)

  return(TopTableSt41.gene.frma.GSEA.sorted.by.score)
}

GenerateRankFile4GSEA(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/Results/","St41-4_GSEA_up_down_300.rnk",300)

#write.csv(TopTableSt41.gene.frma, file="/media/H_driver/2015/Sophia/Results/St41.csv",quote =FALSE)

save.image(file="/media/H_driver/Aimin_project/Sophia/Data_set_3_13_2016.RData")
savehistory(file="/media/H_driver/Aimin_project/Sophia/Data_set_3_13_2016.Rhistory")#load the affy library
.libPaths()
library(gdata)
library(affy)
library(affyio)
library(annotate)
library(plyr)

#biocLite("mogene20sttranscriptcluster.db")
#library(mogene20stprobeset.db)
#library(mogene20sttranscriptcluster.db)

#biocLite("hugene10sttranscriptcluster.db")
#library(hugene10sttranscriptcluster.db)

biocLite("pd.hg.u133.plus.2")
library(pd.hg.u133.plus.2)
biocLite("ggplot2")
library(ggplot2)
biocLite("reshape")
biocLite("plot3D")
biocLite("gplots")
biocLite("ggdendro")
biocLite("RColorBrewer")
biocLite("frma")

library("reshape")
library("plot3D")
library("gplots")
library(ggdendro)
library(RColorBrewer)
library(RColorBrewer)
library(frma)
biocLite("hgu133plus2frmavecs")
#R CMD INSTALL -l R/x86_64-pc-linux-gnu-library/3.2/ /media/H_driver/Aimin_project/hgu133plus2frmavecs_1.5.0.tar.gz
library(hgu133plus2frmavecs)
biocLite("Homo.sapiens")
library("Homo.sapiens")
library(hgu133plus2.db)

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

ReformatCelName<-function(data.set2.cancer.GSM.cel.mapping.3){
data.set2.cancer.GSM.cel.mapping.3.2<-gsub(".cel","",data.set2.cancer.GSM.cel.mapping.3)
data.set2.cancer.GSM.cel.mapping.3.2.2<-gsub(".CEL","",as.character(data.set2.cancer.GSM.cel.mapping.3.2))
data.set2.cancer.GSM.cel.mapping.3.2.2.2<-gsub("-","_",as.character(data.set2.cancer.GSM.cel.mapping.3.2.2))
data.set2.cancer.GSM.cel.mapping.4<-cbind(data.set2.cancer.GSM.cel.mapping.3,data.set2.cancer.GSM.cel.mapping.3.2.2.2,sapply(strsplit(data.set2.cancer.GSM.cel.mapping.3.2.2.2,"_"),"[[",1))
return(data.set2.cancer.GSM.cel.mapping.4)
}

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

#Use frma method to do normalization
setwd("/media/H_driver/2015/Sophia/Cel_file_frma/")
data.set123.raw <- ReadAffy()
data.set123.frma<-frma(data.set123.raw)
ed.set123.frma <- exprs(data.set123.frma)
samp.frma <- sampleNames(data.set123.frma)
probes.frma <- featureNames(data.set123.frma)

samp.frma.1<-gsub(" - ","_",samp.frma)
samp.frma.2<-gsub(".CEL","",samp.frma.1)
samp.frma.3<-gsub(".cel","",samp.frma.2)
samp.frma.4<-gsub("-","_",samp.frma.3)
samp.frma.5<-cbind(samp.frma,samp.frma.4)
samp.frma.cancer.index<-which(samp.frma.5[,2] %in% cel.file.name.key.set.symbol.st134[,1])
samp.frma.6<-samp.frma.5[samp.frma.cancer.index,]

data.set123.frma.cancer<-data.set123.frma[,samp.frma.cancer.index]
ed.set123.frma.cancer<-ed.set123.frma[,samp.frma.cancer.index]
dim(ed.set123.frma.cancer)

samp.frma.7<-cbind(samp.frma.6,colnames(ed.set123.frma.cancer))
samp.frma.8<-samp.frma.7[,2:3]
colnames(samp.frma.8)<-c("key","filename_frma")

cel.file.sample.infor.no.8<-merge(cel.file.name.key.set.symbol.st134,samp.frma.8,by="key",sort=FALSE)

GeneSym.all <- getSYMBOL(rownames(ed.set123.frma.cancer), "hgu133plus2.db")


#<-cbind(names(GeneSym.all),GeneSym.all)

# #Collapse probes into genes
ndata<-ed.set123.frma.cancer
geneSymbol=GeneSym.all
tempdata.byGSym =data.frame(ndata, Symbol = geneSymbol)

colnames(tempdata.byGSym)


tempdata.byGSym.2<-tempdata.byGSym[-which(is.na(tempdata.byGSym$Symbol)),]

test.sample<-ed.set123.frma.cancer[which(GeneSym.all=="ERBB4"),which(colnames(ed.set123.frma.cancer) %in% cel.file.sample.infor.no.8[which(cel.file.sample.infor.no.8[,5]=="st4"),6])]

test.sample.2<-data.frame(test.sample,Symbol="ERBB4")
test.sample.3 = test.sample.2
test.sample.3[,1:5] = apply(test.sample.3[,1:5],2,as.numeric)


data.byGSym = ddply(tempdata.byGSym.2, c("Symbol"),function(h)
  summary = apply(h[,1:(dim(tempdata.byGSym.2)[2]-1)],2,max)
)


rownames(data.byGSym)=data.byGSym$Symbol
colnames(data.byGSym)

data.byGSym.2<-data.byGSym[,-1]
temp = apply(data.byGSym.2,2,as.numeric)
rownames(temp) = rownames(data.byGSym.2)
colnames(temp) = colnames(ndata)

data.byGSym.2 = temp
colnames(data.byGSym.2)
#data.byGSym.2 = ed.set123.frma.cancer

SampleType<-cel.file.sample.infor.no.8$subtype
heatmap_wPCA(data.byGSym.2, "heatmap_allsample_frma_all_probes.pdf","PCA_allsample_frma_all_probes.pdf","/media/H_driver/2015/Sophia/Results/", SampleType)

f.st134.frma <- factor(cel.file.sample.infor.no.8$subtype)
design.st134.frma <- model.matrix(~0+f.st134.frma)
colnames(design.st134.frma) <- levels(f.st134.frma)

reorder.index<-match(cel.file.sample.infor.no.8[,6],colnames(data.byGSym.2))

data.byGSym.3<-data.byGSym.2[,reorder.index]

colnames(data.byGSym.3)


fit.st134.gene.frma <- lmFit(data.byGSym.3, design.st134.frma)
cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134.frma)
fit2.st134.frma  <- contrasts.fit(fit.st134.gene.frma, cont.matrix.st134.frma)
fit2.st134.frma  <- eBayes(fit2.st134.frma)

TopTableSt34.gene.frma<-topTable(fit2.st134.frma,coef=1,n=dim(data.byGSym.2)[1])
TopTableSt41.gene.frma<-topTable(fit2.st134.frma,coef=2,n=dim(data.byGSym.2)[1])

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

OutPut2HtmlTable(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/Results/","St41-7-frma.html","st4-vs-st1_HTML_report")
OutPut2HtmlTable(TopTableSt34.gene.frma,"/media/H_driver/2015/Sophia/Results/","St34-6-frma.html","st3-vs-st4_HTML_report")

#Generate the rank file for GSEA
GenerateRankFile4GSEA<-function(TopTableSt41.gene.frma,out_dir,out_file_name,up_down_top){

  FC.sign=sign(TopTableSt41.gene.frma[,1])
  p.inve=1/TopTableSt41.gene.frma[,4]
  Gene.Rank.Score=FC.sign*p.inve

  TopTableSt41.gene.frma.GSEA<-cbind(TopTableSt41.gene.frma,Gene.Rank.Score)

  TopTableSt41.gene.frma.GSEA.sorted.by.score<-TopTableSt41.gene.frma.GSEA[order(TopTableSt41.gene.frma.GSEA[,7],decreasing = TRUE),]

  Re.GSEA<-cbind(rownames(TopTableSt41.gene.frma.GSEA.sorted.by.score),TopTableSt41.gene.frma.GSEA.sorted.by.score[,7])
  colnames(Re.GSEA)=c("GeneName","Score")

  Re.GSEA.2<-rbind(head(Re.GSEA,up_down_top),tail(Re.GSEA,up_down_top))

  write.table(Re.GSEA.2, file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = FALSE)

  return(TopTableSt41.gene.frma.GSEA.sorted.by.score)
}

GenerateRankFile4GSEA(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/Results/","St41-4_GSEA_up_down_300.rnk",300)

#write.csv(TopTableSt41.gene.frma, file="/media/H_driver/2015/Sophia/Results/St41.csv",quote =FALSE)

#Use all data to perform GSEA analysis
  GenerateFiles4GSEA<-function(ed.set123.frma.cancer,out_dir,out_file_name,subtypeA,subtypeB){

    index.st4<-which(colnames(ed.set123.frma.cancer) %in% cel.file.name.key.set.symbol.st134[which(cel.file.name.key.set.symbol.st134[,5]==subtypeA),2])
    index.st1<-which(colnames(ed.set123.frma.cancer) %in% cel.file.name.key.set.symbol.st134[which(cel.file.name.key.set.symbol.st134[,5]==subtypeB),2])

    ed.set123.frma.cancer.st41<-cbind(ed.set123.frma.cancer[,index.st4],ed.set123.frma.cancer[,index.st1])
    cat(subtypeA,"\t",length(index.st4),subtypeB,"\t",length(index.st1),"\t",dim(ed.set123.frma.cancer.st41),"\n")

    ed.set123.frma.cancer.st41.reformat<-cbind(rownames(ed.set123.frma.cancer.st41),rep("NA",dim(ed.set123.frma.cancer.st41)[1]),ed.set123.frma.cancer.st41)
    colnames(ed.set123.frma.cancer.st41.reformat)[c(1,2)]=c("NAME","Description")

    write.table(ed.set123.frma.cancer.st41.reformat,file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = TRUE)

  }




GenerateFiles4GSEA(ed.set123.frma.cancer,"/media/H_driver/2015/Sophia/Results/","St41-4_GSEA.xls","st4","st1")
GenerateFiles4GSEA(ed.set123.frma.cancer,"/media/H_driver/2015/Sophia/Results/","St34-4_GSEA.xls","st3","st4")

save.image(file="/media/H_driver/Aimin_project/Sophia/Data_set_3_13_2016.RData")
#/Volumes/Bioinformatics$/2015/Sophia/

savehistory(file="/media/H_driver/Aimin_project/Sophia/Data_set_3_13_2016.Rhistory")
#load the affy library
.libPaths()
library(gdata)
library(affy)
library(affyio)
library(annotate)
library(plyr)

#biocLite("mogene20sttranscriptcluster.db")
#library(mogene20stprobeset.db)
#library(mogene20sttranscriptcluster.db)

#biocLite("hugene10sttranscriptcluster.db")
#library(hugene10sttranscriptcluster.db)

biocLite("pd.hg.u133.plus.2")
library(pd.hg.u133.plus.2)
biocLite("ggplot2")
library(ggplot2)
biocLite("reshape")
biocLite("plot3D")
biocLite("gplots")
biocLite("ggdendro")
biocLite("RColorBrewer")
biocLite("frma")

library("reshape")
library("plot3D")
library("gplots")
library(ggdendro)
library(RColorBrewer)
library(RColorBrewer)
library(frma)
biocLite("hgu133plus2frmavecs")
#R CMD INSTALL -l R/x86_64-pc-linux-gnu-library/3.2/ /media/H_driver/Aimin_project/hgu133plus2frmavecs_1.5.0.tar.gz
library(hgu133plus2frmavecs)
biocLite("Homo.sapiens")
library("Homo.sapiens")
library(hgu133plus2.db)

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

ReformatCelName<-function(data.set2.cancer.GSM.cel.mapping.3){
data.set2.cancer.GSM.cel.mapping.3.2<-gsub(".cel","",data.set2.cancer.GSM.cel.mapping.3)
data.set2.cancer.GSM.cel.mapping.3.2.2<-gsub(".CEL","",as.character(data.set2.cancer.GSM.cel.mapping.3.2))
data.set2.cancer.GSM.cel.mapping.3.2.2.2<-gsub("-","_",as.character(data.set2.cancer.GSM.cel.mapping.3.2.2))
data.set2.cancer.GSM.cel.mapping.4<-cbind(data.set2.cancer.GSM.cel.mapping.3,data.set2.cancer.GSM.cel.mapping.3.2.2.2,sapply(strsplit(data.set2.cancer.GSM.cel.mapping.3.2.2.2,"_"),"[[",1))
return(data.set2.cancer.GSM.cel.mapping.4)
}

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

#Use frma method to do normalization
setwd("/media/H_driver/2015/Sophia/Cel_file_frma/")
data.set123.raw <- ReadAffy()
data.set123.frma<-frma(data.set123.raw)
ed.set123.frma <- exprs(data.set123.frma)
samp.frma <- sampleNames(data.set123.frma)
probes.frma <- featureNames(data.set123.frma)

samp.frma.1<-gsub(" - ","_",samp.frma)
samp.frma.2<-gsub(".CEL","",samp.frma.1)
samp.frma.3<-gsub(".cel","",samp.frma.2)
samp.frma.4<-gsub("-","_",samp.frma.3)
samp.frma.5<-cbind(samp.frma,samp.frma.4)
samp.frma.cancer.index<-which(samp.frma.5[,2] %in% cel.file.name.key.set.symbol.st134[,1])
samp.frma.6<-samp.frma.5[samp.frma.cancer.index,]

data.set123.frma.cancer<-data.set123.frma[,samp.frma.cancer.index]
ed.set123.frma.cancer<-ed.set123.frma[,samp.frma.cancer.index]
dim(ed.set123.frma.cancer)

samp.frma.7<-cbind(samp.frma.6,colnames(ed.set123.frma.cancer))
samp.frma.8<-samp.frma.7[,2:3]
colnames(samp.frma.8)<-c("key","filename_frma")

cel.file.sample.infor.no.8<-merge(cel.file.name.key.set.symbol.st134,samp.frma.8,by="key",sort=FALSE)

GeneSym.all <- getSYMBOL(rownames(ed.set123.frma.cancer), "hgu133plus2.db")


#<-cbind(names(GeneSym.all),GeneSym.all)

# #Collapse probes into genes
ndata<-ed.set123.frma.cancer
geneSymbol=GeneSym.all
tempdata.byGSym =data.frame(ndata, Symbol = geneSymbol)

colnames(tempdata.byGSym)

tempdata.byGSym.2<-tempdata.byGSym[-which(is.na(tempdata.byGSym$Symbol)),]

test.sample<-ed.set123.frma.cancer[which(GeneSym.all=="ERBB4"),which(colnames(ed.set123.frma.cancer) %in% cel.file.sample.infor.no.8[which(cel.file.sample.infor.no.8[,5]=="st4"),6])]

test.sample.2<-data.frame(test.sample,Symbol="ERBB4")
test.sample.3 = test.sample.2
test.sample.3[,1:5] = apply(test.sample.3[,1:5],2,as.numeric)


data.byGSym = ddply(tempdata.byGSym.2, c("Symbol"),function(h)
  summary = apply(h[,1:(dim(tempdata.byGSym.2)[2]-1)],2,max)
)


rownames(data.byGSym)=data.byGSym$Symbol
colnames(data.byGSym)

data.byGSym.2<-data.byGSym[,-1]
temp = apply(data.byGSym.2,2,as.numeric)
rownames(temp) = rownames(data.byGSym.2)
colnames(temp) = colnames(ndata)

data.byGSym.2 = temp
colnames(data.byGSym.2)
#data.byGSym.2 = ed.set123.frma.cancer

SampleType<-cel.file.sample.infor.no.8$subtype
heatmap_wPCA(data.byGSym.2, "heatmap_allsample_frma_all_probes.pdf","PCA_allsample_frma_all_probes.pdf","/media/H_driver/2015/Sophia/Results/", SampleType)

f.st134.frma <- factor(cel.file.sample.infor.no.8$subtype)
design.st134.frma <- model.matrix(~0+f.st134.frma)
colnames(design.st134.frma) <- levels(f.st134.frma)

reorder.index<-match(cel.file.sample.infor.no.8[,6],colnames(data.byGSym.2))

data.byGSym.3<-data.byGSym.2[,reorder.index]

#colnames(data.byGSym.3)

fit.st134.gene.frma <- lmFit(data.byGSym.3, design.st134.frma)
cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134.frma)
fit2.st134.frma  <- contrasts.fit(fit.st134.gene.frma, cont.matrix.st134.frma)
fit2.st134.frma  <- eBayes(fit2.st134.frma)

TopTableSt34.gene.frma<-topTable(fit2.st134.frma,coef=1,n=dim(data.byGSym.2)[1])
TopTableSt41.gene.frma<-topTable(fit2.st134.frma,coef=2,n=dim(data.byGSym.2)[1])

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

OutPut2HtmlTable(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/Results/","St41-8-frma.html","st4-vs-st1_HTML_report")
OutPut2HtmlTable(TopTableSt34.gene.frma,"/media/H_driver/2015/Sophia/Results/","St34-8-frma.html","st3-vs-st4_HTML_report")

#Generate the rank file for GSEA
GenerateRankFile4GSEA<-function(TopTableSt41.gene.frma,out_dir,out_file_name,up_down_top){

  FC.sign=sign(TopTableSt41.gene.frma[,1])
  p.inve=1/TopTableSt41.gene.frma[,4]
  Gene.Rank.Score=FC.sign*p.inve

  TopTableSt41.gene.frma.GSEA<-cbind(TopTableSt41.gene.frma,Gene.Rank.Score)

  TopTableSt41.gene.frma.GSEA.sorted.by.score<-TopTableSt41.gene.frma.GSEA[order(TopTableSt41.gene.frma.GSEA[,7],decreasing = TRUE),]

  Re.GSEA<-cbind(rownames(TopTableSt41.gene.frma.GSEA.sorted.by.score),TopTableSt41.gene.frma.GSEA.sorted.by.score[,7])
  colnames(Re.GSEA)=c("GeneName","Score")

  Re.GSEA.2<-rbind(head(Re.GSEA,up_down_top),tail(Re.GSEA,up_down_top))

  write.table(Re.GSEA.2, file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = FALSE)

  return(TopTableSt41.gene.frma.GSEA.sorted.by.score)
}

GenerateRankFile4GSEA(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/Results/","St41-4_GSEA_up_down_300.rnk",300)

#write.csv(TopTableSt41.gene.frma, file="/media/H_driver/2015/Sophia/Results/St41.csv",quote =FALSE)

#Use all data to perform GSEA analysis
GenerateFiles4GSEA<-function(ed.set123.frma.cancer,out_dir,out_file_name,subtypeA,subtypeB){

  index.st4<-which(colnames(ed.set123.frma.cancer) %in% cel.file.name.key.set.symbol.st134[which(cel.file.name.key.set.symbol.st134[,5]==subtypeA),2])
  index.st1<-which(colnames(ed.set123.frma.cancer) %in% cel.file.name.key.set.symbol.st134[which(cel.file.name.key.set.symbol.st134[,5]==subtypeB),2])

  ed.set123.frma.cancer.st41<-cbind(ed.set123.frma.cancer[,index.st4],ed.set123.frma.cancer[,index.st1])
  cat(subtypeA,"\t",length(index.st4),subtypeB,"\t",length(index.st1),"\t",dim(ed.set123.frma.cancer.st41),"\n")

  ed.set123.frma.cancer.st41.reformat<-cbind(rownames(ed.set123.frma.cancer.st41),rep("NA",dim(ed.set123.frma.cancer.st41)[1]),ed.set123.frma.cancer.st41)
  colnames(ed.set123.frma.cancer.st41.reformat)[c(1,2)]=c("NAME","Description")

  write.table(ed.set123.frma.cancer.st41.reformat,file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = TRUE)

}

GenerateFiles4GSEA(ed.set123.frma.cancer,"/media/H_driver/2015/Sophia/Results/","St41-4_GSEA.xls","st4","st1")
GenerateFiles4GSEA(ed.set123.frma.cancer,"/media/H_driver/2015/Sophia/Results/","St34-4_GSEA.xls","st3","st4")

save.image(file="/media/H_driver/Aimin_project/Sophia/Data_set_3_13_2016.RData")
savehistory(file="/media/H_driver/Aimin_project/Sophia/Data_set_3_13_2016.Rhistory")#load the affy library
.libPaths()
library(gdata)
library(affy)
library(affyio)
library(annotate)
library(plyr)

#biocLite("mogene20sttranscriptcluster.db")
#library(mogene20stprobeset.db)
#library(mogene20sttranscriptcluster.db)

#biocLite("hugene10sttranscriptcluster.db")
#library(hugene10sttranscriptcluster.db)

biocLite("pd.hg.u133.plus.2")
library(pd.hg.u133.plus.2)
biocLite("ggplot2")
library(ggplot2)
biocLite("reshape")
biocLite("plot3D")
biocLite("gplots")
biocLite("ggdendro")
biocLite("RColorBrewer")
biocLite("frma")

library("reshape")
library("plot3D")
library("gplots")
library(ggdendro)
library(RColorBrewer)
library(RColorBrewer)
library(frma)
biocLite("hgu133plus2frmavecs")
#R CMD INSTALL -l R/x86_64-pc-linux-gnu-library/3.2/ /media/H_driver/Aimin_project/hgu133plus2frmavecs_1.5.0.tar.gz
library(hgu133plus2frmavecs)
biocLite("Homo.sapiens")
library("Homo.sapiens")
library(hgu133plus2.db)

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

ReformatCelName<-function(data.set2.cancer.GSM.cel.mapping.3){
data.set2.cancer.GSM.cel.mapping.3.2<-gsub(".cel","",data.set2.cancer.GSM.cel.mapping.3)
data.set2.cancer.GSM.cel.mapping.3.2.2<-gsub(".CEL","",as.character(data.set2.cancer.GSM.cel.mapping.3.2))
data.set2.cancer.GSM.cel.mapping.3.2.2.2<-gsub("-","_",as.character(data.set2.cancer.GSM.cel.mapping.3.2.2))
data.set2.cancer.GSM.cel.mapping.4<-cbind(data.set2.cancer.GSM.cel.mapping.3,data.set2.cancer.GSM.cel.mapping.3.2.2.2,sapply(strsplit(data.set2.cancer.GSM.cel.mapping.3.2.2.2,"_"),"[[",1))
return(data.set2.cancer.GSM.cel.mapping.4)
}

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

#Use frma method to do normalization
setwd("/media/H_driver/2015/Sophia/Cel_file_frma/")
data.set123.raw <- ReadAffy()
data.set123.frma<-frma(data.set123.raw)
ed.set123.frma <- exprs(data.set123.frma)
samp.frma <- sampleNames(data.set123.frma)
probes.frma <- featureNames(data.set123.frma)

samp.frma.1<-gsub(" - ","_",samp.frma)
samp.frma.2<-gsub(".CEL","",samp.frma.1)
samp.frma.3<-gsub(".cel","",samp.frma.2)
samp.frma.4<-gsub("-","_",samp.frma.3)
samp.frma.5<-cbind(samp.frma,samp.frma.4)
samp.frma.cancer.index<-which(samp.frma.5[,2] %in% cel.file.name.key.set.symbol.st134[,1])
samp.frma.6<-samp.frma.5[samp.frma.cancer.index,]

data.set123.frma.cancer<-data.set123.frma[,samp.frma.cancer.index]
ed.set123.frma.cancer<-ed.set123.frma[,samp.frma.cancer.index]
dim(ed.set123.frma.cancer)

samp.frma.7<-cbind(samp.frma.6,colnames(ed.set123.frma.cancer))
samp.frma.8<-samp.frma.7[,2:3]
colnames(samp.frma.8)<-c("key","filename_frma")

cel.file.sample.infor.no.8<-merge(cel.file.name.key.set.symbol.st134,samp.frma.8,by="key",sort=FALSE)

GeneSym.all <- getSYMBOL(rownames(ed.set123.frma.cancer), "hgu133plus2.db")


#<-cbind(names(GeneSym.all),GeneSym.all)

# #Collapse probes into genes
ndata<-ed.set123.frma.cancer

geneSymbol=GeneSym.all
tempdata.byGSym =data.frame(ndata, Symbol = geneSymbol)

colnames(tempdata.byGSym)

tempdata.byGSym.2<-tempdata.byGSym[-which(is.na(tempdata.byGSym$Symbol)),]

# test.sample<-ed.set123.frma.cancer[which(GeneSym.all=="ERBB4"),which(colnames(ed.set123.frma.cancer) %in% cel.file.sample.infor.no.8[which(cel.file.sample.infor.no.8[,5]=="st4"),6])]
#
# test.sample.2<-data.frame(test.sample,Symbol="ERBB4")
# test.sample.3 = test.sample.2
# test.sample.3[,1:5] = apply(test.sample.3[,1:5],2,as.numeric)

data.byGSym = ddply(tempdata.byGSym.2, c("Symbol"),function(h)
  summary = apply(h[,1:(dim(tempdata.byGSym.2)[2]-1)],2,max)
)


rownames(data.byGSym)=data.byGSym$Symbol
colnames(data.byGSym)

data.byGSym.2<-data.byGSym[,-1]
temp = apply(data.byGSym.2,2,as.numeric)
rownames(temp) = rownames(data.byGSym.2)
colnames(temp) = colnames(ndata)

data.byGSym.2 = temp
colnames(data.byGSym.2)
#data.byGSym.2 = ed.set123.frma.cancer

SampleType<-cel.file.sample.infor.no.8$subtype
heatmap_wPCA(data.byGSym.2, "heatmap_allsample_frma_all_probes.pdf","PCA_allsample_frma_all_probes.pdf","/media/H_driver/2015/Sophia/Results/", SampleType)

f.st134.frma <- factor(cel.file.sample.infor.no.8$subtype)
design.st134.frma <- model.matrix(~0+f.st134.frma)
colnames(design.st134.frma) <- levels(f.st134.frma)

reorder.index<-match(cel.file.sample.infor.no.8[,6],colnames(data.byGSym.2))

data.byGSym.3<-data.byGSym.2[,reorder.index]

#colnames(data.byGSym.3)

fit.st134.gene.frma <- lmFit(data.byGSym.3, design.st134.frma)
cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134.frma)
fit2.st134.frma  <- contrasts.fit(fit.st134.gene.frma, cont.matrix.st134.frma)
fit2.st134.frma  <- eBayes(fit2.st134.frma)

TopTableSt34.gene.frma<-topTable(fit2.st134.frma,coef=1,n=dim(data.byGSym.2)[1])
TopTableSt41.gene.frma<-topTable(fit2.st134.frma,coef=2,n=dim(data.byGSym.2)[1])

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

OutPut2HtmlTable(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/Results/","St41-8-frma.html","st4-vs-st1_HTML_report")
OutPut2HtmlTable(TopTableSt34.gene.frma,"/media/H_driver/2015/Sophia/Results/","St34-8-frma.html","st3-vs-st4_HTML_report")

#Generate the rank file for GSEA
GenerateRankFile4GSEA<-function(TopTableSt41.gene.frma,out_dir,out_file_name,up_down_top){

  FC.sign=sign(TopTableSt41.gene.frma[,1])
  p.inve=1/TopTableSt41.gene.frma[,4]
  Gene.Rank.Score=FC.sign*p.inve

  TopTableSt41.gene.frma.GSEA<-cbind(TopTableSt41.gene.frma,Gene.Rank.Score)

  TopTableSt41.gene.frma.GSEA.sorted.by.score<-TopTableSt41.gene.frma.GSEA[order(TopTableSt41.gene.frma.GSEA[,7],decreasing = TRUE),]

  Re.GSEA<-cbind(rownames(TopTableSt41.gene.frma.GSEA.sorted.by.score),TopTableSt41.gene.frma.GSEA.sorted.by.score[,7])
  colnames(Re.GSEA)=c("GeneName","Score")

  Re.GSEA.2<-rbind(head(Re.GSEA,up_down_top),tail(Re.GSEA,up_down_top))

  write.table(Re.GSEA.2, file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = FALSE)

  return(TopTableSt41.gene.frma.GSEA.sorted.by.score)
}

GenerateRankFile4GSEA(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/Results/","St41-4_GSEA_up_down_300.rnk",300)

#write.csv(TopTableSt41.gene.frma, file="/media/H_driver/2015/Sophia/Results/St41.csv",quote =FALSE)


##Use all probes for DE
reorder.index.all.probes<-match(cel.file.sample.infor.no.8[,6],colnames(ed.set123.frma.cancer))
ed.set123.frma.cancer.reorder<-ed.set123.frma.cancer[,reorder.index.all.probes]

check.match<-cbind(colnames(ed.set123.frma.cancer.reorder),design.st134.frma,cel.file.sample.infor.no.8)

fit.st134.all.probes.frma <- lmFit(ed.set123.frma.cancer.reorder, design.st134.frma)
cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134.frma)
fit2.st134.all.probes.frma  <- contrasts.fit(fit.st134.all.probes.frma, cont.matrix.st134.frma)
fit2.st134.all.probes.frma  <- eBayes(fit2.st134.all.probes.frma)

TopTableSt34.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=1,n=dim(ed.set123.frma.cancer.reorder)[1])
TopTableSt41.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=2,n=dim(ed.set123.frma.cancer.reorder)[1])

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

RE41<-OutPut2HtmlTableAllProbes(TopTableSt41.all.probes.frma,"/media/H_driver/2015/Sophia/Results/","St41-all-probes-frma-3.html","st4-vs-st1_HTML_report")
RE34<-OutPut2HtmlTableAllProbes(TopTableSt34.all.probes.frma,"/media/H_driver/2015/Sophia/Results/","St34-all-probes-frma-3.html","st3-vs-st4_HTML_report")

#head(RE41)

GenerateRank4AllProbes<-function(TopTableSt41.all.probes.frma,up_down_top){

  FC.sign=sign(TopTableSt41.all.probes.frma[,1])
  p.inve=1/TopTableSt41.all.probes.frma[,4]
  Gene.Rank.Score=FC.sign*p.inve

  TopTableSt41.gene.frma.GSEA<-cbind(TopTableSt41.all.probes.frma,Gene.Rank.Score)

  TopTableSt41.gene.frma.GSEA.sorted.by.score<-TopTableSt41.gene.frma.GSEA[order(TopTableSt41.gene.frma.GSEA[,7],decreasing = TRUE),]

  Re.GSEA<-cbind(rownames(TopTableSt41.gene.frma.GSEA.sorted.by.score),TopTableSt41.gene.frma.GSEA.sorted.by.score[,7])
  colnames(Re.GSEA)=c("Probe","Score")

  Re.GSEA.2<-rbind(head(Re.GSEA,up_down_top),tail(Re.GSEA,up_down_top))

  #write.table(Re.GSEA.2, file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = FALSE)

  return(Re.GSEA.2)
}


TopTableSt41.most.DE.probes.frma<-GenerateRank4AllProbes(TopTableSt41.all.probes.frma,200)
TopTableSt34.most.DE.probes.frma<-GenerateRank4AllProbes(TopTableSt34.all.probes.frma,200)

ed.set123.frma.cancer.reorder.st1<-ed.set123.frma.cancer.reorder[,which(colnames(ed.set123.frma.cancer.reorder) %in% check.match[which(check.match$st1==1),1])]
ed.set123.frma.cancer.reorder.st2<-ed.set123.frma.cancer.reorder[,which(colnames(ed.set123.frma.cancer.reorder) %in% check.match[which(check.match$st2==1),1])]
ed.set123.frma.cancer.reorder.st3<-ed.set123.frma.cancer.reorder[,which(colnames(ed.set123.frma.cancer.reorder) %in% check.match[which(check.match$st3==1),1])]
ed.set123.frma.cancer.reorder.st4<-ed.set123.frma.cancer.reorder[,which(colnames(ed.set123.frma.cancer.reorder) %in% check.match[which(check.match$st4==1),1])]

dim(ed.set123.frma.cancer.reorder.st1)
dim(ed.set123.frma.cancer.reorder.st2)
dim(ed.set123.frma.cancer.reorder.st3)
dim(ed.set123.frma.cancer.reorder.st4)

ed.set123.frma.cancer.reorder.st41<-cbind(ed.set123.frma.cancer.reorder.st4,ed.set123.frma.cancer.reorder.st1)
ed.set123.frma.cancer.reorder.st34<-cbind(ed.set123.frma.cancer.reorder.st3,ed.set123.frma.cancer.reorder.st4)

dim(ed.set123.frma.cancer.reorder.st41)
dim(ed.set123.frma.cancer.reorder.st34)

ed.set123.frma.cancer.reorder.st41.MDE.100<-ed.set123.frma.cancer.reorder.st41[which(rownames(ed.set123.frma.cancer.reorder.st41) %in%  TopTableSt41.most.DE.probes.frma[,1]),]
ed.set123.frma.cancer.reorder.st34.MDE.100<-ed.set123.frma.cancer.reorder.st34[which(rownames(ed.set123.frma.cancer.reorder.st34) %in%  TopTableSt34.most.DE.probes.frma[,1]),]

subtype.s41<-c(as.character(check.match[which(check.match[,9]=="st4"),9]),as.character(check.match[which(check.match[,9]=="st1"),9]))
subtype.s34<-c(as.character(check.match[which(check.match[,9]=="st3"),9]),as.character(check.match[which(check.match[,9]=="st4"),9]))

Draw_heatmap(ed.set123.frma.cancer.reorder.st41.MDE.100,"heatmap_allsample_frma_MDE_200_probes_41.pdf","/media/H_driver/2015/Sophia/Results/",subtype.s41)
Draw_heatmap(ed.set123.frma.cancer.reorder.st34.MDE.100,"heatmap_allsample_frma_MDE_200_probes_34.pdf","/media/H_driver/2015/Sophia/Results/",subtype.s34)


Draw_PCA(ed.set123.frma.cancer.reorder,"PCA_58_cancer_sample_frma_all_probes.pdf","/media/H_driver/2015/Sophia/Results/",check.match$subtype)

save(heatmap_wPCA,ed.set123.frma.cancer.reorder,check.match,file="Data_4_heatmap_PCA.RData")

#Use all data to perform GSEA analysis
GenerateFiles4GSEA<-function(ed.set123.frma.cancer,out_dir,out_file_name,subtypeA,subtypeB){

  index.st4<-which(colnames(ed.set123.frma.cancer) %in% cel.file.name.key.set.symbol.st134[which(cel.file.name.key.set.symbol.st134[,5]==subtypeA),2])
  index.st1<-which(colnames(ed.set123.frma.cancer) %in% cel.file.name.key.set.symbol.st134[which(cel.file.name.key.set.symbol.st134[,5]==subtypeB),2])

  ed.set123.frma.cancer.st41<-cbind(ed.set123.frma.cancer[,index.st4],ed.set123.frma.cancer[,index.st1])
  cat(subtypeA,"\t",length(index.st4),subtypeB,"\t",length(index.st1),"\t",dim(ed.set123.frma.cancer.st41),"\n")

  ed.set123.frma.cancer.st41.reformat<-cbind(rownames(ed.set123.frma.cancer.st41),rep("NA",dim(ed.set123.frma.cancer.st41)[1]),ed.set123.frma.cancer.st41)
  colnames(ed.set123.frma.cancer.st41.reformat)[c(1,2)]=c("NAME","Description")

  write.table(ed.set123.frma.cancer.st41.reformat,file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = TRUE)

}

GenerateFiles4GSEA(ed.set123.frma.cancer,"/media/H_driver/2015/Sophia/Results/","St41-4_GSEA.xls","st4","st1")
GenerateFiles4GSEA(ed.set123.frma.cancer,"/media/H_driver/2015/Sophia/Results/","St34-4_GSEA.xls","st3","st4")
GenerateFiles4GSEA(ed.set123.frma.cancer,"/media/H_driver/2015/Sophia/Results/","St13-4_GSEA.xls","st1","st3")

save.image(file="/media/H_driver/Aimin_project/Sophia/Data_set_3_16_2016_all_probes_based_New.RData")
savehistory(file="/media/H_driver/Aimin_project/Sophia/Data_set_3_16_2016_all_probes_based_New.Rhistory")
#load the affy library
.libPaths()
library(gdata)
library(affy)
library(affyio)
library(annotate)
library(plyr)

#biocLite("mogene20sttranscriptcluster.db")
#library(mogene20stprobeset.db)
#library(mogene20sttranscriptcluster.db)

#biocLite("hugene10sttranscriptcluster.db")
#library(hugene10sttranscriptcluster.db)

biocLite("pd.hg.u133.plus.2")
library(pd.hg.u133.plus.2)
biocLite("ggplot2")
library(ggplot2)
biocLite("reshape")
biocLite("plot3D")
biocLite("gplots")
biocLite("ggdendro")
biocLite("RColorBrewer")
biocLite("frma")

library("reshape")
library("plot3D")
library("gplots")
library(ggdendro)
library(RColorBrewer)
library(RColorBrewer)
library(frma)
biocLite("hgu133plus2frmavecs")
#R CMD INSTALL -l R/x86_64-pc-linux-gnu-library/3.2/ /media/H_driver/Aimin_project/hgu133plus2frmavecs_1.5.0.tar.gz
library(hgu133plus2frmavecs)
biocLite("Homo.sapiens")
library("Homo.sapiens")
library(hgu133plus2.db)

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

ReformatCelName<-function(data.set2.cancer.GSM.cel.mapping.3){
data.set2.cancer.GSM.cel.mapping.3.2<-gsub(".cel","",data.set2.cancer.GSM.cel.mapping.3)
data.set2.cancer.GSM.cel.mapping.3.2.2<-gsub(".CEL","",as.character(data.set2.cancer.GSM.cel.mapping.3.2))
data.set2.cancer.GSM.cel.mapping.3.2.2.2<-gsub("-","_",as.character(data.set2.cancer.GSM.cel.mapping.3.2.2))
data.set2.cancer.GSM.cel.mapping.4<-cbind(data.set2.cancer.GSM.cel.mapping.3,data.set2.cancer.GSM.cel.mapping.3.2.2.2,sapply(strsplit(data.set2.cancer.GSM.cel.mapping.3.2.2.2,"_"),"[[",1))
return(data.set2.cancer.GSM.cel.mapping.4)
}

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

#Use frma method to do normalization
setwd("/media/H_driver/2015/Sophia/Cel_file_frma/")
data.set123.raw <- ReadAffy()
data.set123.frma<-frma(data.set123.raw)
ed.set123.frma <- exprs(data.set123.frma)
samp.frma <- sampleNames(data.set123.frma)
probes.frma <- featureNames(data.set123.frma)

samp.frma.1<-gsub(" - ","_",samp.frma)
samp.frma.2<-gsub(".CEL","",samp.frma.1)
samp.frma.3<-gsub(".cel","",samp.frma.2)
samp.frma.4<-gsub("-","_",samp.frma.3)
samp.frma.5<-cbind(samp.frma,samp.frma.4)
samp.frma.cancer.index<-which(samp.frma.5[,2] %in% cel.file.name.key.set.symbol.st134[,1])
samp.frma.6<-samp.frma.5[samp.frma.cancer.index,]

data.set123.frma.cancer<-data.set123.frma[,samp.frma.cancer.index]
ed.set123.frma.cancer<-ed.set123.frma[,samp.frma.cancer.index]
dim(ed.set123.frma.cancer)

samp.frma.7<-cbind(samp.frma.6,colnames(ed.set123.frma.cancer))
samp.frma.8<-samp.frma.7[,2:3]
colnames(samp.frma.8)<-c("key","filename_frma")

cel.file.sample.infor.no.8<-merge(cel.file.name.key.set.symbol.st134,samp.frma.8,by="key",sort=FALSE)

GeneSym.all <- getSYMBOL(rownames(ed.set123.frma.cancer), "hgu133plus2.db")


#<-cbind(names(GeneSym.all),GeneSym.all)

# #Collapse probes into genes
ndata<-ed.set123.frma.cancer

geneSymbol=GeneSym.all
tempdata.byGSym =data.frame(ndata, Symbol = geneSymbol)

colnames(tempdata.byGSym)

tempdata.byGSym.2<-tempdata.byGSym[-which(is.na(tempdata.byGSym$Symbol)),]

# test.sample<-ed.set123.frma.cancer[which(GeneSym.all=="ERBB4"),which(colnames(ed.set123.frma.cancer) %in% cel.file.sample.infor.no.8[which(cel.file.sample.infor.no.8[,5]=="st4"),6])]
#
# test.sample.2<-data.frame(test.sample,Symbol="ERBB4")
# test.sample.3 = test.sample.2
# test.sample.3[,1:5] = apply(test.sample.3[,1:5],2,as.numeric)

data.byGSym = ddply(tempdata.byGSym.2, c("Symbol"),function(h)
  summary = apply(h[,1:(dim(tempdata.byGSym.2)[2]-1)],2,max)
)


rownames(data.byGSym)=data.byGSym$Symbol
colnames(data.byGSym)

data.byGSym.2<-data.byGSym[,-1]
temp = apply(data.byGSym.2,2,as.numeric)
rownames(temp) = rownames(data.byGSym.2)
colnames(temp) = colnames(ndata)

data.byGSym.2 = temp
colnames(data.byGSym.2)
#data.byGSym.2 = ed.set123.frma.cancer

SampleType<-cel.file.sample.infor.no.8$subtype
heatmap_wPCA(data.byGSym.2, "heatmap_allsample_frma_all_probes.pdf","PCA_allsample_frma_all_probes.pdf","/media/H_driver/2015/Sophia/Results/", SampleType)

f.st134.frma <- factor(cel.file.sample.infor.no.8$subtype)
design.st134.frma <- model.matrix(~0+f.st134.frma)
colnames(design.st134.frma) <- levels(f.st134.frma)

reorder.index<-match(cel.file.sample.infor.no.8[,6],colnames(data.byGSym.2))

data.byGSym.3<-data.byGSym.2[,reorder.index]

#colnames(data.byGSym.3)

fit.st134.gene.frma <- lmFit(data.byGSym.3, design.st134.frma)
cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134.frma)
fit2.st134.frma  <- contrasts.fit(fit.st134.gene.frma, cont.matrix.st134.frma)
fit2.st134.frma  <- eBayes(fit2.st134.frma)

TopTableSt34.gene.frma<-topTable(fit2.st134.frma,coef=1,n=dim(data.byGSym.2)[1])
TopTableSt41.gene.frma<-topTable(fit2.st134.frma,coef=2,n=dim(data.byGSym.2)[1])

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

OutPut2HtmlTable(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/Results/","St41-8-frma.html","st4-vs-st1_HTML_report")
OutPut2HtmlTable(TopTableSt34.gene.frma,"/media/H_driver/2015/Sophia/Results/","St34-8-frma.html","st3-vs-st4_HTML_report")

#Generate the rank file for GSEA
GenerateRankFile4GSEA<-function(TopTableSt41.gene.frma,out_dir,out_file_name,up_down_top){

  FC.sign=sign(TopTableSt41.gene.frma[,1])
  p.inve=1/TopTableSt41.gene.frma[,4]
  Gene.Rank.Score=FC.sign*p.inve

  TopTableSt41.gene.frma.GSEA<-cbind(TopTableSt41.gene.frma,Gene.Rank.Score)

  TopTableSt41.gene.frma.GSEA.sorted.by.score<-TopTableSt41.gene.frma.GSEA[order(TopTableSt41.gene.frma.GSEA[,7],decreasing = TRUE),]

  Re.GSEA<-cbind(rownames(TopTableSt41.gene.frma.GSEA.sorted.by.score),TopTableSt41.gene.frma.GSEA.sorted.by.score[,7])
  colnames(Re.GSEA)=c("GeneName","Score")

  Re.GSEA.2<-rbind(head(Re.GSEA,up_down_top),tail(Re.GSEA,up_down_top))

  write.table(Re.GSEA.2, file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = FALSE)

  return(TopTableSt41.gene.frma.GSEA.sorted.by.score)
}

GenerateRankFile4GSEA(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/Results/","St41-4_GSEA_up_down_300.rnk",300)

#write.csv(TopTableSt41.gene.frma, file="/media/H_driver/2015/Sophia/Results/St41.csv",quote =FALSE)


##Use all probes for DE
reorder.index.all.probes<-match(cel.file.sample.infor.no.8[,6],colnames(ed.set123.frma.cancer))
ed.set123.frma.cancer.reorder<-ed.set123.frma.cancer[,reorder.index.all.probes]

check.match<-cbind(colnames(ed.set123.frma.cancer.reorder),design.st134.frma,cel.file.sample.infor.no.8)

fit.st134.all.probes.frma <- lmFit(ed.set123.frma.cancer.reorder, design.st134.frma)
cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134.frma)
fit2.st134.all.probes.frma  <- contrasts.fit(fit.st134.all.probes.frma, cont.matrix.st134.frma)
fit2.st134.all.probes.frma  <- eBayes(fit2.st134.all.probes.frma)

TopTableSt34.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=1,n=dim(ed.set123.frma.cancer.reorder)[1])
TopTableSt41.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=2,n=dim(ed.set123.frma.cancer.reorder)[1])

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

RE41<-OutPut2HtmlTableAllProbes(TopTableSt41.all.probes.frma,"/media/H_driver/2015/Sophia/Results/","St41-all-probes-frma-3.html","st4-vs-st1_HTML_report")
RE34<-OutPut2HtmlTableAllProbes(TopTableSt34.all.probes.frma,"/media/H_driver/2015/Sophia/Results/","St34-all-probes-frma-3.html","st3-vs-st4_HTML_report")

heatmap_wPCA(ed.set123.frma.cancer.reorder,"heatmap_allsample_frma_all_probes_2.pdf","PCA_allsample_frma_all_probes_.pdf","/media/H_driver/2015/Sophia/Results/",check.match$subtype)

save(heatmap_wPCA,ed.set123.frma.cancer.reorder,check.match,file="Data_4_heatmap_PCA.RData")

#Use all data to perform GSEA analysis
GenerateFiles4GSEA<-function(ed.set123.frma.cancer,out_dir,out_file_name,subtypeA,subtypeB){

  index.st4<-which(colnames(ed.set123.frma.cancer) %in% cel.file.name.key.set.symbol.st134[which(cel.file.name.key.set.symbol.st134[,5]==subtypeA),2])
  index.st1<-which(colnames(ed.set123.frma.cancer) %in% cel.file.name.key.set.symbol.st134[which(cel.file.name.key.set.symbol.st134[,5]==subtypeB),2])

  ed.set123.frma.cancer.st41<-cbind(ed.set123.frma.cancer[,index.st4],ed.set123.frma.cancer[,index.st1])
  cat(subtypeA,"\t",length(index.st4),subtypeB,"\t",length(index.st1),"\t",dim(ed.set123.frma.cancer.st41),"\n")

  ed.set123.frma.cancer.st41.reformat<-cbind(rownames(ed.set123.frma.cancer.st41),rep("NA",dim(ed.set123.frma.cancer.st41)[1]),ed.set123.frma.cancer.st41)
  colnames(ed.set123.frma.cancer.st41.reformat)[c(1,2)]=c("NAME","Description")

  write.table(ed.set123.frma.cancer.st41.reformat,file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = TRUE)

}

GenerateFiles4GSEA(ed.set123.frma.cancer,"/media/H_driver/2015/Sophia/Results/","St41-4_GSEA.xls","st4","st1")
GenerateFiles4GSEA(ed.set123.frma.cancer,"/media/H_driver/2015/Sophia/Results/","St34-4_GSEA.xls","st3","st4")

save.image(file="/media/H_driver/Aimin_project/Sophia/Data_set_3_16_2016_all_probes_based_2.RData")
savehistory(file="/media/H_driver/Aimin_project/Sophia/Data_set_3_16_2016_all_probes_based_2.Rhistory")
#load the affy library
.libPaths()
library(gdata)
library(affy)
library(affyio)
library(annotate)
library(plyr)

#biocLite("mogene20sttranscriptcluster.db")
#library(mogene20stprobeset.db)
#library(mogene20sttranscriptcluster.db)

#biocLite("hugene10sttranscriptcluster.db")
#library(hugene10sttranscriptcluster.db)

biocLite("pd.hg.u133.plus.2")
library(pd.hg.u133.plus.2)
biocLite("ggplot2")
library(ggplot2)
biocLite("reshape")
biocLite("plot3D")
biocLite("gplots")
biocLite("ggdendro")
biocLite("RColorBrewer")
biocLite("frma")

library("reshape")
library("plot3D")
library("gplots")
library(ggdendro)
library(RColorBrewer)
library(RColorBrewer)
library(frma)
biocLite("hgu133plus2frmavecs")
#R CMD INSTALL -l R/x86_64-pc-linux-gnu-library/3.2/ /media/H_driver/Aimin_project/hgu133plus2frmavecs_1.5.0.tar.gz
library(hgu133plus2frmavecs)
biocLite("Homo.sapiens")
library("Homo.sapiens")
library(hgu133plus2.db)

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

data.ER.PR.sample.info.new<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx".sheet=2)

data.ER.PR.sample.info.sheet1<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx",sheet=1)
data.ER.PR.sample.info.sheet2<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx",sheet=2)

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

ReformatCelName<-function(data.set2.cancer.GSM.cel.mapping.3){
data.set2.cancer.GSM.cel.mapping.3.2<-gsub(".cel","",data.set2.cancer.GSM.cel.mapping.3)
data.set2.cancer.GSM.cel.mapping.3.2.2<-gsub(".CEL","",as.character(data.set2.cancer.GSM.cel.mapping.3.2))
data.set2.cancer.GSM.cel.mapping.3.2.2.2<-gsub("-","_",as.character(data.set2.cancer.GSM.cel.mapping.3.2.2))
data.set2.cancer.GSM.cel.mapping.4<-cbind(data.set2.cancer.GSM.cel.mapping.3,data.set2.cancer.GSM.cel.mapping.3.2.2.2,sapply(strsplit(data.set2.cancer.GSM.cel.mapping.3.2.2.2,"_"),"[[",1))
return(data.set2.cancer.GSM.cel.mapping.4)
}

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

#Use frma method to do normalization
setwd("/media/H_driver/2015/Sophia/Cel_file_frma/")
data.set123.raw <- ReadAffy()
data.set123.frma<-frma(data.set123.raw)
ed.set123.frma <- exprs(data.set123.frma)
samp.frma <- sampleNames(data.set123.frma)
probes.frma <- featureNames(data.set123.frma)

samp.frma.1<-gsub(" - ","_",samp.frma)
samp.frma.2<-gsub(".CEL","",samp.frma.1)
samp.frma.3<-gsub(".cel","",samp.frma.2)
samp.frma.4<-gsub("-","_",samp.frma.3)
samp.frma.5<-cbind(samp.frma,samp.frma.4)
samp.frma.cancer.index<-which(samp.frma.5[,2] %in% cel.file.name.key.set.symbol.st134[,1])
samp.frma.6<-samp.frma.5[samp.frma.cancer.index,]

data.set123.frma.cancer<-data.set123.frma[,samp.frma.cancer.index]
ed.set123.frma.cancer<-ed.set123.frma[,samp.frma.cancer.index]
dim(ed.set123.frma.cancer)

samp.frma.7<-cbind(samp.frma.6,colnames(ed.set123.frma.cancer))
samp.frma.8<-samp.frma.7[,2:3]
colnames(samp.frma.8)<-c("key","filename_frma")

cel.file.sample.infor.no.8<-merge(cel.file.name.key.set.symbol.st134,samp.frma.8,by="key",sort=FALSE)

GeneSym.all <- getSYMBOL(rownames(ed.set123.frma.cancer), "hgu133plus2.db")


#<-cbind(names(GeneSym.all),GeneSym.all)

# #Collapse probes into genes
ndata<-ed.set123.frma.cancer

geneSymbol=GeneSym.all
tempdata.byGSym =data.frame(ndata, Symbol = geneSymbol)

colnames(tempdata.byGSym)

tempdata.byGSym.2<-tempdata.byGSym[-which(is.na(tempdata.byGSym$Symbol)),]

# test.sample<-ed.set123.frma.cancer[which(GeneSym.all=="ERBB4"),which(colnames(ed.set123.frma.cancer) %in% cel.file.sample.infor.no.8[which(cel.file.sample.infor.no.8[,5]=="st4"),6])]
#
# test.sample.2<-data.frame(test.sample,Symbol="ERBB4")
# test.sample.3 = test.sample.2
# test.sample.3[,1:5] = apply(test.sample.3[,1:5],2,as.numeric)

data.byGSym = ddply(tempdata.byGSym.2, c("Symbol"),function(h)
  summary = apply(h[,1:(dim(tempdata.byGSym.2)[2]-1)],2,max)
)


rownames(data.byGSym)=data.byGSym$Symbol
colnames(data.byGSym)

data.byGSym.2<-data.byGSym[,-1]
temp = apply(data.byGSym.2,2,as.numeric)
rownames(temp) = rownames(data.byGSym.2)
colnames(temp) = colnames(ndata)

data.byGSym.2 = temp
colnames(data.byGSym.2)
#data.byGSym.2 = ed.set123.frma.cancer

SampleType<-cel.file.sample.infor.no.8$subtype
heatmap_wPCA(data.byGSym.2, "heatmap_allsample_frma_all_probes.pdf","PCA_allsample_frma_all_probes.pdf","/media/H_driver/2015/Sophia/Results/", SampleType)

f.st134.frma <- factor(cel.file.sample.infor.no.8$subtype)
design.st134.frma <- model.matrix(~0+f.st134.frma)
colnames(design.st134.frma) <- levels(f.st134.frma)

reorder.index<-match(cel.file.sample.infor.no.8[,6],colnames(data.byGSym.2))

data.byGSym.3<-data.byGSym.2[,reorder.index]

#colnames(data.byGSym.3)

fit.st134.gene.frma <- lmFit(data.byGSym.3, design.st134.frma)
cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134.frma)
fit2.st134.frma  <- contrasts.fit(fit.st134.gene.frma, cont.matrix.st134.frma)
fit2.st134.frma  <- eBayes(fit2.st134.frma)

TopTableSt34.gene.frma<-topTable(fit2.st134.frma,coef=1,n=dim(data.byGSym.2)[1])
TopTableSt41.gene.frma<-topTable(fit2.st134.frma,coef=2,n=dim(data.byGSym.2)[1])

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

OutPut2HtmlTable(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/Results/","St41-8-frma.html","st4-vs-st1_HTML_report")
OutPut2HtmlTable(TopTableSt34.gene.frma,"/media/H_driver/2015/Sophia/Results/","St34-8-frma.html","st3-vs-st4_HTML_report")

#Generate the rank file for GSEA
GenerateRankFile4GSEA<-function(TopTableSt41.gene.frma,out_dir,out_file_name,up_down_top){

  FC.sign=sign(TopTableSt41.gene.frma[,1])
  p.inve=1/TopTableSt41.gene.frma[,4]
  Gene.Rank.Score=FC.sign*p.inve

  TopTableSt41.gene.frma.GSEA<-cbind(TopTableSt41.gene.frma,Gene.Rank.Score)

  TopTableSt41.gene.frma.GSEA.sorted.by.score<-TopTableSt41.gene.frma.GSEA[order(TopTableSt41.gene.frma.GSEA[,7],decreasing = TRUE),]

  Re.GSEA<-cbind(rownames(TopTableSt41.gene.frma.GSEA.sorted.by.score),TopTableSt41.gene.frma.GSEA.sorted.by.score[,7])
  colnames(Re.GSEA)=c("GeneName","Score")

  Re.GSEA.2<-rbind(head(Re.GSEA,up_down_top),tail(Re.GSEA,up_down_top))

  write.table(Re.GSEA.2, file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = FALSE)

  return(TopTableSt41.gene.frma.GSEA.sorted.by.score)
}

GenerateRankFile4GSEA(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/Results/","St41-4_GSEA_up_down_300.rnk",300)

#write.csv(TopTableSt41.gene.frma, file="/media/H_driver/2015/Sophia/Results/St41.csv",quote =FALSE)


##Use all probes for DE
reorder.index.all.probes<-match(cel.file.sample.infor.no.8[,6],colnames(ed.set123.frma.cancer))
ed.set123.frma.cancer.reorder<-ed.set123.frma.cancer[,reorder.index.all.probes]

check.match<-cbind(colnames(ed.set123.frma.cancer.reorder),design.st134.frma,cel.file.sample.infor.no.8)

fit.st134.all.probes.frma <- lmFit(ed.set123.frma.cancer.reorder, design.st134.frma)
#cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134.frma)
cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",st13="st1-st3",levels=design.st134.frma)

fit2.st134.all.probes.frma  <- contrasts.fit(fit.st134.all.probes.frma, cont.matrix.st134.frma)
fit2.st134.all.probes.frma  <- eBayes(fit2.st134.all.probes.frma)

TopTableSt34.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=1,n=dim(ed.set123.frma.cancer.reorder)[1])
TopTableSt41.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=2,n=dim(ed.set123.frma.cancer.reorder)[1])
TopTableSt13.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=3,n=dim(ed.set123.frma.cancer.reorder)[1])

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

RE41<-OutPut2HtmlTableAllProbes(TopTableSt41.all.probes.frma,"/media/H_driver/2015/Sophia/Results/","St41-all-probes-frma-3.html","st4-vs-st1_HTML_report")
RE34<-OutPut2HtmlTableAllProbes(TopTableSt34.all.probes.frma,"/media/H_driver/2015/Sophia/Results/","St34-all-probes-frma-3.html","st3-vs-st4_HTML_report")
RE13<-OutPut2HtmlTableAllProbes(TopTableSt13.all.probes.frma,"/media/H_driver/2015/Sophia/Results/","St13-all-probes-frma.html","st1-vs-st3_HTML_report")

#head(RE41)

GenerateRank4AllProbes<-function(TopTableSt41.all.probes.frma,up_down_top){

  FC.sign=sign(TopTableSt41.all.probes.frma[,1])
  p.inve=1/TopTableSt41.all.probes.frma[,4]
  Gene.Rank.Score=FC.sign*p.inve

  TopTableSt41.gene.frma.GSEA<-cbind(TopTableSt41.all.probes.frma,Gene.Rank.Score)

  TopTableSt41.gene.frma.GSEA.sorted.by.score<-TopTableSt41.gene.frma.GSEA[order(TopTableSt41.gene.frma.GSEA[,7],decreasing = TRUE),]

  Re.GSEA<-cbind(rownames(TopTableSt41.gene.frma.GSEA.sorted.by.score),TopTableSt41.gene.frma.GSEA.sorted.by.score[,7])
  colnames(Re.GSEA)=c("Probe","Score")

  Re.GSEA.2<-rbind(head(Re.GSEA,up_down_top),tail(Re.GSEA,up_down_top))

  #write.table(Re.GSEA.2, file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = FALSE)

  return(Re.GSEA.2)
}


TopTableSt41.most.DE.probes.frma<-GenerateRank4AllProbes(TopTableSt41.all.probes.frma,200)
TopTableSt34.most.DE.probes.frma<-GenerateRank4AllProbes(TopTableSt34.all.probes.frma,200)

ed.set123.frma.cancer.reorder.st1<-ed.set123.frma.cancer.reorder[,which(colnames(ed.set123.frma.cancer.reorder) %in% check.match[which(check.match$st1==1),1])]
ed.set123.frma.cancer.reorder.st2<-ed.set123.frma.cancer.reorder[,which(colnames(ed.set123.frma.cancer.reorder) %in% check.match[which(check.match$st2==1),1])]
ed.set123.frma.cancer.reorder.st3<-ed.set123.frma.cancer.reorder[,which(colnames(ed.set123.frma.cancer.reorder) %in% check.match[which(check.match$st3==1),1])]
ed.set123.frma.cancer.reorder.st4<-ed.set123.frma.cancer.reorder[,which(colnames(ed.set123.frma.cancer.reorder) %in% check.match[which(check.match$st4==1),1])]

dim(ed.set123.frma.cancer.reorder.st1)
dim(ed.set123.frma.cancer.reorder.st2)
dim(ed.set123.frma.cancer.reorder.st3)
dim(ed.set123.frma.cancer.reorder.st4)

ed.set123.frma.cancer.reorder.st41<-cbind(ed.set123.frma.cancer.reorder.st4,ed.set123.frma.cancer.reorder.st1)
ed.set123.frma.cancer.reorder.st34<-cbind(ed.set123.frma.cancer.reorder.st3,ed.set123.frma.cancer.reorder.st4)

dim(ed.set123.frma.cancer.reorder.st41)
dim(ed.set123.frma.cancer.reorder.st34)

ed.set123.frma.cancer.reorder.st41.MDE.100<-ed.set123.frma.cancer.reorder.st41[which(rownames(ed.set123.frma.cancer.reorder.st41) %in%  TopTableSt41.most.DE.probes.frma[,1]),]
ed.set123.frma.cancer.reorder.st34.MDE.100<-ed.set123.frma.cancer.reorder.st34[which(rownames(ed.set123.frma.cancer.reorder.st34) %in%  TopTableSt34.most.DE.probes.frma[,1]),]

subtype.s41<-c(as.character(check.match[which(check.match[,9]=="st4"),9]),as.character(check.match[which(check.match[,9]=="st1"),9]))
subtype.s34<-c(as.character(check.match[which(check.match[,9]=="st3"),9]),as.character(check.match[which(check.match[,9]=="st4"),9]))

Draw_heatmap(ed.set123.frma.cancer.reorder.st41.MDE.100,"heatmap_allsample_frma_MDE_200_probes_41.pdf","/media/H_driver/2015/Sophia/Results/",subtype.s41)
Draw_heatmap(ed.set123.frma.cancer.reorder.st34.MDE.100,"heatmap_allsample_frma_MDE_200_probes_34.pdf","/media/H_driver/2015/Sophia/Results/",subtype.s34)


Draw_PCA(ed.set123.frma.cancer.reorder,"PCA_58_cancer_sample_frma_all_probes.pdf","/media/H_driver/2015/Sophia/Results/",check.match$subtype)

save(heatmap_wPCA,ed.set123.frma.cancer.reorder,check.match,file="Data_4_heatmap_PCA.RData")

#Use all data to perform GSEA analysis
GenerateFiles4GSEA<-function(ed.set123.frma.cancer,out_dir,out_file_name,subtypeA,subtypeB){

  index.st4<-which(colnames(ed.set123.frma.cancer) %in% cel.file.name.key.set.symbol.st134[which(cel.file.name.key.set.symbol.st134[,5]==subtypeA),2])
  index.st1<-which(colnames(ed.set123.frma.cancer) %in% cel.file.name.key.set.symbol.st134[which(cel.file.name.key.set.symbol.st134[,5]==subtypeB),2])

  ed.set123.frma.cancer.st41<-cbind(ed.set123.frma.cancer[,index.st4],ed.set123.frma.cancer[,index.st1])
  cat(subtypeA,"\t",length(index.st4),subtypeB,"\t",length(index.st1),"\t",dim(ed.set123.frma.cancer.st41),"\n")

  ed.set123.frma.cancer.st41.reformat<-cbind(rownames(ed.set123.frma.cancer.st41),rep("NA",dim(ed.set123.frma.cancer.st41)[1]),ed.set123.frma.cancer.st41)
  colnames(ed.set123.frma.cancer.st41.reformat)[c(1,2)]=c("NAME","Description")

  write.table(ed.set123.frma.cancer.st41.reformat,file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = TRUE)

}

GenerateFiles4GSEA(ed.set123.frma.cancer,"/media/H_driver/2015/Sophia/Results/","St41-4_GSEA.xls","st4","st1")
GenerateFiles4GSEA(ed.set123.frma.cancer,"/media/H_driver/2015/Sophia/Results/","St34-4_GSEA.xls","st3","st4")
GenerateFiles4GSEA(ed.set123.frma.cancer,"/media/H_driver/2015/Sophia/Results/","St13-4_GSEA.xls","st1","st3")

save.image(file="/media/H_driver/Aimin_project/Sophia/Data_set_3_16_2016_all_probes_based_New.RData")
savehistory(file="/media/H_driver/Aimin_project/Sophia/Data_set_3_16_2016_all_probes_based_New.Rhistory")
#load the affy library
.libPaths()
library(gdata)
library(affy)
library(affyio)
library(annotate)
library(plyr)

#biocLite("mogene20sttranscriptcluster.db")
#library(mogene20stprobeset.db)
#library(mogene20sttranscriptcluster.db)

#biocLite("hugene10sttranscriptcluster.db")
#library(hugene10sttranscriptcluster.db)

biocLite("pd.hg.u133.plus.2")
library(pd.hg.u133.plus.2)
biocLite("ggplot2")
library(ggplot2)
biocLite("reshape")
biocLite("plot3D")
biocLite("gplots")
biocLite("ggdendro")
biocLite("RColorBrewer")
biocLite("frma")

library("reshape")
library("plot3D")
library("gplots")
library(ggdendro)
library(RColorBrewer)
library(RColorBrewer)
library(frma)
biocLite("hgu133plus2frmavecs")
#R CMD INSTALL -l R/x86_64-pc-linux-gnu-library/3.2/ /media/H_driver/Aimin_project/hgu133plus2frmavecs_1.5.0.tar.gz
library(hgu133plus2frmavecs)
biocLite("Homo.sapiens")
library("Homo.sapiens")
library(hgu133plus2.db)

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

#data.ER.PR.sample.info.new<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx".sheet=2)

data.ER.PR.sample.info.sheet1<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx",sheet=1) # Label: Sheet2
data.ER.PR.sample.info.sheet2<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx",sheet=2) # Label: Sheet1

data.ER.PR.sample.info.sheet11<-read.table("/media/H_driver/2015/Sophia/Results/ER_PR_Subgroup_PT2_v2-2-sheet1.txt",header=T,sep="\t",as.is=T)
data.ER.PR.sample.info.sheet22<-read.table("/media/H_driver/2015/Sophia/Results/ER_PR_Subgroup_PT2_v2-2-sheet2.txt",header=T,sep="\t",as.is=T)

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

ReformatCelName<-function(data.set2.cancer.GSM.cel.mapping.3){
data.set2.cancer.GSM.cel.mapping.3.2<-gsub(".cel","",data.set2.cancer.GSM.cel.mapping.3)
data.set2.cancer.GSM.cel.mapping.3.2.2<-gsub(".CEL","",as.character(data.set2.cancer.GSM.cel.mapping.3.2))
data.set2.cancer.GSM.cel.mapping.3.2.2.2<-gsub("-","_",as.character(data.set2.cancer.GSM.cel.mapping.3.2.2))
data.set2.cancer.GSM.cel.mapping.4<-cbind(data.set2.cancer.GSM.cel.mapping.3,data.set2.cancer.GSM.cel.mapping.3.2.2.2,sapply(strsplit(data.set2.cancer.GSM.cel.mapping.3.2.2.2,"_"),"[[",1))
return(data.set2.cancer.GSM.cel.mapping.4)
}

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

cel.file.name.key.set.symbol.cutoff.50.old<-ClassifyCancerSamplesIntoSubType(cel.file.cancer.set1.set2.set3,data.ER.PR.sample.info,4)

#data.ER.PR.sample.info.sheet11[which(data.ER.PR.sample.info.sheet11[,4]==1),c(1,4)]

cel.file.name.key.set.symbol.cutoff.5<-ClassifyCancerSamplesIntoSubType(cel.file.cancer.set1.set2.set3,data.ER.PR.sample.info.sheet11,4)
cel.file.name.key.set.symbol.cutoff.10<-ClassifyCancerSamplesIntoSubType(cel.file.cancer.set1.set2.set3,data.ER.PR.sample.info.sheet11,7)
cel.file.name.key.set.symbol.cutoff.50<-ClassifyCancerSamplesIntoSubType(cel.file.cancer.set1.set2.set3,data.ER.PR.sample.info.sheet11,10)

#Use subtype 1,3,4
#cel.file.name.key.set.symbol.st134<-cel.file.name.key.set.symbol.subtype[which(cel.file.name.key.set.symbol.subtype[,5]!="st2"),]

cel.file.name.key.set.symbol.st134.2<-cel.file.name.key.set.symbol.cutoff.50.old[which(cel.file.name.key.set.symbol.cutoff.50.old[,5]!="st2"),]

#Use frma method to do normalization
setwd("/media/H_driver/2015/Sophia/Cel_file_frma/")
data.set123.raw <- ReadAffy()
data.set123.frma<-frma(data.set123.raw)
ed.set123.frma <- exprs(data.set123.frma)
samp.frma <- sampleNames(data.set123.frma)
probes.frma <- featureNames(data.set123.frma)

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

Re.cutoff.50.old<-MapingCelFile4ReadFram(samp.frma,data.set123,ed.set123.frma,cel.file.name.key.set.symbol.st134)

Re.cutoff.5<-MapingCelFile4ReadFram(samp.frma,data.set123,ed.set123.frma,cel.file.name.key.set.symbol.cutoff.5)
Re.cutoff.10<-MapingCelFile4ReadFram(samp.frma,data.set123,ed.set123.frma,cel.file.name.key.set.symbol.cutoff.10)
Re.cutoff.50.new<-MapingCelFile4ReadFram(samp.frma,data.set123,ed.set123.frma,cel.file.name.key.set.symbol.cutoff.50)

names(Re.cutoff.50.old)

GeneSym.all <- getSYMBOL(rownames(ed.set123.frma.cancer), "hgu133plus2.db")


#<-cbind(names(GeneSym.all),GeneSym.all)

# #Collapse probes into genes
ndata<-ed.set123.frma.cancer

geneSymbol=GeneSym.all
tempdata.byGSym =data.frame(ndata, Symbol = geneSymbol)

colnames(tempdata.byGSym)

tempdata.byGSym.2<-tempdata.byGSym[-which(is.na(tempdata.byGSym$Symbol)),]

# test.sample<-ed.set123.frma.cancer[which(GeneSym.all=="ERBB4"),which(colnames(ed.set123.frma.cancer) %in% cel.file.sample.infor.no.8[which(cel.file.sample.infor.no.8[,5]=="st4"),6])]
#
# test.sample.2<-data.frame(test.sample,Symbol="ERBB4")
# test.sample.3 = test.sample.2
# test.sample.3[,1:5] = apply(test.sample.3[,1:5],2,as.numeric)

data.byGSym = ddply(tempdata.byGSym.2, c("Symbol"),function(h)
  summary = apply(h[,1:(dim(tempdata.byGSym.2)[2]-1)],2,max)
)


rownames(data.byGSym)=data.byGSym$Symbol
colnames(data.byGSym)

data.byGSym.2<-data.byGSym[,-1]
temp = apply(data.byGSym.2,2,as.numeric)
rownames(temp) = rownames(data.byGSym.2)
colnames(temp) = colnames(ndata)

data.byGSym.2 = temp
colnames(data.byGSym.2)
#data.byGSym.2 = ed.set123.frma.cancer

SampleType<-cel.file.sample.infor.no.8$subtype
heatmap_wPCA(data.byGSym.2, "heatmap_allsample_frma_all_probes.pdf","PCA_allsample_frma_all_probes.pdf","/media/H_driver/2015/Sophia/Results/", SampleType)

f.st134.frma <- factor(cel.file.sample.infor.no.8$subtype)
design.st134.frma <- model.matrix(~0+f.st134.frma)
colnames(design.st134.frma) <- levels(f.st134.frma)

reorder.index<-match(cel.file.sample.infor.no.8[,6],colnames(data.byGSym.2))

data.byGSym.3<-data.byGSym.2[,reorder.index]

#colnames(data.byGSym.3)

fit.st134.gene.frma <- lmFit(data.byGSym.3, design.st134.frma)
cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134.frma)
fit2.st134.frma  <- contrasts.fit(fit.st134.gene.frma, cont.matrix.st134.frma)
fit2.st134.frma  <- eBayes(fit2.st134.frma)

TopTableSt34.gene.frma<-topTable(fit2.st134.frma,coef=1,n=dim(data.byGSym.2)[1])
TopTableSt41.gene.frma<-topTable(fit2.st134.frma,coef=2,n=dim(data.byGSym.2)[1])

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

OutPut2HtmlTable(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/Results/","St41-8-frma.html","st4-vs-st1_HTML_report")
OutPut2HtmlTable(TopTableSt34.gene.frma,"/media/H_driver/2015/Sophia/Results/","St34-8-frma.html","st3-vs-st4_HTML_report")

#Generate the rank file for GSEA
GenerateRankFile4GSEA<-function(TopTableSt41.gene.frma,out_dir,out_file_name,up_down_top){

  FC.sign=sign(TopTableSt41.gene.frma[,1])
  p.inve=1/TopTableSt41.gene.frma[,4]
  Gene.Rank.Score=FC.sign*p.inve

  TopTableSt41.gene.frma.GSEA<-cbind(TopTableSt41.gene.frma,Gene.Rank.Score)

  TopTableSt41.gene.frma.GSEA.sorted.by.score<-TopTableSt41.gene.frma.GSEA[order(TopTableSt41.gene.frma.GSEA[,7],decreasing = TRUE),]

  Re.GSEA<-cbind(rownames(TopTableSt41.gene.frma.GSEA.sorted.by.score),TopTableSt41.gene.frma.GSEA.sorted.by.score[,7])
  colnames(Re.GSEA)=c("GeneName","Score")

  Re.GSEA.2<-rbind(head(Re.GSEA,up_down_top),tail(Re.GSEA,up_down_top))

  write.table(Re.GSEA.2, file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = FALSE)

  return(TopTableSt41.gene.frma.GSEA.sorted.by.score)
}

GenerateRankFile4GSEA(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/Results/","St41-4_GSEA_up_down_300.rnk",300)

#write.csv(TopTableSt41.gene.frma, file="/media/H_driver/2015/Sophia/Results/St41.csv",quote =FALSE)


##Use all probes for DE
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

DE4cut_off_50<-DeAnalysis(Re.cutoff.50.old)
DE4cut_off_5<-DeAnalysis(Re.cutoff.5)
DE4cut_off_10<-DeAnalysis(Re.cutoff.10)
DE4cut_off_50_new<-DeAnalysis(Re.cutoff.50.new)

factor(Re.cutoff.5[[1]]$subtype)
model.matrix(~0+factor(Re.cutoff.5[[1]]$subtype))
colnames() <- levels(f.st134.frma)

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

#5% DE
RE34.cut5<-OutPut2HtmlTableAllProbes(DE4cut_off_5[[1]],"/media/H_driver/2015/Sophia/Results/","cutoff5-St34-all-probes-frma-3.html","cutoff5_st3-vs-st4_HTML_report")
RE41.cut5<-OutPut2HtmlTableAllProbes(DE4cut_off_5[[2]],"/media/H_driver/2015/Sophia/Results/","cutoff5-St41-all-probes-frma-3.html","cutoff5_st4-vs-st1_HTML_report")
RE13.cut5<-OutPut2HtmlTableAllProbes(DE4cut_off_5[[3]],"/media/H_driver/2015/Sophia/Results/","cutoff5-St13-all-probes-frma.html","cutoff5_st1-vs-st3_HTML_report")
RE12.cut5<-OutPut2HtmlTableAllProbes(DE4cut_off_5[[4]],"/media/H_driver/2015/Sophia/Results/","cutoff5-St12-all-probes-frma-3.html","cutoff5_st1-vs-st2_HTML_report")
RE23.cut5<-OutPut2HtmlTableAllProbes(DE4cut_off_5[[5]],"/media/H_driver/2015/Sophia/Results/","cutoff5-St23-all-probes-frma-3.html","cutoff5_st2-vs-st3_HTML_report")
RE24.cut5<-OutPut2HtmlTableAllProbes(DE4cut_off_5[[6]],"/media/H_driver/2015/Sophia/Results/","cutoff5-St24-all-probes-frma.html","cutoff5_st2-vs-st4_HTML_report")

#10% DE
RE34.cut10<-OutPut2HtmlTableAllProbes(DE4cut_off_10[[1]],"/media/H_driver/2015/Sophia/Results/","cutoff10-St34-all-probes-frma-3.html","cutoff10_st3-vs-st4_HTML_report")
RE41.cut10<-OutPut2HtmlTableAllProbes(DE4cut_off_10[[2]],"/media/H_driver/2015/Sophia/Results/","cutoff10-St41-all-probes-frma-3.html","cutoff10_st4-vs-st1_HTML_report")
RE13.cut10<-OutPut2HtmlTableAllProbes(DE4cut_off_10[[3]],"/media/H_driver/2015/Sophia/Results/","cutoff10-St13-all-probes-frma.html","cutoff10_st1-vs-st3_HTML_report")
RE12.cut10<-OutPut2HtmlTableAllProbes(DE4cut_off_10[[4]],"/media/H_driver/2015/Sophia/Results/","cutoff10-St12-all-probes-frma-3.html","cutoff10_st1-vs-st2_HTML_report")
RE23.cut10<-OutPut2HtmlTableAllProbes(DE4cut_off_10[[5]],"/media/H_driver/2015/Sophia/Results/","cutoff10-St23-all-probes-frma-3.html","cutoff10_st2-vs-st3_HTML_report")
RE24.cut10<-OutPut2HtmlTableAllProbes(DE4cut_off_10[[6]],"/media/H_driver/2015/Sophia/Results/","cutoff-St24-all-probes-frma.html","cutoff10_st2-vs-st4_HTML_report")

#50% DE
RE41<-OutPut2HtmlTableAllProbes(TopTableSt41.all.probes.frma,"/media/H_driver/2015/Sophia/Results/","St41-all-probes-frma-3.html","st4-vs-st1_HTML_report")
RE34<-OutPut2HtmlTableAllProbes(TopTableSt34.all.probes.frma,"/media/H_driver/2015/Sophia/Results/","St34-all-probes-frma-3.html","st3-vs-st4_HTML_report")
RE13<-OutPut2HtmlTableAllProbes(TopTableSt13.all.probes.frma,"/media/H_driver/2015/Sophia/Results/","St13-all-probes-frma.html","st1-vs-st3_HTML_report")
RE41<-OutPut2HtmlTableAllProbes(TopTableSt41.all.probes.frma,"/media/H_driver/2015/Sophia/Results/","St41-all-probes-frma-3.html","st4-vs-st1_HTML_report")
RE34<-OutPut2HtmlTableAllProbes(TopTableSt34.all.probes.frma,"/media/H_driver/2015/Sophia/Results/","St34-all-probes-frma-3.html","st3-vs-st4_HTML_report")
RE13<-OutPut2HtmlTableAllProbes(TopTableSt13.all.probes.frma,"/media/H_driver/2015/Sophia/Results/","St13-all-probes-frma.html","st1-vs-st3_HTML_report")

#head(RE41)

GenerateRank4AllProbes<-function(TopTableSt41.all.probes.frma,up_down_top){

  FC.sign=sign(TopTableSt41.all.probes.frma[,1])
  p.inve=1/TopTableSt41.all.probes.frma[,4]
  Gene.Rank.Score=FC.sign*p.inve

  TopTableSt41.gene.frma.GSEA<-cbind(TopTableSt41.all.probes.frma,Gene.Rank.Score)

  TopTableSt41.gene.frma.GSEA.sorted.by.score<-TopTableSt41.gene.frma.GSEA[order(TopTableSt41.gene.frma.GSEA[,7],decreasing = TRUE),]

  Re.GSEA<-cbind(rownames(TopTableSt41.gene.frma.GSEA.sorted.by.score),TopTableSt41.gene.frma.GSEA.sorted.by.score[,7])
  colnames(Re.GSEA)=c("Probe","Score")

  Re.GSEA.2<-rbind(head(Re.GSEA,up_down_top),tail(Re.GSEA,up_down_top))

  #write.table(Re.GSEA.2, file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = FALSE)

  return(Re.GSEA.2)
}


TopTableSt41.most.DE.probes.frma<-GenerateRank4AllProbes(TopTableSt41.all.probes.frma,200)
TopTableSt34.most.DE.probes.frma<-GenerateRank4AllProbes(TopTableSt34.all.probes.frma,200)

ed.set123.frma.cancer.reorder.st1<-ed.set123.frma.cancer.reorder[,which(colnames(ed.set123.frma.cancer.reorder) %in% check.match[which(check.match$st1==1),1])]
ed.set123.frma.cancer.reorder.st2<-ed.set123.frma.cancer.reorder[,which(colnames(ed.set123.frma.cancer.reorder) %in% check.match[which(check.match$st2==1),1])]
ed.set123.frma.cancer.reorder.st3<-ed.set123.frma.cancer.reorder[,which(colnames(ed.set123.frma.cancer.reorder) %in% check.match[which(check.match$st3==1),1])]
ed.set123.frma.cancer.reorder.st4<-ed.set123.frma.cancer.reorder[,which(colnames(ed.set123.frma.cancer.reorder) %in% check.match[which(check.match$st4==1),1])]

dim(ed.set123.frma.cancer.reorder.st1)
dim(ed.set123.frma.cancer.reorder.st2)
dim(ed.set123.frma.cancer.reorder.st3)
dim(ed.set123.frma.cancer.reorder.st4)

ed.set123.frma.cancer.reorder.st41<-cbind(ed.set123.frma.cancer.reorder.st4,ed.set123.frma.cancer.reorder.st1)
ed.set123.frma.cancer.reorder.st34<-cbind(ed.set123.frma.cancer.reorder.st3,ed.set123.frma.cancer.reorder.st4)

dim(ed.set123.frma.cancer.reorder.st41)
dim(ed.set123.frma.cancer.reorder.st34)

ed.set123.frma.cancer.reorder.st41.MDE.100<-ed.set123.frma.cancer.reorder.st41[which(rownames(ed.set123.frma.cancer.reorder.st41) %in%  TopTableSt41.most.DE.probes.frma[,1]),]
ed.set123.frma.cancer.reorder.st34.MDE.100<-ed.set123.frma.cancer.reorder.st34[which(rownames(ed.set123.frma.cancer.reorder.st34) %in%  TopTableSt34.most.DE.probes.frma[,1]),]

subtype.s41<-c(as.character(check.match[which(check.match[,9]=="st4"),9]),as.character(check.match[which(check.match[,9]=="st1"),9]))
subtype.s34<-c(as.character(check.match[which(check.match[,9]=="st3"),9]),as.character(check.match[which(check.match[,9]=="st4"),9]))

Draw_heatmap(ed.set123.frma.cancer.reorder.st41.MDE.100,"heatmap_allsample_frma_MDE_200_probes_41.pdf","/media/H_driver/2015/Sophia/Results/",subtype.s41)
Draw_heatmap(ed.set123.frma.cancer.reorder.st34.MDE.100,"heatmap_allsample_frma_MDE_200_probes_34.pdf","/media/H_driver/2015/Sophia/Results/",subtype.s34)

Draw_PCA(ed.set123.frma.cancer.reorder,"PCA_58_cancer_sample_frma_all_probes.pdf","/media/H_driver/2015/Sophia/Results/",check.match$subtype)

save(heatmap_wPCA,ed.set123.frma.cancer.reorder,check.match,file="Data_4_heatmap_PCA.RData")

#Use all data to perform GSEA analysis
GenerateFiles4GSEA<-function(ed.set123.frma.cancer,cel.file.name.key.set.symbol.st134,cutoff,out_dir,out_file_name,subtypeA,subtypeB){

  index.st4<-which(colnames(ed.set123.frma.cancer) %in% cel.file.name.key.set.symbol.st134[which(cel.file.name.key.set.symbol.st134[,5]==subtypeA),2])
  index.st1<-which(colnames(ed.set123.frma.cancer) %in% cel.file.name.key.set.symbol.st134[which(cel.file.name.key.set.symbol.st134[,5]==subtypeB),2])

  ed.set123.frma.cancer.st41<-cbind(ed.set123.frma.cancer[,index.st4],ed.set123.frma.cancer[,index.st1])
  cat(subtypeA,"\t",length(index.st4),subtypeB,"\t",length(index.st1),"\t",dim(ed.set123.frma.cancer.st41),"\n")

  sink(paste0(out_dir,cutoff,"_",subtypeA,"vs",subtypeB,"_phenotype.cls"))
  cat(dim(ed.set123.frma.cancer.st41)[2],2,1,"\n")
  cat("#",subtypeA,subtypeB,"\n")
  cat(c(rep(0,length(index.st4)),rep(1,length(index.st1))),"\n")
  sink()

  ed.set123.frma.cancer.st41.reformat<-cbind(rownames(ed.set123.frma.cancer.st41),rep("NA",dim(ed.set123.frma.cancer.st41)[1]),ed.set123.frma.cancer.st41)
  colnames(ed.set123.frma.cancer.st41.reformat)[c(1,2)]=c("NAME","Description")

  cat("#1.2","\n",file=paste0(out_dir,cutoff,"_",out_file_name), sep="",append=TRUE)
  cat(dim(ed.set123.frma.cancer.st41)[1]," ",dim(ed.set123.frma.cancer.st41)[2], "\n",file=paste0(out_dir,cutoff,"_",out_file_name), sep="",append=TRUE)
  write.table(ed.set123.frma.cancer.st41.reformat,file = paste0(out_dir,cutoff,"_",out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = TRUE,append=TRUE)

}

GenerateFiles4GSEA(ed.set123.frma.cancer,"/media/H_driver/2015/Sophia/Results/","St41-4_GSEA.xls","st4","st1")
GenerateFiles4GSEA(ed.set123.frma.cancer,"/media/H_driver/2015/Sophia/Results/","St34-4_GSEA.xls","st3","st4")
GenerateFiles4GSEA(ed.set123.frma.cancer,"/media/H_driver/2015/Sophia/Results/","St13-4_GSEA.xls","st1","st3")

#5%
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.5,5,"/media/H_driver/2015/Sophia/Results/","St34.gct","st3","st4")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.5,5,"/media/H_driver/2015/Sophia/Results/","St41.gct","st4","st1")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.5,5,"/media/H_driver/2015/Sophia/Results/","St13.gct","st1","st3")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.5,5,"/media/H_driver/2015/Sophia/Results/","St12.gct","st1","st2")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.5,5,"/media/H_driver/2015/Sophia/Results/","St23.gct","st2","st3")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.5,5,"/media/H_driver/2015/Sophia/Results/","St24.gct","st2","st4")

#10%
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10,10,"/media/H_driver/2015/Sophia/Results/","St34.gct","st3","st4")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10,10,"/media/H_driver/2015/Sophia/Results/","St41.gct","st4","st1")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10,10,"/media/H_driver/2015/Sophia/Results/","St13.gct","st1","st3")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10,10,"/media/H_driver/2015/Sophia/Results/","St12.gct","st1","st2")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10,10,"/media/H_driver/2015/Sophia/Results/","St23.gct","st2","st3")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10,10,"/media/H_driver/2015/Sophia/Results/","St24.gct","st2","st4")

savehistory(file="/media/H_driver/Aimin_project/Sophia/Data_set_5_6_2016_all_probes_based_New.Rhistory")
save.image(file="/media/H_driver/Aimin_project/Sophia/Data_set_5_6_2016_all_probes_based_New.RData")

#savehistory(file="/media/H_driver/Aimin_project/Sophia/Data_set_3_30_2016_all_probes_based_New.Rhistory")
#function.all<-ll(class="function")
#save(function.all,file="Sophia.R")


cel.file.cancer.set3.only<-cel.file.cancer.set1.set2.set3[which(cel.file.cancer.set1.set2.set3[,3]==3),]

cel.file.name.key.set.symbol.cutoff.5.set3<-ClassifyCancerSamplesIntoSubType(cel.file.cancer.set3.only,data.ER.PR.sample.info.sheet11,4)
cel.file.name.key.set.symbol.cutoff.10.set3<-ClassifyCancerSamplesIntoSubType(cel.file.cancer.set3.only,data.ER.PR.sample.info.sheet11,7)
cel.file.name.key.set.symbol.cutoff.50.set3<-ClassifyCancerSamplesIntoSubType(cel.file.cancer.set3.only,data.ER.PR.sample.info.sheet11,10)

cel.file.name.key.set.symbol.cutoff.5.set3[,5]
cel.file.name.key.set.symbol.cutoff.10.set3[,5]
cel.file.name.key.set.symbol.cutoff.50.set3[,5]

Re.cutoff.5.set3<-MapingCelFile4ReadFram(samp.frma,data.set123,ed.set123.frma,cel.file.name.key.set.symbol.cutoff.5.set3)
Re.cutoff.10.set3<-MapingCelFile4ReadFram(samp.frma,data.set123,ed.set123.frma,cel.file.name.key.set.symbol.cutoff.10.set3)
Re.cutoff.50.new.set3<-MapingCelFile4ReadFram(samp.frma,data.set123,ed.set123.frma,cel.file.name.key.set.symbol.cutoff.50.set3)

length(Re.cutoff.5.set3)
head(Re.cutoff.5.set3[[1]])

DE4cut_off_5_set3<-DeAnalysis2(Re.cutoff.5.set3)
DE4cut_off_10_set3<-DeAnalysis2(Re.cutoff.10.set3)
DE4cut_off_50_new_set3<-DeAnalysis2(Re.cutoff.50.new.set3)

names(DE4cut_off_10_set3)
names(DE4cut_off_50_new_set3)

# Output DE analysis results
OuptuPutDE<-function(DE4cut_off_5){

name.object=deparse(substitute(DE4cut_off_5))

value.cutoff=unlist(strsplit(name.object, "_"))[3]

num.DE<-length(DE4cut_off_5)-1

for(i in 1:num.DE){

  cat(names(DE4cut_off_5)[i],"\n")

  name.DE.pairs=names(DE4cut_off_5)[i]

  OutPut2HtmlTableAllProbes(DE4cut_off_5[[i]],"/media/H_driver/2015/Sophia/Results/Results_Data_Set3/",paste0("cutoff-",value.cutoff,"-",name.DE.pairs,"-all-probes-frma-5.html"),
                            paste0("cutoff-",value.cutoff,"-",name.DE.pairs,"_HTML_report"))
}

}

OuptuPutDE(DE4cut_off_5_set3)
OuptuPutDE(DE4cut_off_10_set3)
OuptuPutDE(DE4cut_off_50_new_set3)

#5%
cel.file.name.key.set.symbol.cutoff.5.data.set3.only<-cel.file.name.key.set.symbol.cutoff.5[which(cel.file.name.key.set.symbol.cutoff.5[,3]==3),]

GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.5.data.set3.only,5,"/media/H_driver/2015/Sophia/Results/Results_Data_Set3/","St34.gct","st3","st4")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.5.data.set3.only,5,"/media/H_driver/2015/Sophia/Results/Results_Data_Set3/","St41.gct","st4","st1")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.5.data.set3.only,5,"/media/H_driver/2015/Sophia/Results/Results_Data_Set3/","St13.gct","st1","st3")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.5.data.set3.only,5,"/media/H_driver/2015/Sophia/Results/Results_Data_Set3/","St12.gct","st1","st2")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.5.data.set3.only,5,"/media/H_driver/2015/Sophia/Results/Results_Data_Set3/","St23.gct","st2","st3")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.5.data.set3.only,5,"/media/H_driver/2015/Sophia/Results/Results_Data_Set3/","St24.gct","st2","st4")

#10%
cel.file.name.key.set.symbol.cutoff.10.data.set3.only<-cel.file.name.key.set.symbol.cutoff.10[which(cel.file.name.key.set.symbol.cutoff.10[,3]==3),]

GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10.data.set3.only,10,"/media/H_driver/2015/Sophia/Results/Results_Data_Set3/","St34.gct","st3","st4")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10.data.set3.only,10,"/media/H_driver/2015/Sophia/Results/Results_Data_Set3/","St41.gct","st4","st1")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10.data.set3.only,10,"/media/H_driver/2015/Sophia/Results/Results_Data_Set3/","St13.gct","st1","st3")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10.data.set3.only,10,"/media/H_driver/2015/Sophia/Results/Results_Data_Set3/","St12.gct","st1","st2")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10.data.set3.only,10,"/media/H_driver/2015/Sophia/Results/Results_Data_Set3/","St23.gct","st2","st3")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10.data.set3.only,10,"/media/H_driver/2015/Sophia/Results/Results_Data_Set3/","St24.gct","st2","st4")

#50%
cel.file.name.key.set.symbol.cutoff.50.data.set3.only<-cel.file.name.key.set.symbol.cutoff.50[which(cel.file.name.key.set.symbol.cutoff.50[,3]==3),]

#GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10.data.set3.only,10,"/media/H_driver/2015/Sophia/Results/Results_Data_Set3/","St34.gct","st3","st4")
#GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10.data.set3.only,10,"/media/H_driver/2015/Sophia/Results/Results_Data_Set3/","St41.gct","st4","st1")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.50.data.set3.only,50,"/media/H_driver/2015/Sophia/Results/Results_Data_Set3/","St13.gct","st1","st3")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.50.data.set3.only,50,"/media/H_driver/2015/Sophia/Results/Results_Data_Set3/","St12.gct","st1","st2")
GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.50.data.set3.only,50,"/media/H_driver/2015/Sophia/Results/Results_Data_Set3/","St23.gct","st2","st3")
#GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10.data.set3.only,10,"/media/H_driver/2015/Sophia/Results/Results_Data_Set3/","St24.gct","st2","st4")

# cel.file.name.key.set.symbol.cutoff.50
# GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10,10,"/media/H_driver/2015/Sophia/Results/","St34.gct","st3","st4")
# GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10,10,"/media/H_driver/2015/Sophia/Results/","St41.gct","st4","st1")
# GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10,10,"/media/H_driver/2015/Sophia/Results/","St13.gct","st1","st3")
# GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10,10,"/media/H_driver/2015/Sophia/Results/","St12.gct","st1","st2")
# GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10,10,"/media/H_driver/2015/Sophia/Results/","St23.gct","st2","st3")
# GenerateFiles4GSEA(ed.set123.frma.cancer,cel.file.name.key.set.symbol.cutoff.10,10,"/media/H_driver/2015/Sophia/Results/","St24.gct","st2","st4")
#load the affy library
.libPaths()
library(gdata)
library(affy)
library(affyio)
library(annotate)
library(plyr)
biocLite("mogene20sttranscriptcluster.db")
library(mogene20stprobeset.db)
library(mogene20sttranscriptcluster.db)

biocLite("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)

biocLite("pd.hg.u133.plus.2")
library(pd.hg.u133.plus.2)
biocLite("ggplot2")
library(ggplot2)
biocLite("reshape")
biocLite("plot3D")
biocLite("gplots")
biocLite("ggdendro")
biocLite("RColorBrewer")
biocLite("frma")

library("reshape")
library("plot3D")
library("gplots")
library(ggdendro)
library(RColorBrewer)
library(RColorBrewer)
library(frma)
biocLite("hgu133plus2frmavecs")
#R CMD INSTALL -l R/x86_64-pc-linux-gnu-library/3.2/ /media/H_driver/Aimin_project/hgu133plus2frmavecs_1.5.0.tar.gz
library(hgu133plus2frmavecs)
biocLite("Homo.sapiens")
library("Homo.sapiens")

# Read in the CEL files in the directory, then normalize the data
setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM")
data.set1 <- ReadAffy()
ed.raw.set1 <- exprs(data.set1)
samp.set1 <- sampleNames(data.set1)
probes.set1 <- featureNames(data.set1)
ed.normalized.set1<- rma(data.set1)

setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set2_corrected")
data.set2 <- ReadAffy()
ed.raw.set2 <- exprs(data.set2)
samp.set2 <- sampleNames(data.set2)
probes.set2 <- featureNames(data.set2)
ed.normalized.set2<- rma(data.set2)

setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set3_HGSCancer_only")
data.set3 <- ReadAffy()
ed.raw.set3 <- exprs(data.set3)
samp.set3 <- sampleNames(data.set3)
probes.set3 <- featureNames(data.set3)
ed.normalized.set3<- rma(data.set3)

data.set1.normalized<-ed.normalized.set1
data.set2.normalized<-ed.normalized.set2
data.set3.normalized<-ed.normalized.set3

ll(dim=T)

colnames(data.set1.normalized)
colnames(data.set2.normalized)
colnames(data.set3.normalized)

length(colnames(data.set1.normalized))
length(colnames(data.set2.normalized))
length(colnames(data.set3.normalized))

unique(colnames(data.set1.normalized))
length(unique(colnames(data.set1.normalized)))
length(unique(colnames(data.set2.normalized)))
length(unique(colnames(data.set3.normalized)))

data.set1.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/Microarray_Set1_FTE_HGSC_LCM/CASE_DESCRIPTION_SET1.xlsx")
data.set2.sample.info<-read.xls("/home/aiminyan/Downloads/GSE28044_George.xlsx")
data.ER.PR.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx")
head(data.ER.PR.sample.info)

head(data.set1.sample.info)
data.set1.sample.info[,2]
data.set1.normalize
data.set2.normalized

grep(2367,unique(data.set1.sample.info[,2]))

unique(data.set2.sample.info[,2])
head(data.ER.PR.sample.info)

dim(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==1),])
dim(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==2),])
dim(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==3),])
dim(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==4),])

MapSample2CELfile<-function(data.ER.PR.sample.info,data.set1.normalized,sample_type){

cat(dim(data.ER.PR.sample.info),dim(data.set1.normalized),sample_type,"\n")

num.sample1=length(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==1),][,1])
cel_file_index_4_sample<-array()

for(i in 1:num.sample1) {
tmp<-grep(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==sample_type),][i,1],colnames(data.set1.normalized))
cel_file_index_4_sample<-c(cel_file_index_4_sample,tmp)
}
cel_file_index_4_sample<-cel_file_index_4_sample[-1]

colnames(data.set1.normalized)[cel_file_index_4_sample]

re<-data.set1.normalized[,cel_file_index_4_sample]

return(re)

}

#ER-PR-
ER.PR.sample1.from.set1<-MapSample2CELfile(data.ER.PR.sample.info,data.set1.normalized,1)
ER.PR.sample1.from.set2<-MapSample2CELfile(data.ER.PR.sample.info,data.set2.normalized,1)
ER.PR.sample1.from.set3<-MapSample2CELfile(data.ER.PR.sample.info,data.set3.normalized,1)

#ER-PR+
ER.PR.sample2.from.set1<-MapSample2CELfile(data.ER.PR.sample.info,data.set1.normalized,2)
ER.PR.sample2.from.set2<-MapSample2CELfile(data.ER.PR.sample.info,data.set2.normalized,2)
ER.PR.sample2.from.set3<-MapSample2CELfile(data.ER.PR.sample.info,data.set3.normalized,2)

#ER+PR-
ER.PR.sample3.from.set1<-MapSample2CELfile(data.ER.PR.sample.info,data.set1.normalized,3)
ER.PR.sample3.from.set2<-MapSample2CELfile(data.ER.PR.sample.info,data.set2.normalized,3)
ER.PR.sample3.from.set3<-MapSample2CELfile(data.ER.PR.sample.info,data.set3.normalized,3)

#ER+PR+
ER.PR.sample4.from.set1<-MapSample2CELfile(data.ER.PR.sample.info,data.set1.normalized,4)
ER.PR.sample4.from.set2<-MapSample2CELfile(data.ER.PR.sample.info,data.set2.normalized,4)
ER.PR.sample4.from.set3<-MapSample2CELfile(data.ER.PR.sample.info,data.set3.normalized,4)

head(ER.PR.sample1.from.set1)
head(ER.PR.sample1.from.set2)
head(ER.PR.sample1.from.set3)

head(ER.PR.sample3.from.set1)
head(ER.PR.sample3.from.set2)
head(ER.PR.sample3.from.set3)

head(data.set1.sample.info)
head(data.set2.sample.info)
head(data.ER.PR.sample.info)

data.set1.sample.info[,2]

data.set2.sample.info[,8:9:10]


cel.file.all<-c(unique(colnames(data.set1.normalized)),unique(colnames(data.set2.normalized)),unique(colnames(data.set3.normalized)))

length(cel.file.all)
cel.file.all

cel.file.all.2<-rbind(cbind(unique(colnames(data.set1.normalized)),rep(1,length(unique(colnames(data.set1.normalized))))),
cbind(unique(colnames(data.set2.normalized)),rep(2,length(unique(colnames(data.set2.normalized))))),
cbind(unique(colnames(data.set3.normalized)),rep(3,length(unique(colnames(data.set3.normalized))))))

cel.file.cancer.data.set2<-as.character(data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)][,8])
cel.file.cancer.data.set2.2<-as.character(data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)][,1])
cel.file.cancer.data.set2.2.2<-cel.file.cancer.data.set2.2[-which(cel.file.cancer.data.set2.2=="")]

data.set2.cancer.GSM.cel.mapping<-data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)]
data.set2.cancer.GSM.cel.mapping.2<-data.set2.cancer.GSM.cel.mapping[,c(1,8)]
data.set2.cancer.GSM.cel.mapping.3<-data.set2.cancer.GSM.cel.mapping.2[which(as.character(data.set2.cancer.GSM.cel.mapping.2[,1])!=""),]

data.set2.cancer.GSM.cel.mapping.3.2<-gsub(".cel","",as.character(data.set2.cancer.GSM.cel.mapping.3[,2]))
data.set2.cancer.GSM.cel.mapping.3.2.2<-gsub(".CEL","",as.character(data.set2.cancer.GSM.cel.mapping.3.2))
data.set2.cancer.GSM.cel.mapping.3.2.2.2<-gsub("-","_",as.character(data.set2.cancer.GSM.cel.mapping.3.2.2))
data.set2.cancer.GSM.cel.mapping.4<-cbind(data.set2.cancer.GSM.cel.mapping.3,sapply(strsplit(data.set2.cancer.GSM.cel.mapping.3.2.2.2,"_"),"[[",1))

cel.file.cancer.data.set1<-as.character(data.set1.sample.info[which(data.set1.sample.info[,4]=="cancer"),2])

cel.file.cancer.set1.set2<-rbind(cbind(cel.file.cancer.data.set1,rep(1,length(cel.file.cancer.data.set1))),
cbind(cel.file.cancer.data.set2,rep(2,length(cel.file.cancer.data.set2))),
cbind(cel.file.cancer.data.set2.2.2,rep(2,length(cel.file.cancer.data.set2.2.2))))

cel.file.cancer.set1.set2.set3<-rbind(cel.file.cancer.set1.set2,cel.file.all.2[which(cel.file.all.2[,2]==3),])

cel.file.cancer.set1.set2.set3.name<-gsub(".CEL","",cel.file.cancer.set1.set2.set3[,1])
cel.file.cancer.set1.set2.set3.name.1<-gsub(".cel","",cel.file.cancer.set1.set2.set3.name)
cel.file.cancer.set1.set2.set3.name.2<-gsub("-","_",cel.file.cancer.set1.set2.set3.name.1)

data.set123.normalized<-cbind(data.set1.normalized,data.set2.normalized,data.set3.normalized)

dim(data.set123.normalized)
dim(data.set123.normalized[,which(colnames(data.set123.normalized) %in% cel.file.cancer.set1.set2.set3[,1])])
colnames(data.set123.normalized)[grep(54,colnames(data.set123.normalized))]

original.cel.file.names<-colnames(data.set123.normalized)
original.cel.file.names.1<-gsub("X","",original.cel.file.names)
original.cel.file.names.2<-gsub(".CEL","",original.cel.file.names.1)
original.cel.file.names.3<-gsub(".cel","",original.cel.file.names.2)
original.cel.file.names.4<-gsub("\\.","_",original.cel.file.names.3)

index.4.cancer.sample<-which(original.cel.file.names.4 %in% cel.file.cancer.set1.set2.set3.name.2)
original.cel.file.name.with.mapped.names<-cbind(original.cel.file.names[index.4.cancer.sample],original.cel.file.names.4[index.4.cancer.sample])
#index.4.no.cancer.sample.1<-cbind(original.cel.file.names[-index.4.cancer.sample],original.cel.file.names.4[-index.4.cancer.sample])
setdiff(cel.file.cancer.set1.set2.set3.name.2,original.cel.file.names.4[index.4.cancer.sample])
which(original.cel.file.names.4[-index.4.cancer.sample] %in% cel.file.cancer.set1.set2.set3.name.2)

cancer.data.set123.normalized<-data.set123.normalized[,index.4.cancer.sample]
dim(cancer.data.set123.normalized)

cancer.cel.file.name.72<-colnames(cancer.data.set123.normalized)
cancer.cel.file.name.72.1<-gsub("X","",cancer.cel.file.name.72)
cancer.cel.file.name.72.2<-gsub(".CEL","",cancer.cel.file.name.72.1)
cancer.cel.file.name.72.3<-gsub(".cel","",cancer.cel.file.name.72.2)
cancer.cel.file.name.72.4<-gsub("\\.","_",cancer.cel.file.name.72.3)
cancer.cel.file.name.72.5<-cancer.cel.file.name.72.4

data.set2.cancer.GSM.cel.mapping.44<-as.character(data.set2.cancer.GSM.cel.mapping.4[which(data.set2.cancer.GSM.cel.mapping.4[,1] %in% cancer.cel.file.name.72.4),3])
cancer.cel.file.name.72.5[which(cancer.cel.file.name.72.5 %in% data.set2.cancer.GSM.cel.mapping.4[,1])]<-data.set2.cancer.GSM.cel.mapping.44
cancer.cel.file.name.72.6<-cbind(cancer.cel.file.name.72,cancer.cel.file.name.72.5)
cancer.cel.file.name.72.7<-c(sapply(strsplit(cancer.cel.file.name.72.6[1:29,2],"_"),"[[",1),
sapply(strsplit(cancer.cel.file.name.72.6[30:72,2],"_"),"[[",4))
cancer.cel.file.name.72.8<-cbind(cancer.cel.file.name.72.6,cancer.cel.file.name.72.7)

#Classify cancer patients to 4 subtypes

subtype.1<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==1),1])
subtype.2<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==2),1])
subtype.3<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==3),1])
subtype.4<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==4),1])

cancer.data.set123.normalized.st1<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.1)]
cancer.data.set123.normalized.st2<-as.data.frame(cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.2)])
colnames(cancer.data.set123.normalized.st2)<-colnames(cancer.data.set123.normalized)[which(cancer.cel.file.name.72.8[,3] %in% subtype.2)]
cancer.data.set123.normalized.st3<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.3)]
cancer.data.set123.normalized.st4<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.4)]

dim(cancer.data.set123.normalized.st1)
dim(cancer.data.set123.normalized.st2)
dim(cancer.data.set123.normalized.st3)
dim(cancer.data.set123.normalized.st4)

head(cancer.data.set123.normalized.st1)
head(cancer.data.set123.normalized.st2)
head(cancer.data.set123.normalized.st3)
head(cancer.data.set123.normalized.st4)

colnames(cancer.data.set123.normalized.st1)
colnames(cancer.data.set123.normalized.st2)
colnames(cancer.data.set123.normalized.st3)
colnames(cancer.data.set123.normalized.st4)

#n.st<-length(colnames(cancer.data.set123.normalized.st1))

# #Use subtype1,2,3,4
# cel.file.sample.infor<-as.data.frame(rbind(
# cbind(colnames(cancer.data.set123.normalized.st1),rep("st1",length(colnames(cancer.data.set123.normalized.st1)))),
# cbind(colnames(cancer.data.set123.normalized.st2),rep("st2",length(colnames(cancer.data.set123.normalized.st2)))),
# cbind(colnames(cancer.data.set123.normalized.st3),rep("st3",length(colnames(cancer.data.set123.normalized.st3)))),
# cbind(colnames(cancer.data.set123.normalized.st4),rep("st4",length(colnames(cancer.data.set123.normalized.st4))))
# ))
# colnames(cel.file.sample.infor)=c("filename","subtype")

#Use subtype 1,3,4
cel.file.sample.infor.no.2<-as.data.frame(rbind(
  cbind(colnames(cancer.data.set123.normalized.st1),rep("st1",length(colnames(cancer.data.set123.normalized.st1)))),
  cbind(colnames(cancer.data.set123.normalized.st3),rep("st3",length(colnames(cancer.data.set123.normalized.st3)))),
  cbind(colnames(cancer.data.set123.normalized.st4),rep("st4",length(colnames(cancer.data.set123.normalized.st4))))
))

colnames(cel.file.sample.infor.no.2)=c("filename","subtype")
f.st134 <- factor(cel.file.sample.infor.no.2$subtype)
design.st134 <- model.matrix(~0+f.st134)
colnames(design.st134) <- levels(f.st134)

cancer.data.st134<-cbind(cancer.data.set123.normalized.st1,
cancer.data.set123.normalized.st3,
cancer.data.set123.normalized.st4)
head(cancer.data.st134)

GeneSym.all <- getSYMBOL(rownames(cancer.data.st134), "hgu133plus2.db")
ndata<-cancer.data.st134
geneSymbol=GeneSym.all
tempdata.byGSym = data.frame(ndata, Symbol = geneSymbol)

tempdata.byGSym.2<-tempdata.byGSym[-which(is.na(tempdata.byGSym$Symbol)),]

rownames(tempdata.byGSym.2) = NULL
data.byGSym = ddply(tempdata.byGSym.2, c("Symbol"),function(h)
  summary = apply(h,2,max)
)

#class(data.byGSym)

#which(is.na(data.byGSym$Symbol))
rownames(data.byGSym)=data.byGSym$Symbol

data.byGSym.2<-data.byGSym[,-61]

dim(data.byGSym.2)

#Use all probes

fit.st134 <- lmFit(cancer.data.st134, design.st134)

cont.matrix.st134 <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134)
fit2.st134  <- contrasts.fit(fit.st134, cont.matrix.st134)
fit2.st134  <- eBayes(fit2.st134)

str(fit2.st134)

TopTableSt34.all<-topTable(fit2.st134,coef=1,n=54674)

TopTableSt34.100<-topTable(fit2.st134,coef=1,adjust="fdr",n=100)

genenames <- as.character(rownames(TopTableSt34.54675))

length(genenames)

genenames

annotation(ed.normalized.set1)
annotation(ed.normalized.set2)
annotation(ed.normalized.set3)
library("hgu133plus2.db")

map <- getpAnnMa("ENTREZID", "hgu133plus2", load=TRUE, type=c("env", "db"))
class(map)
ll <- getEG(genenames,"hgu133plus2.db")
GeneSym <- getSYMBOL(genenames, "hgu133plus2.db")

Probe.gene.sym<-(cbind(genenames,GeneSym))

Probe.gene.sym

length(which(is.na(Probe.gene.sym[,2])))


tab <- data.frame(GeneSym, TopTableSt34.100)
tab <- data.frame(rownames(tab), tab)

colnames(tab)[1] <- c("Probe ID")
ll <- list(ll)
htmlpage(ll, filename="/media/H_driver/2015/Sophia/St34-4.html", title="HTML report",
         othernames=tab, table.head=c("Locus ID",colnames(tab)), table.center=TRUE, digits=6)

colnames(fit2.st134)
topTable(fit2.st134,coef=1)
topTable(fit2.st134,coef=1,adjust="fdr")

topTable(fit2.st134,coef=2)
topTable(fit2.st134,coef=2,adjust="fdr")

results.st134 <- decideTests(fit2.st134,p.value=0.99)

summary(results.st134)
vennDiagram(results.st134)

fit2.st134$genes

unigeneTopTableSt34 <- topTable(fit2.st134,coef=1,n=20,genelist=genelist)
unigeneTopTableSt41 <- topTable(fit2.st134,coef=2,n=20,genelist=genelist)

library(xtable)
xtableUnigeneSt34 <- xtable(unigeneTopTableSt34,display=c("s","s","s","s","g","g","g","e","e","g","g"))
xtableUnigeneSt41 <- xtable(unigeneTopTableSt41,display=c("s","s","s","s","g","g","g","e","e","g","g"))

cat(file="/media/H_driver/2015/Sophia/St34.html","<html>\n<body>")
print.xtable(xtableUnigeneSt34,type="html",file="/media/H_driver/2015/Sophia/St34.html",append=TRUE)
cat(file="/media/H_driver/2015/Sophia/St34.html","</body>\n</html>",append=TRUE)

cat(file="/media/H_driver/2015/Sophia/St41.html","<html>\n<body>")
print.xtable(xtableUnigeneSt41,type="html",file="/media/H_driver/2015/Sophia/St41.html",append=TRUE)
cat(file="/media/H_driver/2015/Sophia/St41.html","</body>\n</html>",append=TRUE)


#Use gene name
temp = apply(data.byGSym.2,2,as.numeric)
rownames(temp) = rownames(data.byGSym.2)
data.byGSym.2 = temp


fit.st134.gene <- lmFit(data.byGSym.2, design.st134)

cont.matrix.st134 <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134)
fit2.st134  <- contrasts.fit(fit.st134.gene, cont.matrix.st134)
fit2.st134  <- eBayes(fit2.st134)

str(fit2.st134)

TopTableSt34.gene<-topTable(fit2.st134,coef=1,n=54674)
dim(TopTableSt34.gene)
head(TopTableSt34.gene)
hist(TopTableSt34.gene[,4])

TopTableSt41.gene<-topTable(fit2.st134,coef=2,n=54674)
dim(TopTableSt41.gene)
head(TopTableSt41.gene)
hist(TopTableSt41.gene[,4])


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

SampleType = factor(gsub("(FTY).[1-3]","\\1",colnames(data.byGSym)))

SampleType<-cel.file.sample.infor.no.2$subtype

heatmap_wPCA(data.byGSym.2, "heatmap_allsample.pdf","PCA_allsample.pdf","/media/H_driver/2015/Sophia/", SampleType)


#Use frma method to do normalization
setwd("/media/H_driver/2015/Sophia/Cel_file_frma/")
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
heatmap_wPCA(data.byGSym.2, "heatmap_allsample_frma.pdf","PCA_allsample_frma.pdf","/media/H_driver/2015/Sophia/", SampleType)

f.st134.frma <- factor(cel.file.sample.infor.no.8$subtype)
design.st134.frma <- model.matrix(~0+f.st134.frma)
colnames(design.st134.frma) <- levels(f.st134.frma)

fit.st134.gene.frma <- lmFit(data.byGSym.2, design.st134.frma)
cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134.frma)
fit2.st134.frma  <- contrasts.fit(fit.st134.gene.frma, cont.matrix.st134.frma)
fit2.st134.frma  <- eBayes(fit2.st134.frma)

#str(fit2.st134)

TopTableSt34.gene.frma<-topTable(fit2.st134.frma,coef=1,n=60000)
dim(TopTableSt34.gene.frma)
head(TopTableSt34.gene.frma)
hist(TopTableSt34.gene.frma[,4])

TopTableSt41.gene.frma<-topTable(fit2.st134.frma,coef=2,n=60000)
dim(TopTableSt41.gene.frma)
head(TopTableSt41.gene.frma)
hist(TopTableSt41.gene.frma[,4])

#genenames <- as.character(rownames(TopTableSt34.54675))
#probe.id<-as.character(rownames(ed.set123.frma.cancer))

#length(genenames)

#genenames

#annotation(ed.set123.frma.cancer)
#annotation(ed.normalized.set2)
#annotation(ed.normalized.set3)
#library("hgu133plus2.db")
#map <- getpAnnMa("ENTREZID", "hgu133plus2", load=TRUE, type=c("env", "db"))
#map <- getAnnMa("ENTREZID", "hgu133plus2", load=TRUE, type=c("env", "db"))
#class(map)
#ll <- getEG(probe.id,"hgu133plus2.db")
#GeneSym <- getSYMBOL(probe.id,"hgu133plus2.db")

#getEG(GeneSym[1],"hgu133plus2.db")

#Probe.gene.sym<-(cbind(probe.id,GeneSym))

#Probe.gene.sym

#length(which(is.na(Probe.gene.sym[,2])))


#tab <- data.frame(GeneSym, TopTableSt34.100)
#tab <- data.frame(rownames(tab), tab)

#colnames(tab)[1] <- c("Probe ID")
#ll <- list(ll)


mart = useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org", dataset="mmusculus_gene_ensembl")

all.entrezgene <- unique( getBM(attributes = "entrezgene",values = "*", mart = ensembl) )
#G_list<-getBM(attributes = c("mgi_symbol", "entrezgene","ensembl_gene_id","chromosome_name","start_position","end_position"),
#              filter="ensembl_gene_id", values=genes, mart=mart)

library(biomaRt)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
entrezgene.gene.symbol<-getBM(attributes = c("entrezgene","hgnc_symbol"),filters="hgnc_symbol",values=rownames(TopTableSt41.gene.frma),mart=human)

hgnc.gene.symbol<-getBM(attributes = c("hgnc_id","hgnc_symbol"),filters="hgnc_symbol",values=rownames(TopTableSt41.gene.frma),mart=human)

ensembl.gene.symbol<-getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),filters="hgnc_symbol",values=rownames(TopTableSt41.gene.frma),mart=human)

hgnc.gene.symbol.ENTREZID<-select(Homo.sapiens, keys=rownames(TopTableSt41.gene.frma),columns=c("SYMBOL","ENTREZID"), keytype="SYMBOL")

head(hgnc.gene.symbol.ENTREZID)

TopTableSt41.gene.frma.2<-data.frame(rownames(TopTableSt41.gene.frma),TopTableSt41.gene.frma)

head(TopTableSt41.gene.frma.2)
colnames(TopTableSt41.gene.frma.2)[1]="SYMBOL"

TopTableSt41.gene.frma.3<-merge(TopTableSt41.gene.frma.2,hgnc.gene.symbol.ENTREZID,by="SYMBOL",sort=FALSE)
dim(TopTableSt41.gene.frma.3)
head(TopTableSt41.gene.frma.3)

htmlpage(list(TopTableSt41.gene.frma.3$ENTREZID),filename="/media/H_driver/2015/Sophia/St34-4-frma.html", title="HTML report",
         othernames=TopTableSt41.gene.frma.3, table.head=c("FG",colnames(TopTableSt41.gene.frma.3)), table.center=TRUE, digits=6)

TopTableSt41.gene.frma.4<-merge(TopTableSt41.gene.frma.2,hgnc.gene.symbol.ENTREZID,by="SYMBOL",sort=FALSE)
dim(TopTableSt41.gene.frma.4)
head(TopTableSt41.gene.frma.4)

htmlpage(list(TopTableSt41.gene.frma.4$ENTREZID),filename="/media/H_driver/2015/Sophia/St41-4-frma.html", title="st4-vs-st1_HTML_report",
         othernames=TopTableSt41.gene.frma.4[,-8], table.head=c("ENTREZID",colnames(TopTableSt41.gene.frma.3[,-8])), table.center=TRUE, digits=6)

htmlpage(filename="/media/H_driver/2015/Sophia/St34-4-frma.html", title="HTML report",
         othernames=TopTableSt34.gene.frma, table.head=c(colnames(TopTableSt34.gene.frma)), table.center=TRUE, digits=6)

library("Homo.sapiens")

# clustering plot
#fit#Differential gene expression analysis for st3 vs st4
#probe<-rownames(cancer.data.set123.normalized.st1)
#geneSymbol<-getSYMBOL(probe, "hugene10sttranscriptcluster.db")
#Differential gene expression analysis for st2 vs st1
#Differential gene expression analysis for st4 vs st1
save.image(file="/media/H_driver/Aimin_project/Sophia/Data_set_3.RData")#load the affy library
library(gdata)
library(affy)
library(affyio)
library(annotate)
library(plyr)
biocLite("mogene20sttranscriptcluster.db")
library(mogene20stprobeset.db)
library(mogene20sttranscriptcluster.db)

biocLite("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)

biocLite("pd.hg.u133.plus.2")
library(pd.hg.u133.plus.2)
biocLite("ggplot2")
library(ggplot2)
biocLite("reshape")
biocLite("plot3D")
biocLite("gplots")
biocLite("ggdendro")
biocLite("RColorBrewer")
biocLite("frma")

library("reshape")
library("plot3D")
library("gplots")
library(ggdendro)
library(RColorBrewer)
library(RColorBrewer)
library(frma)

# Read in the CEL files in the directory, then normalize the data
setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM")
data.set1 <- ReadAffy()
ed.raw.set1 <- exprs(data.set1)
samp.set1 <- sampleNames(data.set1)
probes.set1 <- featureNames(data.set1)
ed.normalized.set1<- rma(data.set1)

setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set2_corrected")
data.set2 <- ReadAffy()
ed.raw.set2 <- exprs(data.set2)
samp.set2 <- sampleNames(data.set2)
probes.set2 <- featureNames(data.set2)
ed.normalized.set2<- rma(data.set2)

setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set3_HGSCancer_only")
data.set3 <- ReadAffy()
ed.raw.set3 <- exprs(data.set3)
samp.set3 <- sampleNames(data.set3)
probes.set3 <- featureNames(data.set3)
ed.normalized.set3<- rma(data.set3)

data.set1.normalized<-ed.normalized.set1
data.set2.normalized<-ed.normalized.set2
data.set3.normalized<-ed.normalized.set3

ll(dim=T)

colnames(data.set1.normalized)
colnames(data.set2.normalized)
colnames(data.set3.normalized)

length(colnames(data.set1.normalized))
length(colnames(data.set2.normalized))
length(colnames(data.set3.normalized))

unique(colnames(data.set1.normalized))
length(unique(colnames(data.set1.normalized)))
length(unique(colnames(data.set2.normalized)))
length(unique(colnames(data.set3.normalized)))

data.set1.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/Microarray_Set1_FTE_HGSC_LCM/CASE_DESCRIPTION_SET1.xlsx")
data.set2.sample.info<-read.xls("/home/aiminyan/Downloads/GSE28044_George.xlsx")
data.ER.PR.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx")
head(data.ER.PR.sample.info)

head(data.set1.sample.info)
data.set1.sample.info[,2]
data.set1.normalize
data.set2.normalized

grep(2367,unique(data.set1.sample.info[,2]))

unique(data.set2.sample.info[,2])
head(data.ER.PR.sample.info)

dim(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==1),])
dim(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==2),])
dim(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==3),])
dim(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==4),])

MapSample2CELfile<-function(data.ER.PR.sample.info,data.set1.normalized,sample_type){

cat(dim(data.ER.PR.sample.info),dim(data.set1.normalized),sample_type,"\n")

num.sample1=length(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==1),][,1])
cel_file_index_4_sample<-array()

for(i in 1:num.sample1) {
tmp<-grep(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==sample_type),][i,1],colnames(data.set1.normalized))
cel_file_index_4_sample<-c(cel_file_index_4_sample,tmp)
}
cel_file_index_4_sample<-cel_file_index_4_sample[-1]

colnames(data.set1.normalized)[cel_file_index_4_sample]

re<-data.set1.normalized[,cel_file_index_4_sample]

return(re)

}

#ER-PR-
ER.PR.sample1.from.set1<-MapSample2CELfile(data.ER.PR.sample.info,data.set1.normalized,1)
ER.PR.sample1.from.set2<-MapSample2CELfile(data.ER.PR.sample.info,data.set2.normalized,1)
ER.PR.sample1.from.set3<-MapSample2CELfile(data.ER.PR.sample.info,data.set3.normalized,1)

#ER-PR+
ER.PR.sample2.from.set1<-MapSample2CELfile(data.ER.PR.sample.info,data.set1.normalized,2)
ER.PR.sample2.from.set2<-MapSample2CELfile(data.ER.PR.sample.info,data.set2.normalized,2)
ER.PR.sample2.from.set3<-MapSample2CELfile(data.ER.PR.sample.info,data.set3.normalized,2)

#ER+PR-
ER.PR.sample3.from.set1<-MapSample2CELfile(data.ER.PR.sample.info,data.set1.normalized,3)
ER.PR.sample3.from.set2<-MapSample2CELfile(data.ER.PR.sample.info,data.set2.normalized,3)
ER.PR.sample3.from.set3<-MapSample2CELfile(data.ER.PR.sample.info,data.set3.normalized,3)

#ER+PR+
ER.PR.sample4.from.set1<-MapSample2CELfile(data.ER.PR.sample.info,data.set1.normalized,4)
ER.PR.sample4.from.set2<-MapSample2CELfile(data.ER.PR.sample.info,data.set2.normalized,4)
ER.PR.sample4.from.set3<-MapSample2CELfile(data.ER.PR.sample.info,data.set3.normalized,4)

head(ER.PR.sample1.from.set1)
head(ER.PR.sample1.from.set2)
head(ER.PR.sample1.from.set3)

head(ER.PR.sample3.from.set1)
head(ER.PR.sample3.from.set2)
head(ER.PR.sample3.from.set3)

head(data.set1.sample.info)
head(data.set2.sample.info)
head(data.ER.PR.sample.info)

data.set1.sample.info[,2]

data.set2.sample.info[,8:9:10]

cel.file.all<-c(unique(colnames(data.set1.normalized)),unique(colnames(data.set2.normalized)),unique(colnames(data.set3.normalized)))

length(cel.file.all)
cel.file.all

cel.file.all.2<-rbind(cbind(unique(colnames(data.set1.normalized)),rep(1,length(unique(colnames(data.set1.normalized))))),
cbind(unique(colnames(data.set2.normalized)),rep(2,length(unique(colnames(data.set2.normalized))))),
cbind(unique(colnames(data.set3.normalized)),rep(3,length(unique(colnames(data.set3.normalized))))))

cel.file.cancer.data.set2<-as.character(data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)][,8])
cel.file.cancer.data.set2.2<-as.character(data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)][,1])
cel.file.cancer.data.set2.2.2<-cel.file.cancer.data.set2.2[-which(cel.file.cancer.data.set2.2=="")]

data.set2.cancer.GSM.cel.mapping<-data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)]
data.set2.cancer.GSM.cel.mapping.2<-data.set2.cancer.GSM.cel.mapping[,c(1,8)]
data.set2.cancer.GSM.cel.mapping.3<-data.set2.cancer.GSM.cel.mapping.2[which(as.character(data.set2.cancer.GSM.cel.mapping.2[,1])!=""),]

data.set2.cancer.GSM.cel.mapping.3.2<-gsub(".cel","",as.character(data.set2.cancer.GSM.cel.mapping.3[,2]))
data.set2.cancer.GSM.cel.mapping.3.2.2<-gsub(".CEL","",as.character(data.set2.cancer.GSM.cel.mapping.3.2))
data.set2.cancer.GSM.cel.mapping.3.2.2.2<-gsub("-","_",as.character(data.set2.cancer.GSM.cel.mapping.3.2.2))
data.set2.cancer.GSM.cel.mapping.4<-cbind(data.set2.cancer.GSM.cel.mapping.3,sapply(strsplit(data.set2.cancer.GSM.cel.mapping.3.2.2.2,"_"),"[[",1))

cel.file.cancer.data.set1<-as.character(data.set1.sample.info[which(data.set1.sample.info[,4]=="cancer"),2])

cel.file.cancer.set1.set2<-rbind(cbind(cel.file.cancer.data.set1,rep(1,length(cel.file.cancer.data.set1))),
cbind(cel.file.cancer.data.set2,rep(2,length(cel.file.cancer.data.set2))),
cbind(cel.file.cancer.data.set2.2.2,rep(2,length(cel.file.cancer.data.set2.2.2))))

cel.file.cancer.set1.set2.set3<-rbind(cel.file.cancer.set1.set2,cel.file.all.2[which(cel.file.all.2[,2]==3),])

cel.file.cancer.set1.set2.set3.name<-gsub(".CEL","",cel.file.cancer.set1.set2.set3[,1])
cel.file.cancer.set1.set2.set3.name.1<-gsub(".cel","",cel.file.cancer.set1.set2.set3.name)
cel.file.cancer.set1.set2.set3.name.2<-gsub("-","_",cel.file.cancer.set1.set2.set3.name.1)

data.set123.normalized<-cbind(data.set1.normalized,data.set2.normalized,data.set3.normalized)
dim(data.set123.normalized)
dim(data.set123.normalized[,which(colnames(data.set123.normalized) %in% cel.file.cancer.set1.set2.set3[,1])])
colnames(data.set123.normalized)[grep(54,colnames(data.set123.normalized))]

original.cel.file.names<-colnames(data.set123.normalized)
original.cel.file.names.1<-gsub("X","",original.cel.file.names)
original.cel.file.names.2<-gsub(".CEL","",original.cel.file.names.1)
original.cel.file.names.3<-gsub(".cel","",original.cel.file.names.2)
original.cel.file.names.4<-gsub("\\.","_",original.cel.file.names.3)

index.4.cancer.sample<-which(original.cel.file.names.4 %in% cel.file.cancer.set1.set2.set3.name.2)
original.cel.file.name.with.mapped.names<-cbind(original.cel.file.names[index.4.cancer.sample],original.cel.file.names.4[index.4.cancer.sample])
#index.4.no.cancer.sample.1<-cbind(original.cel.file.names[-index.4.cancer.sample],original.cel.file.names.4[-index.4.cancer.sample])
setdiff(cel.file.cancer.set1.set2.set3.name.2,original.cel.file.names.4[index.4.cancer.sample])
which(original.cel.file.names.4[-index.4.cancer.sample] %in% cel.file.cancer.set1.set2.set3.name.2)

cancer.data.set123.normalized<-data.set123.normalized[,index.4.cancer.sample]
dim(cancer.data.set123.normalized)

cancer.cel.file.name.72<-colnames(cancer.data.set123.normalized)
cancer.cel.file.name.72.1<-gsub("X","",cancer.cel.file.name.72)
cancer.cel.file.name.72.2<-gsub(".CEL","",cancer.cel.file.name.72.1)
cancer.cel.file.name.72.3<-gsub(".cel","",cancer.cel.file.name.72.2)
cancer.cel.file.name.72.4<-gsub("\\.","_",cancer.cel.file.name.72.3)
cancer.cel.file.name.72.5<-cancer.cel.file.name.72.4

data.set2.cancer.GSM.cel.mapping.44<-as.character(data.set2.cancer.GSM.cel.mapping.4[which(data.set2.cancer.GSM.cel.mapping.4[,1] %in% cancer.cel.file.name.72.4),3])
cancer.cel.file.name.72.5[which(cancer.cel.file.name.72.5 %in% data.set2.cancer.GSM.cel.mapping.4[,1])]<-data.set2.cancer.GSM.cel.mapping.44
cancer.cel.file.name.72.6<-cbind(cancer.cel.file.name.72,cancer.cel.file.name.72.5)
cancer.cel.file.name.72.7<-c(sapply(strsplit(cancer.cel.file.name.72.6[1:29,2],"_"),"[[",1),
sapply(strsplit(cancer.cel.file.name.72.6[30:72,2],"_"),"[[",4))
cancer.cel.file.name.72.8<-cbind(cancer.cel.file.name.72.6,cancer.cel.file.name.72.7)

#Classify cancer patients to 4 subtypes

subtype.1<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==1),1])
subtype.2<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==2),1])
subtype.3<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==3),1])
subtype.4<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==4),1])

cancer.data.set123.normalized.st1<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.1)]
cancer.data.set123.normalized.st2<-as.data.frame(cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.2)])
colnames(cancer.data.set123.normalized.st2)<-colnames(cancer.data.set123.normalized)[which(cancer.cel.file.name.72.8[,3] %in% subtype.2)]
cancer.data.set123.normalized.st3<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.3)]
cancer.data.set123.normalized.st4<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.4)]

dim(cancer.data.set123.normalized.st1)
dim(cancer.data.set123.normalized.st2)
dim(cancer.data.set123.normalized.st3)
dim(cancer.data.set123.normalized.st4)

head(cancer.data.set123.normalized.st1)
head(cancer.data.set123.normalized.st2)
head(cancer.data.set123.normalized.st3)
head(cancer.data.set123.normalized.st4)

colnames(cancer.data.set123.normalized.st1)
colnames(cancer.data.set123.normalized.st2)
colnames(cancer.data.set123.normalized.st3)
colnames(cancer.data.set123.normalized.st4)

#n.st<-length(colnames(cancer.data.set123.normalized.st1))

# #Use subtype1,2,3,4
# cel.file.sample.infor<-as.data.frame(rbind(
# cbind(colnames(cancer.data.set123.normalized.st1),rep("st1",length(colnames(cancer.data.set123.normalized.st1)))),
# cbind(colnames(cancer.data.set123.normalized.st2),rep("st2",length(colnames(cancer.data.set123.normalized.st2)))),
# cbind(colnames(cancer.data.set123.normalized.st3),rep("st3",length(colnames(cancer.data.set123.normalized.st3)))),
# cbind(colnames(cancer.data.set123.normalized.st4),rep("st4",length(colnames(cancer.data.set123.normalized.st4))))
# ))
# colnames(cel.file.sample.infor)=c("filename","subtype")

#Use subtype 1,3,4
cel.file.sample.infor.no.2<-as.data.frame(rbind(
  cbind(colnames(cancer.data.set123.normalized.st1),rep("st1",length(colnames(cancer.data.set123.normalized.st1)))),
  cbind(colnames(cancer.data.set123.normalized.st3),rep("st3",length(colnames(cancer.data.set123.normalized.st3)))),
  cbind(colnames(cancer.data.set123.normalized.st4),rep("st4",length(colnames(cancer.data.set123.normalized.st4))))
))

colnames(cel.file.sample.infor.no.2)=c("filename","subtype")
f.st134 <- factor(cel.file.sample.infor.no.2$subtype)
design.st134 <- model.matrix(~0+f.st134)
colnames(design.st134) <- levels(f.st134)

cancer.data.st134<-cbind(cancer.data.set123.normalized.st1,
cancer.data.set123.normalized.st3,
cancer.data.set123.normalized.st4)
head(cancer.data.st134)
GeneSym.all <- getSYMBOL(rownames(cancer.data.st134), "hgu133plus2.db")

ndata<-cancer.data.st134
geneSymbol=GeneSym.all
tempdata.byGSym = data.frame(ndata, Symbol = geneSymbol)

tempdata.byGSym.2<-tempdata.byGSym[-which(is.na(tempdata.byGSym$Symbol)),]

rownames(tempdata.byGSym.2) = NULL
data.byGSym = ddply(tempdata.byGSym.2, c("Symbol"),function(h)
  summary = apply(h,2,max)
)

#class(data.byGSym)

#which(is.na(data.byGSym$Symbol))
rownames(data.byGSym)=data.byGSym$Symbol

data.byGSym.2<-data.byGSym[,-61]

dim(data.byGSym.2)

#Use all probes

fit.st134 <- lmFit(cancer.data.st134, design.st134)

cont.matrix.st134 <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134)
fit2.st134  <- contrasts.fit(fit.st134, cont.matrix.st134)
fit2.st134  <- eBayes(fit2.st134)

str(fit2.st134)

TopTableSt34.all<-topTable(fit2.st134,coef=1,n=54674)

TopTableSt34.100<-topTable(fit2.st134,coef=1,adjust="fdr",n=100)

genenames <- as.character(rownames(TopTableSt34.54675))

length(genenames)

genenames

annotation(ed.normalized.set1)
annotation(ed.normalized.set2)
annotation(ed.normalized.set3)
library("hgu133plus2.db")

map <- getpAnnMa("ENTREZID", "hgu133plus2", load=TRUE, type=c("env", "db"))
class(map)
ll <- getEG(genenames,"hgu133plus2.db")
GeneSym <- getSYMBOL(genenames, "hgu133plus2.db")

Probe.gene.sym<-(cbind(genenames,GeneSym))

Probe.gene.sym

length(which(is.na(Probe.gene.sym[,2])))


tab <- data.frame(GeneSym, TopTableSt34.100)
tab <- data.frame(rownames(tab), tab)

colnames(tab)[1] <- c("Probe ID")
ll <- list(ll)
htmlpage(ll, filename="/media/H_driver/2015/Sophia/St34-4.html", title="HTML report",
         othernames=tab, table.head=c("Locus ID",colnames(tab)), table.center=TRUE, digits=6)

colnames(fit2.st134)
topTable(fit2.st134,coef=1)
topTable(fit2.st134,coef=1,adjust="fdr")

topTable(fit2.st134,coef=2)
topTable(fit2.st134,coef=2,adjust="fdr")

results.st134 <- decideTests(fit2.st134,p.value=0.99)

summary(results.st134)
vennDiagram(results.st134)

fit2.st134$genes

unigeneTopTableSt34 <- topTable(fit2.st134,coef=1,n=20,genelist=genelist)
unigeneTopTableSt41 <- topTable(fit2.st134,coef=2,n=20,genelist=genelist)

library(xtable)
xtableUnigeneSt34 <- xtable(unigeneTopTableSt34,display=c("s","s","s","s","g","g","g","e","e","g","g"))
xtableUnigeneSt41 <- xtable(unigeneTopTableSt41,display=c("s","s","s","s","g","g","g","e","e","g","g"))

cat(file="/media/H_driver/2015/Sophia/St34.html","<html>\n<body>")
print.xtable(xtableUnigeneSt34,type="html",file="/media/H_driver/2015/Sophia/St34.html",append=TRUE)
cat(file="/media/H_driver/2015/Sophia/St34.html","</body>\n</html>",append=TRUE)

cat(file="/media/H_driver/2015/Sophia/St41.html","<html>\n<body>")
print.xtable(xtableUnigeneSt41,type="html",file="/media/H_driver/2015/Sophia/St41.html",append=TRUE)
cat(file="/media/H_driver/2015/Sophia/St41.html","</body>\n</html>",append=TRUE)


#Use gene name
temp = apply(data.byGSym.2,2,as.numeric)
rownames(temp) = rownames(data.byGSym.2)
data.byGSym.2 = temp


fit.st134.gene <- lmFit(data.byGSym.2, design.st134)

cont.matrix.st134 <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134)
fit2.st134  <- contrasts.fit(fit.st134.gene, cont.matrix.st134)
fit2.st134  <- eBayes(fit2.st134)

str(fit2.st134)

TopTableSt34.gene<-topTable(fit2.st134,coef=1,n=54674)
dim(TopTableSt34.gene)
head(TopTableSt34.gene)
hist(TopTableSt34.gene[,4])

TopTableSt41.gene<-topTable(fit2.st134,coef=2,n=54674)
dim(TopTableSt41.gene)
head(TopTableSt41.gene)
hist(TopTableSt41.gene[,4])



SampleType = factor(gsub("(FTY).[1-3]","\\1",colnames(data.byGSym)))

SampleType<-cel.file.sample.infor.no.2$subtype

heatmap_wPCA(data.byGSym.2, "heatmap_allsample.pdf","PCA_allsample.pdf","/media/H_driver/2015/Sophia/", SampleType)

# clustering plot

#fit#Differential gene expression analysis for st3 vs st4
#probe<-rownames(cancer.data.set123.normalized.st1)
#geneSymbol<-getSYMBOL(probe, "hugene10sttranscriptcluster.db")
#Differential gene expression analysis for st2 vs st1
#Differential gene expression analysis for st4 vs st1
save.image(file="/media/H_driver/Aimin_project/Sophia/Data_set_3.RData")




biocLite("frmaExampleData")
biocLite("hgu133afrmavecs")
library(frmaExampleData)
library(hgu133afrmavecs)

data(AffyBatchExample)
object <- frma(AffyBatchExample)
#' Title
#'
#' @param reorder.column.ed.set123.frma.cancer
#' @param GeneSym.all
#' @param data.set123.raw.with.set.label.and.type.cancer
#' @param output_file
#' @param output_dir
#'
#' @return
#' @export
#'
#' @examples


load("/media/H_driver/Aimin_project/Sophia/Data_4_heatmap_PCA.RData")
library("plot3D")
library("gplots")
library(ggdendro)
library(RColorBrewer)
library(RColorBrewer)
library(ggplot2)
library(plyr)
heatmap_wPCA(ed.set123.frma.cancer.reorder,"heatmap_allsample_frma_all_probes_2.pdf","PCA_allsample_frma_all_probes_.pdf","/media/H_driver/2015/Sophia/Results/",check.match$subtype)x = readHTMLTable('http://www.disastercenter.com/crime/iacrime.htm')
colnames(data.set123.raw)

data.set123.raw.with.set.label<-rbind(cbind(colnames(data.set1),rep("DS1",length(colnames(data.set1)))),
cbind(colnames(data.set2),rep("DS2",length(colnames(data.set2)))),
cbind(colnames(data.set3),rep("DS3",length(colnames(data.set3)))))

dim(data.set123.raw.with.set.label)
dim(ed.set123.frma)
colnames(ed.set123.frma)

reorder.column.ed.set123.frma<-ed.set123.frma[,match(data.set123.raw.with.set.label[,1],colnames(ed.set123.frma))]

colnames(reorder.column.ed.set123.frma)
#data.set123.raw.with.set.label[,1]

#colnames(data.byGSym.2)
#SampleType

ndata<-reorder.column.ed.set123.frma

geneSymbol=GeneSym.all
tempdata.byGSym =data.frame(ndata, Symbol = geneSymbol)

colnames(tempdata.byGSym)

tempdata.byGSym.2<-tempdata.byGSym[-which(is.na(tempdata.byGSym$Symbol)),]

# test.sample<-ed.set123.frma.cancer[which(GeneSym.all=="ERBB4"),which(colnames(ed.set123.frma.cancer) %in% cel.file.sample.infor.no.8[which(cel.file.sample.infor.no.8[,5]=="st4"),6])]
#
# test.sample.2<-data.frame(test.sample,Symbol="ERBB4")
# test.sample.3 = test.sample.2
# test.sample.3[,1:5] = apply(test.sample.3[,1:5],2,as.numeric)

data.byGSym = ddply(tempdata.byGSym.2, c("Symbol"),function(h)
  summary = apply(h[,1:(dim(tempdata.byGSym.2)[2]-1)],2,max)
)


rownames(data.byGSym)=data.byGSym$Symbol
colnames(data.byGSym)

data.byGSym.2<-data.byGSym[,-1]
temp = apply(data.byGSym.2,2,as.numeric)
rownames(temp) = rownames(data.byGSym.2)
colnames(temp) = colnames(ndata)

data.byGSym.2 = temp
colnames(data.byGSym.2)
#data.byGSym.2 = ed.set123.frma.cancer

#SampleType<-data.set123.raw.with.set.label[,2]

Draw_PCA(data.byGSym.2,"PCA_all_samples_frma_3_data_sets.pdf","/media/H_driver/2015/Sophia/Results/",data.set123.raw.with.set.label[,2])

#heatmap_wPCA(data.byGSym.2, "heatmap_allsample_frma_all_probes.pdf","PCA_allsample_frma_all_probes.pdf","/media/H_driver/2015/Sophia/Results/", SampleType)


#Based on cancer samples only in data set 1,2,3


data.set123.raw.with.set.label.and.type<-cbind(data.set123.raw.with.set.label,rep("Normal",dim(data.set123.raw.with.set.label)[1]))

data.set123.raw.with.set.label.and.type[which(data.set123.raw.with.set.label.and.type[,1] %in% colnames(ed.set123.frma.cancer)),3]<-"Cancer"

data.set123.raw.with.set.label.and.type

data.set123.raw.with.set.label.and.type.cancer<-data.set123.raw.with.set.label[which(data.set123.raw.with.set.label.and.type[,3]=="Cancer"),]

reorder.column.ed.set123.frma.cancer<-ed.set123.frma.cancer[,match(data.set123.raw.with.set.label.and.type.cancer[,1],colnames(ed.set123.frma.cancer))]

dim(reorder.column.ed.set123.frma.cancer)

ndata<-reorder.column.ed.set123.frma.cancer

geneSymbol=GeneSym.all
tempdata.byGSym =data.frame(ndata, Symbol = geneSymbol)

colnames(tempdata.byGSym)

tempdata.byGSym.2<-tempdata.byGSym[-which(is.na(tempdata.byGSym$Symbol)),]

# test.sample<-ed.set123.frma.cancer[which(GeneSym.all=="ERBB4"),which(colnames(ed.set123.frma.cancer) %in% cel.file.sample.infor.no.8[which(cel.file.sample.infor.no.8[,5]=="st4"),6])]
#
# test.sample.2<-data.frame(test.sample,Symbol="ERBB4")
# test.sample.3 = test.sample.2
# test.sample.3[,1:5] = apply(test.sample.3[,1:5],2,as.numeric)

data.byGSym = ddply(tempdata.byGSym.2, c("Symbol"),function(h)
  summary = apply(h[,1:(dim(tempdata.byGSym.2)[2]-1)],2,max)
)


rownames(data.byGSym)=data.byGSym$Symbol
colnames(data.byGSym)

data.byGSym.2<-data.byGSym[,-1]
temp = apply(data.byGSym.2,2,as.numeric)
rownames(temp) = rownames(data.byGSym.2)
colnames(temp) = colnames(ndata)

data.byGSym.2 = temp
colnames(data.byGSym.2)


Draw_PCA(data.byGSym.2,"PCA_cancer_samples_only_frma_3_data_sets.pdf","/media/H_driver/2015/Sophia/Results/",data.set123.raw.with.set.label.and.type.cancer[,2])

#colnames(reorder.column.ed.set123.frma.cancer)

#colnames(ed.set123.frma.cancer)

#Use cancer samples in data set 3 only
data.set3.raw.with.set.label.and.type.cancer<-data.set123.raw.with.set.label.and.type[which(data.set123.raw.with.set.label.and.type[,3]=="Cancer"
                                                                                          &data.set123.raw.with.set.label.and.type[,2]=="DS3"),]

savehistory(file="/media/H_driver/Aimin_project/Sophia/Redo_analysis_PCA_based_on_Data_set.Rhistory")
save.image(file="/media/H_driver/Aimin_project/Sophia/Redo_analysis_PCA_based_on_Data_set.RData")
