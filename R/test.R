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
