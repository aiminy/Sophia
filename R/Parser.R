#' Title
#'
#' @return
#' @export
#'
#' @examples
#' Parser("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM",
#' "/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set2_corrected",
#' "/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set3_HGSCancer_only",
#' "/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/Microarray_Set1_FTE_HGSC_LCM/CASE_DESCRIPTION_SET1.xlsx",
#' "/home/aiminyan/Downloads/GSE28044_George.xlsx",
#' "/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx"
#' )
#'
#'
Parser <- function(dir_set1,dir_set2,dir_set3,dir_sample1,dir_sample2,dir_sample3) {
  # Read in the CEL files in the directory, then normalize the data
  setwd(dir_set1)
  data.set1 <- ReadAffy()
  ed.raw.set1 <- exprs(data.set1)
  samp.set1 <- sampleNames(data.set1)
  probes.set1 <- featureNames(data.set1)
  ed.normalized.set1<- rma(data.set1)

  setwd(dir_set2)
  data.set2 <- ReadAffy()
  ed.raw.set2 <- exprs(data.set2)
  samp.set2 <- sampleNames(data.set2)
  probes.set2 <- featureNames(data.set2)
  ed.normalized.set2<- rma(data.set2)

  setwd(dir_set3)
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

  data.set1.sample.info<-read.xls(dir_sample1)
  data.set2.sample.info<-read.xls(dir_sample2)
  data.ER.PR.sample.info<-read.xls(dir_sample3)

  re<-list(set1=data.set1.normalized,set2=data.set2.normalized,set3=data.set3.normalized,sample1=data.set1.sample.info,sample2=data.set2.sample.info,sample3=data.ER.PR.sample.info)

  return(re)

  }
