# Output DE analysis results
#' Title
#'
#' @param DE4cut_off_5
#'
#' @return
#' @export
#'
#' @examples
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
