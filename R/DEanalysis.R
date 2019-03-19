#' Specify structure of DESeq2_analysis_object
#'
#' @export
setClass(Class = "DESeq2_analysis_object",
         slots = c(results="data.frame", DE_genes="list", Number_DE_genes="list"))


#' DESeq2 differential testing
#'
#' Wrapper Function to perform DESeq2 differential testing
#' @export
DEAnalysis <- function(condition,
                       alpha = 0.05,
                       lfcThreshold = 0,
                       sigFC = 2,
                       multiple_testing = "IHW",
                       independentFiltering="TRUE",
                       shrinkage = TRUE,
                       shrinkType = "normal"){
  # create results_list
  results_list <- list()
  # print parameters
  results_list$parameters <-list(multiple_testing = multiple_testing,
                                 p_value_threshold = alpha,
                                 log2_FC_threshold = lfcThreshold,
                                 shrinkage = shrinkage,
                                 shrinkage_type = shrinkType)
  # Run results() function on comparisons defined in comparison table
  for (i in 1:nrow(comparison_table)){
    # create DE_object
    DE_object <- new(Class = "DESeq2_analysis_object")
    # IHW
    if (multiple_testing=="IHW") {
      res_deseq_lfc <- results(dds,
                               contrast = c(condition,
                                            paste(comparison_table$comparison[i]),
                                            paste(comparison_table$control[i])),
                               lfcThreshold = lfcThreshold,
                               alpha = alpha,
                               filterFun = ihw,
                               altHypothesis = "greaterAbs")
      # Independent Filtering
    }else {
      res_deseq_lfc <- results(dds,
                               contrast = c(condition,
                                            paste(comparison_table$comparison[i]),
                                            paste(comparison_table$control[i])),
                               lfcThreshold = lfcThreshold,
                               alpha = alpha,
                               independentFiltering = independentFiltering,
                               altHypothesis = "greaterAbs",
                               pAdjustMethod= multiple_testing)
    }
    if(shrinkage == TRUE){
      res_deseq_lfc <- lfcShrink(dds,
                                 contrast = c(condition,
                                              paste(comparison_table$comparison[i]),
                                              paste(comparison_table$control[i])),
                                 res=res_deseq_lfc,
                                 type = shrinkType)
    }
    res_deseq_lfc <- as.data.frame(res_deseq_lfc)
    # indicate significant DE genes
    res_deseq_lfc$regulation <- ifelse(!is.na(res_deseq_lfc$padj)&
                                         res_deseq_lfc$padj <= alpha&
                                         res_deseq_lfc$log2FoldChange > log(sigFC,2),
                                       "up",
                                       ifelse(!is.na(res_deseq_lfc$padj)&
                                                res_deseq_lfc$padj <= alpha&
                                                res_deseq_lfc$log2FoldChange < -log(sigFC,2),
                                              "down",
                                              "n.s."))
    # add gene annotation to results table
    res_deseq_lfc$GENEID <- row.names(res_deseq_lfc) # ensembl-IDs as row names
    res_deseq_lfc <- merge(res_deseq_lfc,
                           norm_anno[,c("GENEID",
                                        "SYMBOL",
                                        "GENETYPE",
                                        "DESCRIPTION",
                                        "CHR")],
                           by = "GENEID")
    row.names(res_deseq_lfc) <- res_deseq_lfc$GENEID
    res_deseq_lfc$comparison<-paste(comparison_table$comparison[i]," vs ",comparison_table$control[i],
                                    sep="")
    # re-order results table
    if (multiple_testing=="IHW") {
      res_deseq_lfc<-res_deseq_lfc[,c("GENEID",
                                      "SYMBOL",
                                      "GENETYPE",
                                      "DESCRIPTION",
                                      "CHR",
                                      "comparison",
                                      "regulation",
                                      "baseMean",
                                      "log2FoldChange",
                                      "lfcSE",
                                      "stat",
                                      "pvalue",
                                      "padj",
                                      "weight")]
    }else{
      res_deseq_lfc<-res_deseq_lfc[,c("GENEID",
                                      "SYMBOL",
                                      "GENETYPE",
                                      "DESCRIPTION",
                                      "CHR",
                                      "comparison",
                                      "regulation",
                                      "baseMean",
                                      "log2FoldChange",
                                      "lfcSE",
                                      "stat",
                                      "pvalue",
                                      "padj")]
    }
    # print result table
    DE_object@results <- res_deseq_lfc
    # print DE genes in seperate tables
    DE_object@DE_genes <- list(up_regulated_Genes = res_deseq_lfc[res_deseq_lfc$regulation =="up",],
                               down_regulated_Genes= res_deseq_lfc[res_deseq_lfc$regulation =="down",])
    # print the numbers of DE genes
    DE_object@Number_DE_genes <- list(up_regulated_Genes = nrow(DE_object@DE_genes$up_regulated_Genes),
                                      down_regulated_Genes= nrow(DE_object@DE_genes$down_regulated_Genes))
    # write DE_object into results_list
    results_list[[paste(comparison_table$comparison[i], "vs", comparison_table$control[i], sep="_")]] <- DE_object
  }
  return(results_list)
}


#'  Union of DEgenes
#'
#'  Gives the intersection of DE gene lists (comparisons)
#'
#'  @export

uDEG <- function(comparisons){
  uDEGs <- NULL
  tmp <- DEresults[names(DEresults) %in% comparisons]
  for(i in 1:length(comparisons)){
    DEGs <- as.data.frame(tmp[[i]]@results[tmp[[i]]@results$regulation %in% c("up","down"),])
    uDEGs <- unique(c(uDEGs, DEGs$GENEID))
  }
  uDEGs
}

