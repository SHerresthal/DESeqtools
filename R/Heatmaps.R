#' Plot a heatmap
#'
#' This function provides a wrapper for the \link[pheatmap:pheatmap]{pheatmap}.
#' It can plot a subset of genes, either annotated as EnsemblIDs or Gene Symbol
#'
#' @param geneset A character vector containing the genes to plot or "all" for plotting all genes of the dataset (hardcoded to norm_anno). By default, all genes are plotted.
#' @param title The title of the plot
#' @param keyType Annotation of the genes in geneset. Either "Ensembl" or "Symbol"
#' @param show_rownames Shall the gene names be shown in the heatmap? Default is FALSE.
#' @param cluster_cols Shall the colums be clustered? Default is FALSE
#' @param column_annotation column annotation of the heatmap
#' @return Returns a heatmap
#' @export

plotHeatmap <- function (input = norm_anno,
                         geneset = "all",
                         title = "",
                         keyType = "Ensembl",
                         show_rownames = FALSE,
                         cluster_cols = FALSE,
                         column_annotation = plot_annotation,
                         sample_annotation = sample_table)
{
  if (geneset[1] != "all") {
    if (keyType == "Ensembl") {
      input <- input[input$GENEID %in% geneset, ]
    }
    else if (keyType == "Symbol") {
      input <- input[input$SYMBOL %in% geneset, ]
    }
    else {
      print("Wrong keyType. Choose Ensembl or Symbol!")
    }
  }
  rownames(input) <- paste(input$GENEID, ":", input$SYMBOL,
                           sep = "")
  input <- input[, colnames(input) %in% sample_annotation$ID]
  input_scale <- t(scale(t(input)))
  input_scale <- input_scale[, order(sample_annotation[[plot_order]],
                                     decreasing = FALSE)]
  pheatmap(input_scale, main = title,
           show_rownames = show_rownames,
           show_colnames = TRUE,
           cluster_cols = cluster_cols,
           fontsize = 7,
           annotation_col = column_annotation,
           annotation_colors = ann_colors,
           breaks = scaleColors(data = input_scale, maxvalue = 2)[["breaks"]],
           color = scaleColors(data = input_scale, maxvalue = 2)[["color"]])
}



#' Heatmap of genes of specified GO, KEGG or HALLMARK gene sets
#'
#' This function provides a wrapper for the \link[pheatmap:pheatmap]{pheatmap}.
#' It can plot a subset of genes from a certain gene set based on the GO, KEGG or HALLMARK database.
#'
#'#############
#' @param input dataset to be plotted, default is norm_anno
#' @param cat Category of the geneset. Either KEGG, GO or HALLMARK
#' @param term Term / ID of the geneset. Either a GO term, a KEGG pathway or a set of Hallmark genes
#' @param organism either "mouse" or "human"
#' @param show_rownames Shall the gene names be shown in the heatmap? Default is FALSE.
#' @param cluster_cols Shall the colums be clustered? Default is TRUE
#' @param column_annotation. column annotation of the heatmap
#' @return Returns a heatmap
#' @export

plotGeneSetHeatmap <- function(input. = norm_anno,
                               sample_annotation. = sample_table,
                               cat,
                               term,
                               organism,
                               show_rownames =TRUE,
                               cluster_cols = FALSE,
                               column_annotation. =  plot_annotation){
  if(organism == "mouse"){
    GO <- GO_mm
    KEGG <- KEGG_mm
  } else if(organism == "human"){
    GO <- GO_hs
    KEGG <- KEGG_hs
  } else (stop("Wrong organism specified!"))

  xterm <- paste("^", term, "$", sep="")
  if(cat=="GO"){
    genes <- unique(GO[grep(xterm,GO$TERM),"SYMBOL"])
  }
  if(cat=="KEGG"){
    genes <- unique(KEGG[grep(xterm,KEGG$PATHWAY),"SYMBOL"])
  }
  if(cat=="HALLMARK"){
    genes <- unique(hallmark_genes[grep(xterm,hallmark_genes$ont),"gene"])
    if(organism == "mouse"){
      genes <- getLDS(attributes = c("entrezgene"),
                      filters = "entrezgene",
                      values = genes,
                      mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
                      attributesL = c("mgi_symbol"),
                      martL = useMart("ensembl", dataset = "mmusculus_gene_ensembl"),
                      uniqueRows=T)[,2]
    }
  }
  plotHeatmap(input = input.,
              sample_annotation = sample_annotation.,
              geneset = genes,
              keyType = "Symbol",
              title = paste("Heatmap of present genes annotated to: ",term, sep=""),
              show_rownames = show_rownames,
              cluster_cols = cluster_cols,
              column_annotation = column_annotation.)
}

#' Function to plot a heatmap of genes responsible for gene set enrichment
#' @export
plotGSEAHeatmap<-function(GSEA_result,
                          GeneSet,
                          term,
                          regulation){
  xterm <- paste("^", term, "$", sep="")
  tmp <- GSEA_result[grep(xterm,GSEA_result$Description),]
  gene.list <- unique(unlist(strsplit(tmp$geneID, split = "/")))

  if(GeneSet == "KEGG"){
    gene.list <- bitr(gene.list,
                      fromType = "ENTREZID",
                      toType="SYMBOL",
                      OrgDb="org.Mm.eg.db")[,2]
  }

  if(GeneSet == "HALLMARK" | GeneSet == "ImmunoSignatures" | GeneSet == "Motifs"){
    gene.list <- getLDS(attributes = c("hgnc_symbol"),
                        filters = "hgnc_symbol",
                        values = gene.list,
                        mart = human,
                        attributesL = c("mgi_symbol"),
                        martL = mouse,
                        uniqueRows=T)[,2]
  }

  plotHeatmap(geneset = gene.list,
              keyType = "Symbol",
              title = paste("Heatmap of genes responsible for enrichment of term: ",term,", in ",deparse(substitute(GSEA_result)),sep=""),
              show_rownames = TRUE,
              cluster_cols = FALSE)
}

#' Function to plot heatmaps of DE genes based on pheatmap
#' @export
plotDEHeatmap <- function(comparison,
                          factor,
                          conditions="all",
                          show_rownames = FALSE,
                          cluster_cols = FALSE){

  geneset <- DEresults[[comparison]]@results[DEresults[[comparison]]@results$regulation %in% c("up","down"),"GENEID"]

  input <- norm_anno[norm_anno$GENEID %in% geneset,]
  rownames(input) <- paste(input$GENEID, ":", input$SYMBOL, sep="")

  if(conditions[1] == "all"){
    input <- input[,colnames(input) %in% sample_table$ID]
    input_scale <- t(scale(t(input)))
    input_scale <- input_scale[,order(sample_table[[plot_order]], decreasing = FALSE)]
  } else {
    input <- input[,colnames(input) %in% sample_table[as.vector(sample_table[[factor]]) %in% conditions,]$ID,]
    input_scale <- t(scale(t(input)))
  }

  pheatmap(input_scale,
           main=paste("Heatmap of significant DE genes in: ",comparison,sep=""),
           show_rownames=show_rownames,
           show_colnames=TRUE,
           cluster_cols = cluster_cols,
           fontsize = 7,
           annotation_col = plot_annotation,
           annotation_colors = ann_colors,
           breaks = scaleColors(data = input_scale, maxvalue = 2)[["breaks"]],
           color = scaleColors(data = input_scale, maxvalue = 2)[["color"]])
}

