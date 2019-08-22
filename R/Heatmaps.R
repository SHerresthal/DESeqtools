#' Plot a heatmap
#'
#' This function provides a wrapper for the \link[pheatmap:pheatmap]{pheatmap}.
#' It can plot a subset of genes, either annotated as EnsemblIDs or Gene Symbol
#'
#' @param input a dataframe for the read counts, norm_anno as the default
#' @param geneset A character vector containing the genes to plot or "all" for plotting all genes of the dataset (hardcoded to norm_anno). By default, all genes are plotted.
#' @param title The title of the plot
#' @param smp_table A sample table for processing the heatmap, by default smp_table
#' @param gene_type A character from gene_annotation$gene_type defineing the genes set to use. "all" is set by default and will use all genes.
#' @param keyType Annotation of the genes in geneset. Either "Ensembl" or "Symbol"
#' @param show_rownames Shall the gene names be shown in the heatmap? Default is FALSE.
#' @param cluster_cols Shall "test" the colums be clusteres? Default is FALSE
#' @return Returns a heatmap
#' @export

plotHeatmap <- function(input=norm_anno,
                        smp_table=sample_table,
                        geneset,
                        title="",
                        keyType = "Ensembl",
                        show_rownames = FALSE,
                        cluster_cols = FALSE,
                        gene_type="all",
                        plot_anno=plot_annotation){
  
  if(gene_type=="all"){
    if(geneset[1] =="all"){
      input <- input
    }else{#
      if(keyType == "Ensembl"){
        input <- input[input$GENEID %in% geneset,]
      } else if(keyType == "Symbol"){
        input <- input[input$SYMBOL %in% geneset,]
      } else{
        print("Wrong keyType. Choose Ensembl or Symbol!")
      }
    }
    rownames(input) <- paste(input$GENEID, ":", input$SYMBOL, sep="")
    input <- input[,colnames(input) %in% smp_table[["ID"]]]
    
    input_scale <- t(scale(t(input)))
    input_scale <- input_scale[,order(smp_table[[plot_order]], decreasing = FALSE)]
    
    pheatmap(input_scale,
             main=title,
             show_rownames=show_rownames,
             show_colnames=TRUE,
             cluster_cols = cluster_cols, 
             fontsize = 7,
             annotation_col = plot_anno,
             annotation_colors = ann_colors,
             breaks = scaleColors(data = input_scale, maxvalue = 2)[["breaks"]], 
             color = scaleColors(data = input_scale, maxvalue = 2)[["color"]])
  }else{
    if(geneset[1] =="all"){
      input <- input
    }else{#
      if(keyType == "Ensembl"){
        input <- input[input$GENEID %in% geneset,]
      } else if(keyType == "Symbol"){
        input <- input[input$SYMBOL %in% geneset,]
      } else{
        print("Wrong keyType. Choose Ensembl or Symbol!")
      }
    }
    rownames(input) <- paste(input$GENEID, ":", input$SYMBOL, sep="")
    input<-input[input[["GENETYPE"]]==gene_type,]
    input <- input[,colnames(input) %in% smp_table[["ID"]]]
    
    input_scale <- t(scale(t(input)))
    input_scale <- input_scale[,order(smp_table[[plot_order]], decreasing = FALSE)]
    
    pheatmap(input_scale,
             main=title,
             show_rownames=show_rownames,
             show_colnames=TRUE,
             cluster_cols = cluster_cols, 
             fontsize = 7,
             annotation_col = plot_anno,
             annotation_colors = ann_colors,
             breaks = scaleColors(data = input_scale, maxvalue = 2)[["breaks"]], 
             color = scaleColors(data = input_scale, maxvalue = 2)[["color"]]) 
  }
  
  
  
}



#' Heatmap of genes of specified GO, KEGG or HALLMARK gene sets
#'
#' This function provides a wrapper for the \link[pheatmap:pheatmap]{pheatmap}.
#' It can plot a subset of genes from a certain gene set based on the GO, KEGG or HALLMARK database.
#'
#' @param cat Category of the geneset. Either KEGG, GO or HALLMARK
#' @param term Term / ID of the geneset. Either a GO term, a KEGG pathway or a set of Hallmark genes
#' @param organism either "mouse" or "human"
#' @param gene_type A character from gene_annotation$gene_type defineing the genes set to use. "all" is set by default and will use all genes.
#' @param show_rownames Shall the gene names be shown in the heatmap? Default is FALSE.
#' @param cluster_cols Shall the colums be clustered? Default is TRUE
#' @return Returns a heatmap
#' @export

plotGeneSetHeatmap <- function(input=norm_anno,
                               smp_table=sample_table,
                               plot_anno=plot_annotation,
                             
                               cat,
                               term,
                               organism,
                               show_rownames =TRUE, 
                               cluster_cols = FALSE,
                               gene_type="all"){
  
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
      genes <- getLDS(attributes = c("entrezgene_id"), 
                      filters = "entrezgene_id", 
                      values = genes, 
                      mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl"), 
                      attributesL = c("mgi_symbol"), 
                      martL = useMart("ensembl", dataset = "mmusculus_gene_ensembl"), 
                      uniqueRows=T)[,2]
    }
  }
  
  plotHeatmap(input=input, 
              smp_table=smp_table,
              geneset = genes,
              keyType = "Symbol",
              title = paste("Heatmap of present genes annotated to: ",
                            term, sep=""),
              show_rownames = show_rownames,
           
              cluster_cols = cluster_cols,
              plot_anno=plot_anno,
              gene_type=gene_type)
}

#' Function to plot a heatmap of genes responsible for gene set enrichment
#' @export
plotGSEAHeatmap<-function(input=norm_anno,
                          smp_table = sample_table,
                          plot_anno = plot_annotation, 
                          GSEA_result,
                          GeneSet, 
                          term,
                          regulation,
                          show_rownames = TRUE,
                          cluster_cols = F,
                          gene_type="all"){
  
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
  
  plotHeatmap(input=input,
              smp_table = smp_table,
              plot_anno = plot_anno,
              geneset = gene.list,
              keyType = "Symbol",
              title = paste("Heatmap of genes responsible for enrichment of term:",
                            term,", in ",deparse(substitute(GSEA_result)),sep=""),
              show_rownames = show_rownames,
              cluster_cols = cluster_cols,
              gene_type=gene_type)
}

#' Function to plot heatmaps of DE genes based on pheatmap
#' @export
plotDEHeatmap <- function(input=norm_anno,
                          smp_table=sample_table,
                          plot_anno=plot_annotation,
                          comparison,
                          factor,
                          gene_anno=gene_annotation,
                          conditions="all",
                          gene_type="all",
                          show_rownames = FALSE,
                          cluster_cols = FALSE){
  
  geneset <- DEresults[[comparison]]@results[DEresults[[comparison]]@results$regulation %in% c("up","down"),"GENEID"]
  
  input <- input[input$GENEID %in% geneset,]
  
  if(conditions[1] == "all"){
    input <- input[,colnames(input) %in% smp_table$ID]
    input_scale <- t(scale(t(input)))
  } else {
    input <- input[,colnames(input) %in% smp_table[as.vector(smp_table[[factor]]) %in% conditions,]$ID,]
    input_scale <- t(scale(t(input)))
    smp_table<-subset(smp_table,smp_table[[factor]] %in% conditions)
  }
  
  input_scale<-as.data.frame(input_scale)
  input_scale$GENEID <- rownames(input_scale)
  gene_anno <- gene_anno[match(rownames(input_scale), gene_anno$GENEID),]
  input_scale <- merge(input_scale,
                        gene_anno,
                        by = "GENEID")
  rownames(input_scale) <- input_scale$GENEID
  
  title=paste("Heatmap of significant DE genes in: ",comparison,sep="")
  
  
  
  plotHeatmap(input=input_scale, 
              smp_table=smp_table,
              geneset = geneset,
              title = title,
              keyType = "Ensembl",#Deniz: cannot be "Symbol", no matter what!!!!!!
              show_rownames = show_rownames,
              cluster_cols = cluster_cols,
              plot_anno=plot_anno,
              gene_type=gene_type)
}
 
