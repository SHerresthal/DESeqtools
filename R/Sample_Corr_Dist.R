#' Correlation & Distance
#'
#' This function provides the sample-to-sample correlation .
#' 
#'
#' @param input A dataframe containing the read counts of the genes for each sample. By default, "norm_anno" normalized read counts is used.
#' @param anno A sample table containing info for the samples
#' @param condition A character from sample table df defining the name of column in the table to use. "condition" is set by default.
#' @return Returns a dataframe
#' @export

corr_function<-function(sampleCor = cd_input,
                        gene_anno=gene_annotation,
                        plot_anno=plot_annotation,
                        title=title,
                        gene_type="all",
                        cluster_rows = F,
                        cluster_cols = F,
                        mean=F){
  
  
  if(mean==T){
    sampleCor$GENEID <- row.names(sampleCor)
    gene_anno <- gene_anno[match(rownames(sampleCor), gene_anno$GENEID),]
    sampleCor <- merge(sampleCor,
                       gene_anno,
                       by = "GENEID")
    
    rownames(sampleCor) <- sampleCor$GENEID
    sampleCor<-mean_function(input=sampleCor,
                             anno=sample_table,
                             condition="condition")
    
    sampleCor <- sampleCor[,colnames(sampleCor) %in% sample_table[["condition"]]]
    
  }else{
    sampleCor<-sampleCor
  }
  
  if(gene_type=="all"){
    
    if(mean==T){
      sampleCor <- sampleCor[,colnames(sampleCor) %in% sample_table[["condition"]]]
      sampleCor <- as.matrix(cor(sampleCor, use="all.obs", method="pearson"))
      rownames(sampleCor)<- unique(sample_table$condition)
      colnames(sampleCor)<- unique(sample_table$condition)
      
    }else{
      sampleCor <- sampleCor[,colnames(sampleCor) %in% sample_table[["ID"]]]
      sampleCor <- as.matrix(cor(sampleCor, use="all.obs", method="pearson"))
      rownames(sampleCor)<- sample_table$ID
      colnames(sampleCor)<- sample_table$ID
      
    }
    
    pheatmap(sampleCor,
             main="Sample Correlation based on variance-stabilized counts",
             annotation_row = plot_anno,
             annotation_col = plot_anno,
             annotation_colors = ann_colors,
             cluster_rows = cluster_rows,
             cluster_cols = cluster_cols,
             fontsize = 8)
    
  }else{
    #Deniz:filtering if gene_type is not "all"
    sampleCor$GENEID <- row.names(sampleCor)
    gene_anno <- gene_anno[match(rownames(sampleCor), gene_anno$GENEID),]
    sampleCor <- merge(sampleCor,
                       gene_anno,
                       by = "GENEID")
    
    rownames(sampleCor) <- sampleCor$GENEID
    sampleCor<-sampleCor[sampleCor[["GENETYPE"]]==gene_type,]
    
    if(mean==T){
      sampleCor <- sampleCor[,colnames(sampleCor) %in% sample_table[["condition"]]]
      sampleCor <- as.matrix(cor(sampleCor, use="all.obs", method="pearson"))
      rownames(sampleCor)<- unique(sample_table$condition)
      colnames(sampleCor)<- unique(sample_table$condition)
      
      
    }else{
      sampleCor <- sampleCor[,colnames(sampleCor) %in% sample_table[["ID"]]]
      sampleCor <- as.matrix(cor(sampleCor, use="all.obs", method="pearson"))
      rownames(sampleCor)<- sample_table$ID
      colnames(sampleCor)<- sample_table$ID
      
      
    }
    
    pheatmap(sampleCor,
             main=title,
             annotation_row = plot_anno,
             annotation_col = plot_anno,
             annotation_colors = ann_colors,
             cluster_rows = cluster_rows,
             cluster_cols = cluster_cols,
             fontsize = 8)
  }
}

#' This function provides the sample-to-sample distances .
#' 
#'
#' @param sampleDist input
#' @param sample_table A sample table containing info for the samples
#' @param gene_anno 
#' @return Returns a dataframe
#' @export
#' 
 
dist_function<-function(sampleDist = cd_input,
                        gene_anno=gene_annotation,
                        plot_anno=plot_annotation,
                        title=title,
                        gene_type="all",
                        mean=F){

  
  if(mean==T){
    sampleDist$GENEID <- row.names(sampleDist)
    gene_anno <- gene_anno[match(rownames(sampleDist), gene_anno$GENEID),]
    sampleDist <- merge(sampleDist,
                        gene_anno,
                        by = "GENEID")
    
    rownames(sampleDist) <- sampleDist$GENEID
    sampleDist<-mean_function(input=sampleDist,
                              anno=sample_table,
                              condition="condition")
    
    sampleDist <- sampleDist[,colnames(sampleDist) %in% sample_table[["condition"]]]
  }else{
    sampleDist<-sampleDist
  }
  if(gene_type=="all"){
    
    if(mean==T){
      sampleDist <- sampleDist[,colnames(sampleDist) %in% sample_table[["condition"]]]
      sampleDist <- as.matrix(dist(t(sampleDist)))
      rownames(sampleDist)<- unique(sample_table$condition)
      colnames(sampleDist)<- unique(sample_table$condition)
      
    }else{
      sampleDist <- sampleDist[,colnames(sampleDist) %in% sample_table[["ID"]]]
      sampleDist <- as.matrix(dist(t(sampleDist)))
      rownames(sampleDist)<- sample_table$ID
      colnames(sampleDist)<- sample_table$ID
    }
    
    pheatmap(sampleDist,
             clustering_distance_rows =as.dist(sampleDist),
             clustering_distance_cols =as.dist(sampleDist),
             main="Sample distances based on variance-stabilized counts per sample",
             annotation_row = plot_anno, 
             annotation_col = plot_anno,
             annotation_colors = ann_colors,
             fontsize = 8)
  }else{
    #filtering if gene_type is not "all"
    sampleDist$GENEID <- row.names(sampleDist)
    gene_anno <- gene_anno[match(rownames(sampleDist), gene_anno$GENEID),]
    sampleDist <- merge(sampleDist,
                        gene_anno,
                        by = "GENEID")
    
    rownames(sampleDist) <- sampleDist$GENEID
    sampleDist<-sampleDist[sampleDist[["GENETYPE"]]==gene_type,]
    
    if(mean==T){
      sampleDist <- sampleDist[,colnames(sampleDist) %in% sample_table[["condition"]]]
      sampleDist <- as.matrix(dist(t(sampleDist)))
      rownames(sampleDist)<- unique(sample_table$condition)
      colnames(sampleDist)<- unique(sample_table$condition)
    }else{
      sampleDist <- sampleDist[,colnames(sampleDist) %in% sample_table[["ID"]]]
      sampleDist <- as.matrix(dist(t(sampleDist)))
      rownames(sampleDist)<- sample_table$ID
      colnames(sampleDist)<- sample_table$ID
    }
    
    pheatmap(sampleDist,
             clustering_distance_rows =as.dist(sampleDist),
             clustering_distance_cols =as.dist(sampleDist),
             main=title,
             annotation_row = plot_anno, 
             annotation_col = plot_anno,
             annotation_colors = ann_colors,
             fontsize = 8)
  }
}
