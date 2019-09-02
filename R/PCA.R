#' Function to plot a PCA
#'
#'
#' @param pca_input input dataset, default: dds_vst
#' @param ntop number of highest variable genes used for PCA
#' @param xPC Number of principal component on the x axis
#' @param yPC Number of principal component on the y axis
#' @param color Column of sample_table that defines the colour
#' @param anno_colour annotation of the colour values
#' @param shape Column of sample_table that defines the shape
#' @param point_size Size of points in plot
#' @param title Plot title
#' @param label Column of sample_table that defines the label
#' @param label_subset A character vector with a subset of labels to plot
#' @export


plotPCA <- function(pca_input = dds_vst,
                    ntop=500, 
                    xPC=1, 
                    yPC=2,
                    color,
                    anno_colour,
                    shape="NULL",
                    point_size=3,
                    title="PCA",
                    label=NULL,
                    label_subset=NULL,
                    gene_anno=gene_annotation,
                    gene_type="all"){
  
  if(gene_type=="all"){
    if(!is.data.frame(pca_input)){
      vst_matrix <-as.matrix(assay(pca_input))
    }else{
      vst_matrix <- pca_input
    }
    
    if(ntop=="all"){
      pca <- prcomp(t(vst_matrix)) 
    }else{
      # select the ntop genes by variance
      select <- order(rowVars(vst_matrix), decreasing=TRUE)[c(1:ntop)]
      pca <- prcomp(t(vst_matrix[select,]))
    }
    
    #calculate explained variance per PC
    explVar <- pca$sdev^2/sum(pca$sdev^2)
    # transform variance to percent
    percentVar <- round(100 * explVar[c(xPC,yPC)], digits=1)
    
    # Define data for plotting  
    pcaData <- data.frame(xPC=pca$x[,xPC], 
                          yPC=pca$x[,yPC], 
                          color = sample_table[[color]],
                          name= as.character(sample_table$ID),
                          stringsAsFactors = F)
    
    #plot PCA
    if(is.factor(pcaData$color) || is.character(pcaData$color)|| is.integer(pcaData$color)){
      if(shape == "NULL"){
        pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color)) +
          geom_point(size =point_size)
      }else{
        pcaData$shape = sample_table[[shape]]
        pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color, shape=shape)) +
          geom_point(size =point_size) +
          scale_shape_discrete(name=shape)
      }
      
      if(anno_colour[1] == "NULL"){
        pca_plot <- pca_plot + scale_color_discrete(name=color)
      }else{
        pca_plot <- pca_plot + scale_color_manual(values=anno_colour, name=color)
      }
      
    }else if(is.numeric(pcaData$color)){
      if(shape == "NULL"){
        pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color)) +
          geom_point(size =point_size) +
          scale_color_gradientn(colours = bluered(100),name=color)
      }else{
        pcaData$shape = sample_table[[shape]]
        pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color, shape=shape)) +
          geom_point(size =point_size) +
          scale_color_gradientn(colours = bluered(100),name=color)+
          scale_shape_discrete(name=shape)
      }
    }
    
    # adds a label to the plot. To label only specific points, put them in the arument label_subset
    if (!is.null(label) == TRUE){
      pcaData$label <- sample_table[[label]]
      if(!is.null(label_subset) == TRUE){
        pcaData_labeled <- pcaData[pcaData$label %in% label_subset,]
      } else {
        pcaData_labeled <- pcaData
      }
      pca_plot <- pca_plot +
        geom_text_repel(data = pcaData_labeled, aes(label = label), nudge_x = 2, nudge_y = 2, colour = "black")
    }
    
    pca_plot <- pca_plot+
      xlab(paste0("PC ",xPC, ": ", percentVar[1], "% variance")) +
      ylab(paste0("PC ",yPC,": ", percentVar[2], "% variance")) +
      coord_fixed()+
      theme_classic()+        
      theme(aspect.ratio = 1)+
      ggtitle(title)
    
    pca_plot
    #else for gene_type
  }else{
    if(!is.data.frame(pca_input)){
      vst_matrix <-as.matrix(assay(pca_input))
    }else{
      vst_matrix <- pca_input
    }
    
    #Deniz:filtering for gene_type 
    vst_matrix<-as.data.frame(vst_matrix)
    vst_matrix$GENEID <- row.names(vst_matrix)
    gene_anno <- gene_anno[match(rownames(vst_matrix), gene_anno$GENEID),]
    
    vst_matrix <- merge(vst_matrix,
                        gene_anno,
                        by = "GENEID")
    rownames(vst_matrix) <- vst_matrix$GENEID
    vst_matrix<-vst_matrix[vst_matrix[["GENETYPE"]]==gene_type,]
    vst_matrix <- vst_matrix[,colnames(vst_matrix) %in% sample_table[["ID"]]]
    vst_matrix<-as.matrix(vst_matrix)
    
    if(ntop=="all"){
      pca <- prcomp(t(vst_matrix)) 
    }else{
      # select the ntop genes by variance
      select <- order(rowVars(vst_matrix), decreasing=TRUE)[c(1:ntop)]
      pca <- prcomp(t(vst_matrix[select,]))
    }
    
    
    #calculate explained variance per PC
    explVar <- pca$sdev^2/sum(pca$sdev^2)
    # transform variance to percent
    percentVar <- round(100 * explVar[c(xPC,yPC)], digits=1)
    
    # Define data for plotting  
    pcaData <- data.frame(xPC=pca$x[,xPC], 
                          yPC=pca$x[,yPC], 
                          color = sample_table[[color]],
                          name= as.character(sample_table$ID),
                          stringsAsFactors = F)
    
    #plot PCA
    if(is.factor(pcaData$color) || is.character(pcaData$color)|| is.integer(pcaData$color)){
      if(shape == "NULL"){
        pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color)) +
          geom_point(size =point_size)
      }else{
        pcaData$shape = sample_table[[shape]]
        pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color, shape=shape)) +
          geom_point(size =point_size) +
          scale_shape_discrete(name=shape)
      }
      
      if(anno_colour[1] == "NULL"){
        pca_plot <- pca_plot + scale_color_discrete(name=color)
      }else{
        pca_plot <- pca_plot + scale_color_manual(values=anno_colour, name=color)
      }
      
    }else if(is.numeric(pcaData$color)){
      if(shape == "NULL"){
        pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color)) +
          geom_point(size =point_size) +
          scale_color_gradientn(colours = bluered(100),name=color)
      }else{
        pcaData$shape = sample_table[[shape]]
        pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color, shape=shape)) +
          geom_point(size =point_size) +
          scale_color_gradientn(colours = bluered(100),name=color)+
          scale_shape_discrete(name=shape)
      }
    }
    
    # adds a label to the plot. To label only specific points, put them in the argument label_subset
    if (!is.null(label) == TRUE){
      pcaData$label <- sample_table[[label]]
      if(!is.null(label_subset) == TRUE){
        pcaData_labeled <- pcaData[pcaData$label %in% label_subset,]
      } else {
        pcaData_labeled <- pcaData
      }
      pca_plot <- pca_plot +
        geom_text_repel(data = pcaData_labeled, aes(label = label), nudge_x = 2, nudge_y = 2, colour = "black")
    }
    
    pca_plot <- pca_plot+
      xlab(paste0("PC ",xPC, ": ", percentVar[1], "% variance")) +
      ylab(paste0("PC ",yPC,": ", percentVar[2], "% variance")) +
      coord_fixed()+
      theme_classic()+        
      theme(aspect.ratio = 1)+
      ggtitle(title)
    
    pca_plot
  }
}

#' Function to plot heatmaps of PC loadings
#'
#' @param PC principal component of which the loadings should be plotted
#' @param ntop number of genes to plot, sorted by variance
#' @export
#'
plotLoadings <- function(input=dds_vst_df, #or removedbatch_dds_vst
                         PC, 
                         ntop){
  if(ntop=="all"){
    pca <- prcomp(t(input)) 
  }else{
    select <- order(rowVars(input), decreasing=TRUE)[c(1:ntop)]
    pca <- prcomp(t(input[select,]))
  }
  
  Loadings <- pca$rotation[,PC]
  Loadings <- Loadings[order(Loadings, decreasing = T)]
  Loadings <- names(Loadings[c(1:20,(length(Loadings)-19):length(Loadings))])
  
  heatmap <- norm_anno[norm_anno$GENEID %in% Loadings,]
  rownames(heatmap) <- paste(heatmap$GENEID,": ",heatmap$SYMBOL,sep="")  
  heatmap <- heatmap[,colnames(heatmap) %in% sample_table$ID]
  heatmap_scale <- as.matrix(t(scale(t(heatmap))))
  
  # Heatmap
  pheatmap(heatmap_scale,
           main=paste("Hierarchical Clustering of top20 ",PC, " loadings in both directions",sep=""),
           show_rownames=TRUE,
           show_colnames = TRUE,
           annotation_col = plot_annotation,
           annotation_colors = ann_colors,
           breaks = scaleColors(heatmap_scale, 2)[["breaks"]], 
           color = scaleColors(heatmap_scale, 2)[["color"]],
           cluster_cols = T,
           fontsize=6)
}

# 3D PCA
plot3D_pca<-function(pca3d_input=dds_vst,
                     sample_table=sample_table,
                     gene_anno=gene_annotation,
                     gene_type="all",
                     title="3D Scatter plot_PCA",
                     xPC=1,
                     yPC=2,
                     zPC=3,
                     ntop=500,
                     anno_colour=col_condition,
                     point_size=3){
  
  samplePCA_3d<-as.matrix(assay(pca3d_input))
  
  if(gene_type=="all"){
    if(ntop=="all"){
      pca <- prcomp(t(samplePCA_3d)) 
    }else{
      # select the ntop genes by variance
      select <- order(rowVars(samplePCA_3d), decreasing=TRUE)[c(1:ntop)]
      pca <- prcomp(t(samplePCA_3d[select,]))
    }
    
    #calculate explained variance per PC
    explVar <- pca$sdev^2/sum(pca$sdev^2)
    # transform variance to percent
    percentVar <- round(100 * explVar[c(xPC,yPC,zPC)], digits=1)
    
    # Define data for plotting  
    pcaData_3D <- data.frame(xPC=pca$x[,xPC], 
                             yPC=pca$x[,yPC],
                             zPC=pca$x[,zPC],
                             condition = sample_table$condition,
                             ID= as.character(sample_table$ID),
                             stringsAsFactors = F)
    
    
    pcaData_3D$condition <- as.factor(pcaData_3D$condition)
    
    
    p <- plot_ly(pcaData_3D, x = ~xPC, y = ~yPC, z = ~zPC, color = ~condition, colors = anno_colour,
                 text = ~paste('ID:', ID), marker = list(size = point_size)) %>% #, '<br>Genotype_Stim:', Genotype_Stim, '<br>Preparation:', Preparation
      layout(title=title,
             scene = list(xaxis = list(title = paste0("PC ",xPC,": ", percentVar[1], "% variance")),
                          yaxis = list(title = paste0("PC ",yPC,": ", percentVar[2], "% variance")),
                          zaxis = list(title = paste0("PC ",zPC,": ", percentVar[3], "% variance"))))
    
  }else{
    #filtering
    samplePCA_3d<-as.data.frame(samplePCA_3d)
    samplePCA_3d$GENEID <- row.names(samplePCA_3d)
    gene_anno <- gene_anno[match(rownames(samplePCA_3d), gene_anno$GENEID),]
    
    samplePCA_3d <- merge(samplePCA_3d,
                          gene_anno,
                          by = "GENEID")
    rownames(samplePCA_3d) <- samplePCA_3d$GENEID
    samplePCA_3d<-samplePCA_3d[samplePCA_3d[["GENETYPE"]]==gene_type,]
    samplePCA_3d <- samplePCA_3d[,colnames(samplePCA_3d) %in% sample_table[["ID"]]]
    samplePCA_3d<-as.matrix(samplePCA_3d)
    ####
    if(ntop=="all"){
      pca <- prcomp(t(samplePCA_3d)) 
    }else{
      # select the ntop genes by variance
      select <- order(rowVars(samplePCA_3d), decreasing=TRUE)[c(1:ntop)]
      pca <- prcomp(t(samplePCA_3d[select,]))
    }
    
    #calculate explained variance per PC
    explVar <- pca$sdev^2/sum(pca$sdev^2)
    # transform variance to percent
    percentVar <- round(100 * explVar[c(xPC,yPC,zPC)], digits=1)
    
    # Define data for plotting  
    pcaData_3D <- data.frame(xPC=pca$x[,xPC], 
                             yPC=pca$x[,yPC],
                             zPC=pca$x[,zPC],
                             condition = sample_table$condition,
                             ID= as.character(sample_table$ID),
                             stringsAsFactors = F)
    
    
    pcaData_3D$condition <- as.factor(pcaData_3D$condition)
    
    p <- plot_ly(pcaData_3D, x = ~xPC, y = ~yPC, z = ~zPC, color = ~condition, colors = anno_colour,
                 text = ~paste('ID:', ID),marker = list(size = point_size)) %>% #, '<br>Genotype_Stim:', Genotype_Stim, '<br>Preparation:', Preparation
      layout(title=title,
             scene = list(xaxis = list(title = paste0("PC ",xPC,": ", percentVar[1], "% variance")),
                          yaxis = list(title = paste0("PC ",yPC,": ", percentVar[2], "% variance")),
                          zaxis = list(title = paste0("PC ",zPC,": ", percentVar[3], "% variance"))))
    
  }
  p
}



