#' Function to plot a PCA
#'
#'
#' @param pca_input input dataset, default: dds_vst
#' @param pca_sample_table sample table with meta annotation that should be used
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
                    pca_sample_table = sample_table,
                    ntop=500,
                    xPC=1,
                    yPC=2,
                    color,
                    anno_colour,
                    shape="NULL",
                    point_size=3,
                    title="PCA",
                    label = NULL,
                    label_subset = NULL){

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
                        color = pca_sample_table[[color]],
                        name= as.character(pca_sample_table$ID),
                        stringsAsFactors = F)

  #plot PCA
  if(is.factor(pcaData$color) || is.character(pcaData$color)|| is.integer(pcaData$color)){
    if(shape == "NULL"){
      pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color)) +
        geom_point(size =point_size)
    }else{
      pcaData$shape = pca_sample_table[[shape]]
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
      pcaData$shape = pca_sample_table[[shape]]
      pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color, shape=shape)) +
        geom_point(size =point_size) +
        scale_color_gradientn(colours = bluered(100),name=color)+
        scale_shape_discrete(name=shape)
    }
  }

  # adds a label to the plot. To label only specific points, put them in the arument label_subset
  if (!is.null(label) == TRUE){
    pcaData$label <- pca_sample_table[[label]]
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

#' Function to plot heatmaps of PC loadings
#'
#' @param PC principal component of which the loadings should be plotted
#' @param pca_input data to calculate the PCA on (usually the variance stabilized counts)
#' @param heatmap_input data to show in heatmap (usually norm_anno)
#' @param ntop number of genes to plot, sorted by variance
#' @export
#'
plotLoadings <- function(pca_input = dds_vst, heatmap_input = norm_anno, sample_annotation = sample_table, PC, ntop){
  if(ntop=="all"){
    pca <- prcomp(t(assay(pca_input)))
  }else{
    select <- order(rowVars(assay(pca_input)), decreasing=TRUE)[c(1:ntop)]
    pca <- prcomp(t(assay(pca_input)[select,]))
  }

  Loadings <- pca$rotation[,PC]
  Loadings <- Loadings[order(Loadings, decreasing = T)]
  Loadings <- names(Loadings[c(1:20,(length(Loadings)-19):length(Loadings))])

  heatmap <- heatmap_input[heatmap_input$GENEID %in% Loadings,]
  rownames(heatmap) <- paste(heatmap$GENEID,": ",heatmap$SYMBOL,sep="")
  heatmap <- heatmap[,colnames(heatmap) %in% sample_annotation$ID]
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


