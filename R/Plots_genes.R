#' Function to plot the normalized or batch-corrected counts for a single gene in a box plot
#'
#' @export

plotSingleGene <-function(data=norm_anno,
                          symbol,
                          condition="Genotype_Age",
                          anno_colour=col_genotype_age,
                          shape = NULL) {

  input<-as.data.frame(data)
  rownames(input)<- input$GENEID

  if(sum(input$SYMBOL == symbol) == 0){
    stop("Gene not present")
  }else{
    plots<-list()
    for (i in 1:sum(input$SYMBOL == symbol)) {
      geneCounts <- as.data.frame(t(input[input$SYMBOL == symbol, colnames(input) %in% sample_table$ID]))
      geneCounts$condition <- sample_table[[condition]]
      GENEID<-colnames(geneCounts)[i]
      colnames(geneCounts)[i]<-"y"

      if(!is.null(anno_colour)){
        if (is.null(shape)){
          plot<-ggplot(geneCounts, aes(x = condition, y = y, colour=condition)) +
            scale_color_manual(values=anno_colour)+
            geom_beeswarm(cex = 3, na.rm=T)
        }else{
          geneCounts$shape <- sample_table[[shape]]
          legend_shape<-paste0(shape)
          plot<-ggplot(geneCounts, aes(x = condition, y = y, colour=condition)) +
            scale_color_manual(values=anno_colour)+
            geom_beeswarm(cex = 3, na.rm=T, aes(shape=shape)) +
            scale_shape(name=legend_shape)
        }
      }else{
        if (is.null(shape)){
          plot<-ggplot(geneCounts, aes(x = condition, y = y, colour=condition)) +
            scale_color_brewer(palette = "Spectral")+
            geom_beeswarm(cex = 3, na.rm=T, aes(size=3))
        }else{
          geneCounts$shape <- sample_table[[shape]]
          legend_shape<-paste0(shape)
          plot<-ggplot(geneCounts, aes(x = condition, y = y, colour=condition)) +
            scale_color_brewer(palette = "Spectral")+
            geom_beeswarm(cex = 3, na.rm=T, aes(shape=shape)) +
            scale_shape(name=legend_shape)
        }
      }
      plots[[i]]<-plot+
        geom_boxplot(width=.5,alpha=0) +
        stat_boxplot(geom ='errorbar',width=.25) +
        ylab("Normalized counts") +
        scale_y_continuous(expand=c(0.05,0.25)) +
        expand_limits(y=0) +
        labs(title=paste(symbol, GENEID, sep=": "),colour=condition)+
        theme_classic()+
        theme(plot.title = element_text(hjust=0.5))
    }
    if(sum(input$SYMBOL== symbol)>1){
      print("Selected gene symbol assigned to more than one gene (Ensembl ID)")
      multiplot(plots)
    }else{
      print("Selected gene symbol assigned to one gene (Ensembl ID)")
      multiplot(plots)
    }
  }
}


#' Plot BoxPlot of highest expressed genes, ranked by mean expression over all samples
#'
#' @param numGenes Number of genes to show
#' @param data The data to plot, should be in the format of norm_anno
#' @param plot_title The title of the plot (default is "Expression of *numgenes* highest expressed genes
#' @return Prints a ggplot object
#'
#' @export

highestGenes <- function(numGenes = 10, data = norm_anno, plot_title = NULL)
{
  tmp <- data[, colnames(data) %in% sample_table$ID]
  tmp <- tmp[order(rowMeans(tmp), decreasing = T), ] # order according to maximal mean expression value
  tmp <- tmp[1:numGenes, ]
  tmp <- melt(t(tmp))
  colnames(tmp) <- c("sample", "gene", "value")
  idx <- match(tmp$gene, data$GENEID)
  tmp$symbol <- as.factor(data$SYMBOL[idx])
  tmp$symbol <- factor(tmp$symbol, levels = rev(unique(tmp$symbol)))
  tmp$description <- data$DESCRIPTION[idx]
  tmp$genetype <- data$GENETYPE[idx]

  if(is.null(plot_title)){
    title <- paste("Expression of", numGenes, "highest expressed genes")
  } else {
    title <- plot_title
  }

  ggplot(tmp, aes(x = tmp$symbol, y = value)) +
    geom_quasirandom(size = 0.8) +
    geom_boxplot(aes(fill = genetype), alpha = 0.6, outlier.shape = NA) +
    scale_fill_brewer(palette = "Paired", name = "Gene Type") +
    xlab("") +
    ylab("Normalized Expression") +
    ggtitle(paste(title)) +
    theme_bw() +
    coord_flip() +
    theme(axis.text.x = element_text(size = 18, angle = 90, hjust = 1),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          plot.title = element_text(size = 18, face = "bold")) +
    guides(fill = guide_legend(override.aes = list(alpha = 1)))
}

