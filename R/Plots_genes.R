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
        if (is.null(shape_opt)){
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


#' Plot BoxPlot of highest expressed genes
#'
#' @export

highestGenes <- function(numGenes=10){
  tmp <- norm_anno[,colnames(norm_anno) %in% sample_table$ID]
  tmp <- tmp[order(rowMeans(tmp), decreasing = T),]
  tmp <- tmp[1:numGenes,]
  tmp <- melt(t(tmp))
  colnames(tmp)<- c("sample","gene","value")

  idx <- match(tmp$gene,norm_anno$GENEID)
  tmp$symbol <- as.factor(norm_anno$SYMBOL[idx])
  tmp$symbol <- factor(tmp$symbol, levels = rev(unique(tmp$symbol)))

  ggplot(tmp, aes(x = tmp$symbol, y = value)) +
    geom_boxplot()+
    xlab("Gene")+
    ylab("Normalized Expression")+
    ggtitle(paste("Expression of", numGenes, "highest expressed genes")) +
    theme_bw() +
    coord_flip() +
    theme(axis.text.x = element_text(size=8, angle = 90, hjust = 1),
          plot.title = element_text(size = 8, face = "bold"))
}

