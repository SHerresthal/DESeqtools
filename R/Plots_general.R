#' Function for plotting multiple plots in grid
#' @export

multiplot<-function(plots=plots,
                    cols=1){

  layout <- matrix(seq(1, cols * length(plots)/cols),
                   ncol = cols,
                   nrow = length(plots)/cols)


  if (length(plots)==1) {
    print(plots[[1]])
  }else{
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:length(plots)){
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#'Venn Diagram
#'@export

plotVenn <- function(comparisons,
                     regulation=NULL){
  venn <- NULL
  for(i in 1:length(comparisons)){
    res <- DEresults[names(DEresults) %in% comparisons]
    comp <- as.data.frame(res[[i]]@results)
    if(is.null(regulation)){
      DE <- ifelse(comp$regulation %in% c("up","down"), 1, 0)
      venn <- cbind(venn, DE)
      colnames(venn)[i]<- paste(names(res)[[i]], "up&down", sep=": ")
    } else {
      DE <- ifelse(comp$regulation == regulation, 1, 0)
      venn <- cbind(venn, DE)
      colnames(venn)[i]<- paste(names(res)[[i]], regulation, sep=": ")
    }

  }
  vennDiagram(venn,cex = 1, counts.col = "blue")
}

#' Ratio plot function
#' @export
plotRatios <- function(comp1, comp2){
  U <- NULL
  c <- c(comp1,comp2)
  U <- uDEG(comparisons = c, keyType = "Ensembl")
  Ratio <- NULL
  for(i in 1:length(c)){
    tmp <- DEresults[names(DEresults) %in% c]
    comp <- as.data.frame(tmp[[i]]@results)
    DE <- as.data.frame(comp[rownames(comp) %in% U,])
    Ratio <- as.data.frame(cbind(Ratio,DE$log2FoldChange))
  }
  colnames(Ratio)<- c
  rownames(Ratio) <- U
  ggplot(Ratio, aes(x=Ratio[,1], y=Ratio[,2])) +
    geom_point(colour = "grey", size = 1.5) +
    theme_bw() +
    xlab(comp1)+
    ylab(comp2) +
    geom_abline(slope = c(-1,1),intercept = 0, colour="grey") +
    geom_hline(yintercept = c(0,log(2,2),-(log(2,2))))+
    geom_hline(yintercept = c(log(2,2),-(log(2,2))), colour="firebrick1")+
    geom_vline(xintercept = c(log(2,2),-(log(2,2))))+
    geom_vline(xintercept = c(log(2,2),-(log(2,2))), colour="firebrick1")+
    theme(text = element_text(size=10))+
    ggtitle(paste(comp1," vs ",comp2,": ",length(U)," DE genes",sep=""))
}


#' Volcano Plot function
#' @export
plotVolcano <-  function(comparison,
                         labelnum=20){

  # specify labeling
  upDE <-  as.data.frame(DEresults[[comparison]]@results[DEresults[[comparison]]@results$regulation =="up",])
  FClabel_up <- upDE[order(abs(upDE$log2FoldChange), decreasing = TRUE),]
  if(nrow(FClabel_up)>labelnum){
    FClabel_up <- as.character(FClabel_up[c(1:labelnum),"GENEID"])
  } else {
    FClabel_up <- as.character(FClabel_up$GENEID)}
  plabel_up <- upDE[order(upDE$padj, decreasing = FALSE),]
  if(nrow(plabel_up)>labelnum){
    plabel_up <- as.character(plabel_up[c(1:labelnum),"GENEID"])
  } else {
    plabel_up <- as.character(plabel_up$GENEID)}

  downDE <-  as.data.frame(DEresults[[comparison]]@results[DEresults[[comparison]]@results$regulation =="down",])
  FClabel_down <- downDE[order(abs(downDE$log2FoldChange), decreasing = TRUE),]
  if(nrow(FClabel_down)>labelnum){
    FClabel_down <- as.character(FClabel_down[c(1:labelnum),"GENEID"])
  } else {
    FClabel_down <- as.character(FClabel_down$GENEID)}
  plabel_down <- downDE[order(downDE$padj, decreasing = FALSE),]
  if(nrow(plabel_down)>labelnum){
    plabel_down <- as.character(plabel_down[c(1:labelnum),"GENEID"])
  } else {
    plabel_down <- as.character(plabel_down$GENEID)}


  label<- unique(c(FClabel_up, plabel_up, FClabel_down, plabel_down))

  data <- DEresults[[comparison]]@results
  data$label<- ifelse(data$GENEID %in% label == "TRUE",as.character(data$SYMBOL), "")
  data <- data[,colnames(data) %in% c("label", "log2FoldChange", "padj", "regulation")]

  # Volcano Plot
  ggplot(data=na.omit(data), aes(x=log2FoldChange, y=-log10(padj), colour=regulation)) +
    geom_point(alpha=0.4, size=1.75) +
    scale_color_manual(values=c("cornflowerblue","grey", "firebrick"))+
    scale_x_continuous() +
    scale_y_continuous() +
    xlab("log2(FoldChange)") +
    ylab("-log10(padj)") +
    geom_vline(xintercept = 0, colour="black")+
    geom_vline(xintercept = c(-log(2,2),log(2,2)), colour="red")+
    geom_hline(yintercept=-log(0.05,10),colour="red")+
    geom_text_repel(data=na.omit(data[!data$label =="",]),aes(label=label), size=3)+
    guides(colour=FALSE) +
    ggtitle(paste("Volcano Plot of: ",comparison,sep="")) +
    theme_bw()
}

#' Custom function to plot baseMean versus fold change
#' @export
plotMA <- function(comparison,
                   ylim=c(-2,2),
                   padjThreshold=0.05,
                   xlab = "mean of normalized counts",
                   ylab = expression(log[2]~fold~change),
                   log = "x",
                   cex=0.45){
  x <- as.data.frame(DEresults[[comparison]]@results)
  if (!(is.data.frame(x) && all(c("baseMean", "log2FoldChange") %in% colnames(x)))){
    stop("'x' must be a data frame with columns named 'baseMean', 'log2FoldChange'.")
  }
  col = ifelse(x$padj>=padjThreshold, "gray32", "red3")
  py = x$log2FoldChange
  if(missing(ylim)){
    ylim = c(-1,1) * quantile(abs(py[is.finite(py)]), probs=0.99) * 1.1
  }
  plot(x=x$baseMean,
       y=pmax(ylim[1], pmin(ylim[2], py)),
       log=log,
       pch=ifelse(py<ylim[1], 6, ifelse(py>ylim[2], 2, 16)),
       cex=cex,
       col=col,
       xlab=xlab,
       ylab=ylab,
       ylim=ylim,
       main=comparison)
  abline(h=0, lwd=4, col="#ff000080")
  abline(h=c(-1,1), lwd=2, col="dodgerblue")
}

#' Plot p value distribution
#'
#' @export

plotPvalues <- function(comparison){
  res <- as.data.frame(DEresults[[comparison]]@results)
  ggplot(na.omit(res), aes(x=pvalue)) +
    geom_histogram(aes(y=..count..),
                   binwidth = 0.01) +
    theme_bw()+
    ggtitle(paste("p value histogram of: ",comparison,sep=""))
}


#' Ranked Fold Change plot
#'
#' @export

plotFCrank <- function(comp1,
                       comp2){
  rank <- na.omit(DEresults[names(DEresults) == comp1][[1]]@results)
  rank <- rank[rank$padj < 0.05 , c("GENEID","comparison","log2FoldChange")]
  rank <- rank[order(rank$log2FoldChange,decreasing = TRUE),]
  rank$rank <- c(1:nrow(rank))
  rank2 <- DEresults[names(DEresults) == comp2][[1]]@results
  rank2 <- rank2[rownames(rank),c("GENEID","comparison","log2FoldChange")]
  rank2$rank <- rank$rank
  rank <- rbind(rank, rank2)
  ggplot(rank,aes(x=rank,y=log2FoldChange,color=comparison)) +
    geom_point(alpha=0.5) +
    geom_line(aes(group=GENEID),color="grey",alpha=0.2)+
    theme_bw() +
    ylab("log2(FoldChange)")+
    xlab(paste("FC rank of " ,comp1, sep=""))+
    geom_hline(yintercept = 0)+
    geom_hline(yintercept = c(log(2,2),-(log(2,2))), colour="firebrick1")+
    theme(text = element_text(size=10))+
    ggtitle("Comparison of fold changes (comp1 padj<0.05")+
    theme(legend.position="bottom")
}

