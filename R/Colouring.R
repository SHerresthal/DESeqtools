#' Creating colour annotations
#'
#' This functions creates colour annotations for categorical and continuous variables in a data frame to be plotted in a heatmap.
#'
#' @param data The input dataframe
#' @param factorcolours A vector of names of colour palettes from the RColourbrewer package. These palettes will be used to colour categorical variables
#' @param continuouscolors A vector of colour names that will be chosen for continuous data
#'
#' @return Returns a list with colours that can be used as an input for pheatmap
#' @export

setColours <- function(data = plot_annotation,
                       factorcolours = c("Set1","Greys"),
                       continouscolours = c("forestgreen")){
  ann_colors <- list()
  for(i in 1:ncol(plot_annotation)){
    tmp <- plot_annotation[,i]
    if(is.factor(tmp) == T){
      col <- brewer.pal(n = length(levels(tmp)), factorcolours[i])
      if(length(levels(tmp)) < 3) {
        col <- col[1:length(levels(tmp))]
      }
      names(col) <- levels(tmp)
    } else {
      tmp <- tmp[order(tmp)]
      col <- colorRampPalette(c("white", "forestgreen"))(length(tmp))
      names(col) <- tmp
    }
    ann_colors[[colnames(plot_annotation)[i]]] <- col
  }
  return(ann_colors)
}

#' Scaling colours in a heatmap
#'
#
#' Returns a palette of colours from blue to white to red with white ranging
#' from the lowest to the highest value. The white colour is set to 0.
#'
#'
#' @param data The data to use
#' @param maxvalue The value at which the colour should be fully red or blue. By
#'   default, this value is set to the abbsolute largest value
#'
#' @return Returns a list of breaks and colour values that can be used inside the \link[pheatmap:pheatmap]{pheatmap} function.
#' @export


scaleColors <- function(data = input_scale, # data to use
                        maxvalue = NULL # value at which the color is fully red / blue
){
  if(is.null(maxvalue)){
    maxvalue <- floor(min(abs(min(data)), max(data)))
  }
  if(max(data) > abs(min(data))){
    if(ceiling(max(data)) == maxvalue){
      myBreaks <- c(floor(-max(data)), seq(-maxvalue+0.2, maxvalue-0.2, 0.2),  ceiling(max(data)))
    } else{
      myBreaks <- c(floor(-max(data)), seq(-maxvalue, maxvalue, 0.2),  ceiling(max(data)))
    }
    paletteLength <- length(myBreaks)
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  } else {
    if(-floor(min(data)) == maxvalue){
      myBreaks <- c(floor(min(data)), seq(-maxvalue+0.2, maxvalue-0.2, 0.2),  ceiling(min(data)))
    } else{
      myBreaks <- c(floor(min(data)), seq(-maxvalue, maxvalue, 0.2),  ceiling(abs(min(data))))
    }
    paletteLength <- length(myBreaks)
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  }
  return(list(breaks = myBreaks, color = myColor))
}
