#' Mean of replicates  
#'
#' This function provides the calculation of the average of replicates per condition.
#' 
#'
#' @param input A dataframe containing the read counts of the genes for each sample. By default, "norm_anno" normalized read counts is used.
#' @param anno A sample table containing info for the samples
#' @param condition A character from sample table df defining the name of column in the table to use. "condition" is set by default.
#' @return Returns a dataframe
#' @export


mean_function<- function(input=norm_anno, 
                         anno=sample_table,
                         condition="condition"){
  
  conditions <- unique(anno[,colnames(anno) == condition])
  
  df <- data.frame(matrix(nrow = nrow(input),
                          ncol= length(conditions)))
  colnames(df) <- as.character(conditions)
  rownames(df) <- norm_anno$GENEID
  
  for(i in conditions){
    i <- as.character(i)
    #print(i)
    tmp <- input[,colnames(input) %in%
                   anno$ID[anno[,colnames(anno)==condition]==i]]
    if(class(tmp)=="numeric"){
      tmp<-as.data.frame(tmp)
      colnames(tmp)<-i
      df[,i]<-tmp
    }else{
      df[,i]<- rowMeans(tmp)
    }
  }
  
  df$GENEID <- row.names(df)
  rownames(df) <- df$GENEID
  df <- merge(df, gene_annotation, by = "GENEID")
  rownames(df) <- df$GENEID
  
  return(df)
}


#' Generate mean sample table 
#'
#' This function provides the calculation of mean of the total reads, percent aligned etc. 
#' in sample table to make mean_sample_table of the average of replicates per condition.
#' 
#'
#' @param input A dataframe containing the read counts of the genes for each sample. By default, "norm_anno" normalized read counts is used.
#' @param anno A sample table containing info for the samples
#' @param condition A character from sample table df defining the name of column in the table to use. "condition" is set by default.
#' @return Returns a dataframe
#' @export

mean_sample_definition<-function(input=mean_sample_table,
                                 anno=sample_table,
                                 condition="condition"){
  
  
  conditions <- unique(anno[,colnames(anno) == condition])
  anno$ID<-as.character(anno$ID)
  
  for(i in conditions){
    i <- as.character(i)
    print(i)
    num<-select_if(anno,is.numeric)
    num$condition<-anno[,colnames(anno) == condition]
    tmp <- anno[rownames(anno) %in%
                  rownames(num)[num[,colnames(num)==condition]==i],]
    tmp<-select_if(tmp,is.numeric)
    mean_sample_table[i,colnames(mean_sample_table) %in% colnames(tmp)]<- t(as.matrix(colMeans(tmp)))
  }
  mean_sample_table$ID<-mean_sample_table[[condition]]
  return(mean_sample_table)
}

