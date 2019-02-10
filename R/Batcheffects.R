#' Removes batch effects
#'
#' Function to remove potential batch effects using the removeBatchEffect function from limma
#' @export

limmaBatchEffectRemoval <- function(input=dds_vst,
                                    batchfactor, # name of batch effect column in sample_table
                                    batchfactor_2=NULL,
                                    modelfactor){ # name of model effect column in sample_table

  # rlog-transformed input
  x <- as.matrix(assay(input))

  # design matrix
  model <- model.matrix(~sample_table[,c(modelfactor)])

  # run batch remocal function
  if(is.numeric(sample_table[,colnames(sample_table) == batchfactor[1]])==T){
    as.data.frame(removeBatchEffect(x,
                                    covariates = sample_table[,colnames(sample_table) %in% batchfactor],
                                    design = model))
  }else{
    if(is.null(batchfactor_2)){
      as.data.frame(removeBatchEffect(x=x,
                                      batch = sample_table[,colnames(sample_table) == batchfactor],
                                      design = model))
    }else{
      as.data.frame(removeBatchEffect(x=x,
                                      batch = sample_table[,colnames(sample_table) == batchfactor],
                                      batch2 = sample_table[,colnames(sample_table) == batchfactor_2],
                                      design = model))
    }
  }
}
