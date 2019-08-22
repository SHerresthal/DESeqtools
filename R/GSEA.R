#' GSEA function
#' @export

GSEA <-  function(comparison,
                  organism,
                  GeneSets =c("GO","KEGG","Hallmark","cannonicalPathways","Motifs","ImmunoSignatures"),
                  GOntology = "BP",
                  pCorrection = "bonferroni", # choose the p-value adjustment method
                  pvalueCutoff = 0.05, # set the unadj. or adj. p-value cutoff (depending on correction method)
                  qvalueCutoff = 0.05, # set the q-value cutoff (FDR corrected)
                  showMax = 20,
                  font.size = 8){
  
  results <- list()
  
  if(organism == "mouse") {
    OrgDb = org.Mm.eg.db
  } else if(organism == "human"){
    OrgDb = org.Hs.eg.db
  } else {print("Wrong Organism. Select mouse or human.")}
  
  res <- DEresults[[comparison]]
  DE_up <- as.data.frame(res@DE_genes$up_regulated_Genes)$SYMBOL
  entrez_up <- bitr(DE_up, fromType = "SYMBOL", toType="ENTREZID", OrgDb=OrgDb)$ENTREZID
  DE_down <- as.data.frame(res@DE_genes$down_regulated_Genes)$SYMBOL
  entrez_down <- bitr(DE_down, fromType = "SYMBOL", toType="ENTREZID", OrgDb=OrgDb)$ENTREZID  
  
  # GO enrichment
  if("GO" %in% GeneSets){
    print("Performing GO enrichment")
    if(length(entrez_up)<20){
      print("Too few upregulated genes for GO enrichment (<20)")
      results$GO_up <- "Too few upregulated genes for GO enrichment (<20)"
    }else{
      eGO_up <- enrichGO(gene = entrez_up,
                         universe = universe_Entrez,
                         OrgDb = OrgDb,
                         ont = GOntology,
                         pAdjustMethod = pCorrection,
                         pvalueCutoff  = pvalueCutoff,
                         qvalueCutoff  = qvalueCutoff,
                         readable      = T)
      
      results$GOup <- as.data.frame(eGO_up)
      if(nrow(results$GOup)<1){
        results$GOup_plot <- "No GO enrichment for upregulated genes"
      }else{
        results$GOup_plot <- clusterProfiler::dotplot(eGO_up, 
                                                      showCategory = showMax, 
                                                      font.size= font.size, 
                                                      title = paste("GO enrichment for genes upregulated in: ", comparison,sep="")
        )
      }
    }
    if(length(entrez_down)<20){
      print("Too few downregulated genes for GO enrichment (<20)")
      results$GO_down <- "Too few downregulated genes for GO enrichment (<20)"
    }else{
      eGO_down <- enrichGO(gene = entrez_down,
                           universe = universe_Entrez,
                           OrgDb = OrgDb,
                           ont = GOntology,
                           pAdjustMethod = pCorrection,
                           pvalueCutoff  = pvalueCutoff,
                           qvalueCutoff  = qvalueCutoff,
                           readable      = T)
      
      results$GOdown <- as.data.frame(eGO_down)
      if(nrow(results$GOdown)<1){
        results$GOdown_plot <- "No GO enrichment for downregulated genes"
      }else{
        results$GOdown_plot <- clusterProfiler::dotplot(eGO_down, 
                                                        showCategory = showMax, 
                                                        font.size= font.size, 
                                                        title = paste("GO enrichment for genes downregulated in: ", comparison,sep="")
        )
      }
    }
  }
  
  # KEGG enrichment
  if("KEGG" %in% GeneSets){
    print("Performing KEGG enrichment")
    
    if(organism == "mouse") {org = "mmu"} 
    if(organism == "human"){org = "hsa"}
    
    if(length(entrez_up)<20){
      print("Too few upregulated genes for KEGG enrichment (<20)")
      results$KEGG_up <- "Too few upregulated genes for KEGG enrichment (<20)"
    }else{
      eKEGG_up <- enrichKEGG(gene = entrez_up, 
                             organism = org,
                             universe = universe_Entrez, 
                             pAdjustMethod = pCorrection,
                             pvalueCutoff  = pvalueCutoff,
                             qvalueCutoff = qvalueCutoff)
      
      results$KEGGup <- as.data.frame(eKEGG_up)
      if(nrow(results$KEGGup)<1){
        results$KEGGup_plot <- "No KEGG enrichment for upregulated genes"
      }else{
        results$KEGGup_plot <- clusterProfiler::dotplot(eKEGG_up,  
                                                        showCategory = showMax, 
                                                        font.size= font.size, 
                                                        title = paste("KEGG enrichment for genes upregulated in: ",comparison,sep="")
        )
      }
    }
    if(length(entrez_down)<20){
      print("Too few downregulated genes for KEGG enrichment (<20)")
      results$KEGG_down <- "Too few downregulated genes for KEGG enrichment (<20)"
    } else{
      eKEGG_down <- enrichKEGG(gene = entrez_down, 
                               organism = org,
                               universe = universe_Entrez, 
                               pAdjustMethod = pCorrection,
                               pvalueCutoff  = pvalueCutoff,
                               qvalueCutoff = qvalueCutoff)
      
      results$KEGGdown <- as.data.frame(eKEGG_down)
      if(nrow(results$KEGGdown)<1){
        results$KEGGdown_plot <- "No KEGG enrichment for downregulated genes"
      }else{
        results$KEGGdown_plot <- clusterProfiler::dotplot(eKEGG_down,
                                                          showCategory = showMax,
                                                          font.size= font.size,
                                                          title = paste("KEGG enrichment for genes upregulated in: ",comparison,sep="")
        )
      }
    }
  }
  
  if("Hallmark" %in% GeneSets |"cannonicalPathways" %in% GeneSets| "ImmunoSignatures" %in% GeneSets | "Motifs" %in% GeneSets){
    genes_up_hsa <- getLDS(attributes = c("mgi_symbol"), 
                           filters = "mgi_symbol", 
                           values = DE_up, 
                           mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl"), 
                           attributesL = c("hgnc_symbol"), 
                           martL = useMart("ensembl", dataset = "hsapiens_gene_ensembl"), 
                           uniqueRows=T)[,2]
    entrez_up_hsa <- bitr(genes_up_hsa, fromType = "SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
    genes_down_hsa <- getLDS(attributes = c("mgi_symbol"), 
                             filters = "mgi_symbol", 
                             values = DE_down, 
                             mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl"), 
                             attributesL = c("hgnc_symbol"), 
                             martL = useMart("ensembl", dataset = "hsapiens_gene_ensembl"), 
                             uniqueRows=T)[,2]
    entrez_down_hsa <- bitr(genes_down_hsa, fromType = "SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  }
  
  # Hallmark enrichment
  if("Hallmark" %in% GeneSets){
    print("Performing Hallmark enrichment")
    if(length(entrez_up_hsa)<20){
      print("Too few upregulated genes for Hallmark enrichment (<20)")
      results$Hallmark_up <- "Too few upregulated genes for Hallmark enrichment (<20)"
    }else{
      Hallmark_up <- enricher(entrez_up_hsa,
                              TERM2GENE=hallmark_genes,
                              universe = universe_mouse2human_Entrez,  
                              pAdjustMethod = pCorrection,
                              pvalueCutoff  = pvalueCutoff,
                              qvalueCutoff = qvalueCutoff)
      
      results$HALLMARKup <- as.data.frame(Hallmark_up)
      if(nrow(results$HALLMARKup)<1){
        results$HALLMARKup_plot <- "No Hallmark enrichment for upregulated genes"
      }else{
        results$HALLMARKup_plot <- clusterProfiler::dotplot(Hallmark_up,
                                                            showCategory = showMax,
                                                            font.size= font.size,
                                                            title = paste("Hallmark enrichment for genes upregulated in: ",comparison,sep="")
        )
      }
    }
    if(length(entrez_down_hsa)<20){
      print("Too few downregulated genes for Hallmark enrichment (<20)")
      results$Hallmark_down <- "Too few downregulated genes for Hallmark enrichment (<20)"
    }else{
      Hallmark_down <- enricher(entrez_down_hsa,
                                TERM2GENE=hallmark_genes,
                                universe = universe_mouse2human_Entrez,  
                                pAdjustMethod = pCorrection,
                                pvalueCutoff  = pvalueCutoff,
                                qvalueCutoff = qvalueCutoff)
      
      results$HALLMARKdown <- as.data.frame(Hallmark_down)
      if(nrow(results$HALLMARKdown)<1){
        results$HALLMARKdown_plot <-"No Hallmark enrichment for downregulated genes"
      }else{
        results$HALLMARKdown_plot <- clusterProfiler::dotplot(Hallmark_down,
                                                              showCategory = showMax,
                                                              font.size= font.size,
                                                              title = paste("Hallmark enrichment for genes downregulated in: ",comparison,sep="")
        )
      }
    }
  }
  
  # Cannonical Pathway enrichment
  if("cannonicalPathways" %in% GeneSets){
    print("Performing Cannonical Pathway (C2) enrichment")
    if(length(entrez_up_hsa)<20){
      print("Too few upregulated genes for Cannonical Pathway enrichment (<20)")
      results$cannonicalPathways_up <- "Too few upregulated genes for Motif enrichment (<20)"
    }else{
      cannonicalPathways_up <- enricher(entrez_up_hsa,
                                        TERM2GENE=cannonicalPathway_genes,
                                        universe = universe_mouse2human_Entrez,  
                                        pAdjustMethod = pCorrection,
                                        pvalueCutoff  = pvalueCutoff,
                                        qvalueCutoff = qvalueCutoff)
      
      results$cannonicalPathwaysup <- as.data.frame(cannonicalPathways_up)
      if(nrow(results$cannonicalPathwaysup)<1){
        results$cannonicalPathwaysup_plot <- "No cannonical pathway enrichment for upregulated genes"
      }else{
        results$cannonicalPathwaysup_plot <- clusterProfiler::dotplot(cannonicalPathways_up,
                                                                      showCategory = showMax,
                                                                      font.size= font.size,
                                                                      title = paste("Cannonical pathway  enrichment for genes upregulated in: ",comparison,sep="")
        )
      }
    }
    
    if(length(entrez_down_hsa)<20){
      print("Too few downregulated genes for cannonical pathway  enrichment (<20)")
      results$cannonicalPathways_down <- "Too few downregulated genes for cannonical pathway enrichment (<20)"
    }else{
      cannonicalPathways_down <- enricher(entrez_down_hsa,
                                          TERM2GENE=cannonicalPathway_genes,
                                          universe = universe_mouse2human_Entrez,  
                                          pAdjustMethod = pCorrection,
                                          pvalueCutoff  = pvalueCutoff,
                                          qvalueCutoff = qvalueCutoff)
      
      results$cannonicalPathwaysdown <- as.data.frame(cannonicalPathways_down)
      if(nrow(results$cannonicalPathwaysdown)<1){
        results$cannonicalPathwaysdown_plot <-"No cannonical pathway enrichment for downregulated genes"
      }else{
        results$cannonicalPathwaysdown_plot <- clusterProfiler::dotplot(cannonicalPathways_down,
                                                                        showCategory = showMax,
                                                                        font.size= font.size,
                                                                        title = paste("Cannonical pathway enrichment for genes downregulated in: ",comparison,sep="")
        )
      }
    }
  }
  
  # Motif enrichment
  if("Motifs" %in% GeneSets){
    print("Performing Motif enrichment")
    if(length(entrez_up_hsa)<20){
      print("Too few upregulated genes for Motif enrichment (<20)")
      results$Motif_up <- "Too few upregulated genes for Motif enrichment (<20)"
    }else{
      Motif_up <- enricher(entrez_up_hsa,
                           TERM2GENE=motifs,
                           universe = universe_mouse2human_Entrez,  
                           pAdjustMethod = pCorrection,
                           pvalueCutoff  = pvalueCutoff,
                           qvalueCutoff = qvalueCutoff)
      
      results$Motifup <- as.data.frame(Motif_up)
      if(nrow(results$Motifup)<1){
        results$Motifup_plot <- "No Motif enrichment for upregulated genes"
      }else{
        results$Motifup_plot <- clusterProfiler::dotplot(Motif_up,
                                                         showCategory = showMax,
                                                         font.size= font.size,
                                                         title = paste("Motif enrichment for genes upregulated in: ",comparison,sep="")
        )
      }
    }
    
    if(length(entrez_down_hsa)<20){
      print("Too few downregulated genes for Motif enrichment (<20)")
      results$Motif_down <- "Too few downregulated genes for Motif enrichment (<20)"
    }else{
      Motif_down <- enricher(entrez_down_hsa,
                             TERM2GENE=motifs,
                             universe = universe_mouse2human_Entrez,  
                             pAdjustMethod = pCorrection,
                             pvalueCutoff  = pvalueCutoff,
                             qvalueCutoff = qvalueCutoff)
      
      results$Motifdown <- as.data.frame(Motif_down)
      if(nrow(results$Motifdown)<1){
        results$Motifdown_plot <-"No Motif enrichment for downregulated genes"
      }else{
        results$Motifdown_plot <- clusterProfiler::dotplot(Motif_down,
                                                           showCategory = showMax,
                                                           font.size= font.size,
                                                           title = paste("Motif enrichment for genes downregulated in: ",comparison,sep="")
        )
      }
    }
  }
  
  # Immunosignatures enrichment
  if("ImmunoSignatures" %in% GeneSets){
    print("Performing immunesignature enrichment")
    if(length(entrez_up_hsa)<20){
      print("Too few upregulated genes for Immunosignature enrichment (<20)")
      results$ImmSig_up <- "Too few upregulated genes for Immunosignature enrichment (<20)"
    }else{
      ImmSig_up <- enricher(entrez_up_hsa,
                            TERM2GENE=immuno_genes,
                            universe = universe_mouse2human_Entrez,  
                            pAdjustMethod = pCorrection,
                            pvalueCutoff  = pvalueCutoff,
                            qvalueCutoff = qvalueCutoff)
      
      results$ImmSigup <- as.data.frame(ImmSig_up)
      if(nrow(results$ImmSigup)<1){
        results$ImmSigup_plot <- "No Immunosignature enrichment for upregulated genes"
      }else{
        results$ImmSigup_plot <- clusterProfiler::dotplot(ImmSig_up,
                                                          showCategory = showMax,
                                                          font.size= font.size,
                                                          title = paste("Immunosignature enrichment for genes upregulated in: ",comparison,sep="")
        )
      }
    }
    if(length(entrez_down_hsa)<20){
      print("Too few downregulated genes for Immunosignature enrichment (<20)")
      results$ImmSig_down <- "Too few downregulated genes for Immunosignature enrichment (<20)"
    }else{
      ImmSig_down <- enricher(entrez_down_hsa,
                              TERM2GENE=immuno_genes,
                              universe = universe_mouse2human_Entrez,  
                              pAdjustMethod = pCorrection,
                              pvalueCutoff  = pvalueCutoff,
                              qvalueCutoff = qvalueCutoff)
      
      results$ImmSigdown <- as.data.frame(ImmSig_down)
      if(nrow(results$ImmSigdown)<1){
        results$ImmSigdown_plot <- "No Immunosignature enrichment for downregulated genes"
      }else{
        results$ImmSigdown_plot <- clusterProfiler::dotplot(ImmSig_down,
                                                            showCategory = showMax,
                                                            font.size= font.size,
                                                            title = paste("Immunosignature enrichment for genes downregulated in: ",comparison,sep="")
        )
      }
    }
  }
  results
}

#' GO & KEGG enrichment across comparisons
#' @export
compareGSEA <- function(comparisons,
                        organism,
                        GeneSets =c("GO","KEGG"),
                        ontology= "BP",
                        pCorrection = "bonferroni", # choose the p-value adjustment method
                        pvalueCutoff = 0.05, # set the unadj. or adj. p-value cutoff (depending on correction method)
                        qvalueCutoff = 0.05, # set the q-value cutoff (FDR corrected)
                        showMax = 20){

  if(organism == "mouse") {
    OrgDb = org.Mm.eg.db
  } else if(organism == "human"){
    OrgDb = org.Hs.eg.db
  } else {stop("Wrong Organism. Select mouse or human.")}

  ENTREZlist <-  list()
  for(i in 1:length(comparisons)){
    res <- DEresults[names(DEresults) %in% comparisons]
    DE_up <- as.data.frame(res[[i]]@DE_genes$up_regulated_Genes)$SYMBOL
    entrez_up <- bitr(DE_up, fromType = "SYMBOL", toType="ENTREZID", OrgDb=OrgDb)$ENTREZID
    DE_down <- as.data.frame(res[[i]]@DE_genes$down_regulated_Genes)$SYMBOL
    entrez_down <- bitr(DE_down, fromType = "SYMBOL", toType="ENTREZID", OrgDb=OrgDb)$ENTREZID
    x <- setNames(list(entrez_up, entrez_down),
                  c(paste(names(res[i]),"_up",sep=""),
                    paste(names(res[i]),"_down",sep="")))
    ENTREZlist <- c(ENTREZlist,x)
  }

  list <- list()

  # Compare the Clusters regarding their GO enrichment
  if("GO" %in% GeneSets){
    print("Performing GO enrichment")
    CompareClusters_GO <- compareCluster(geneCluster = ENTREZlist,
                                         fun = "enrichGO",
                                         universe = universe_Entrez,
                                         OrgDb = OrgDb,
                                         ont = ontology,
                                         pvalueCutoff  = pvalueCutoff,
                                         pAdjustMethod = pCorrection,
                                         qvalueCutoff  = pvalueCutoff,
                                         readable      = T)
    list$GOresults <- as.data.frame(CompareClusters_GO)
    list$GOplot <- clusterProfiler::dotplot(CompareClusters_GO, showCategory = showMax, by = "geneRatio", font.size=10)
  }

  if("KEGG" %in% GeneSets){
    print("Performing KEGG enrichment")

    if(organism == "mouse"){org = "mmu"}
    if(organism == "human"){org = "hsa"}

    # Compare the Clusters regarding their KEGG enrichment
    CompareClusters_KEGG <- compareCluster(geneCluster = ENTREZlist,
                                           fun = "enrichKEGG",
                                           universe = universe_Entrez,
                                           organism = org,
                                           pvalueCutoff  = pvalueCutoff,
                                           pAdjustMethod = pCorrection,
                                           qvalueCutoff  = pvalueCutoff)
    list$KEGGresults <- as.data.frame(CompareClusters_KEGG)
    list$KEGGplot <- clusterProfiler::dotplot(CompareClusters_KEGG, showCategory = showMax, by = "geneRatio", font.size=10)
  }
  list
}

