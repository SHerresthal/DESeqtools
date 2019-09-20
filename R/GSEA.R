#' GSEA function
#' @export

GSEA <-  function(comparison,
                  organism,
                  DE_results = DEresults,
                  GeneSets =c("GO","KEGG","DO","Hallmark","canonicalPathways","Motifs","ImmunoSignatures"),
                  GOntology = "BP",
                  pCorrection = "bonferroni", # choose the p-value adjustment method
                  pvalueCutoff = 0.05, # set the unadj. or adj. p-value cutoff (depending on correction method)
                  qvalueCutoff = 0.05 # set the q-value cutoff (FDR corrected)
){

  results <- list()

  if(organism == "mouse") {
    OrgDb = org.Mm.eg.db
  } else if(organism == "human"){
    OrgDb = org.Hs.eg.db
  } else {print("Wrong Organism. Select mouse or human.")}

  res <- DE_results[[comparison]]
  DE_up <- as.data.frame(res@DE_genes$up_regulated_Genes)$SYMBOL
  entrez_up <- bitr(DE_up, fromType = "SYMBOL", toType="ENTREZID", OrgDb=OrgDb)$ENTREZID
  DE_down <- as.data.frame(res@DE_genes$down_regulated_Genes)$SYMBOL
  entrez_down <- bitr(DE_down, fromType = "SYMBOL", toType="ENTREZID", OrgDb=OrgDb)$ENTREZID

  # GO enrichment ###############################
  if("GO" %in% GeneSets){
    print("Performing GO enrichment")
    if(length(entrez_up)<20){
      print("Too few upregulated genes for GO enrichment (<20)")
      results$GOup <- "Too few upregulated genes for GO enrichment (<20)"
    }else{
      results$GOup <- as.data.frame(enrichGO(gene = entrez_up,
                                             universe = universe_Entrez,
                                             OrgDb = OrgDb,
                                             ont = GOntology,
                                             pAdjustMethod = pCorrection,
                                             pvalueCutoff  = pvalueCutoff,
                                             qvalueCutoff  = qvalueCutoff,
                                             readable      = T))

      if(nrow(results$GOup)>0){results$GOup$Enrichment <- paste("GO enrichment for genes upregulated in ",comparison,sep="")}
    }
    if(length(entrez_down)<20){
      print("Too few downregulated genes for GO enrichment (<20)")
      results$GOdown <- "Too few downregulated genes for GO enrichment (<20)"
    }else{
      results$GOdown <- as.data.frame(enrichGO(gene = entrez_down,
                                               universe = universe_Entrez,
                                               OrgDb = OrgDb,
                                               ont = GOntology,
                                               pAdjustMethod = pCorrection,
                                               pvalueCutoff  = pvalueCutoff,
                                               qvalueCutoff  = qvalueCutoff,
                                               readable      = T))
      if(nrow(results$GOdown)>0){results$GOdown$Enrichment <- paste("GO enrichment for genes downregulated in ",comparison,sep="")}
    }
  }

  # KEGG enrichment ##########################################
  if("KEGG" %in% GeneSets){
    print("Performing KEGG enrichment")

    if(organism == "mouse") {org = "mmu"}
    if(organism == "human"){org = "hsa"}

    if(length(entrez_up)<20){
      print("Too few upregulated genes for KEGG enrichment (<20)")
      results$KEGGup <- "Too few upregulated genes for KEGG enrichment (<20)"
    }else{
      results$KEGGup <- as.data.frame(enrichKEGG(gene = entrez_up,
                                                 organism = org,
                                                 universe = universe_Entrez,
                                                 pAdjustMethod = pCorrection,
                                                 pvalueCutoff  = pvalueCutoff,
                                                 qvalueCutoff = qvalueCutoff))
      if(nrow(results$KEGGup)>0){results$KEGGup$Enrichment <- paste("KEGG enrichment for genes upregulated in ",comparison,sep="")}
    }
    if(length(entrez_down)<20){
      print("Too few downregulated genes for KEGG enrichment (<20)")
      results$KEGGdown <- "Too few downregulated genes for KEGG enrichment (<20)"
    } else{
      results$KEGGdown <- as.data.frame(enrichKEGG(gene = entrez_down,
                                                   organism = org,
                                                   universe = universe_Entrez,
                                                   pAdjustMethod = pCorrection,
                                                   pvalueCutoff  = pvalueCutoff,
                                                   qvalueCutoff = qvalueCutoff))
      if(nrow(results$KEGGdown)>0){results$KEGGdown$Enrichment <- paste("KEGG enrichment for genes downregulated in ",comparison,sep="")}
    }
  }

  if("Hallmark" %in% GeneSets |
     "DO" %in% GeneSets |
     "canonicalPathways" %in% GeneSets|
     "ImmunoSignatures" %in% GeneSets |
     "Motifs" %in% GeneSets){
    if(organism == "mouse"){

    entrez_up_hsa <- as.character(getLDS(attributes = c("mgi_symbol"),
                                         filters = "mgi_symbol",
                                         values = DE_up,
                                         mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl"),
                                         attributesL = c("entrezgene_id"),
                                         martL = useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
                                         uniqueRows=T)[,2])
    entrez_down_hsa <- getLDS(attributes = c("mgi_symbol"),
                              filters = "mgi_symbol",
                              values = DE_down,
                              mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl"),
                              attributesL = c("entrezgene_id"),
                              martL = useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
                              uniqueRows=T)[,2]

    } else if(organism == "human"){
      entrez_up_hsa <- entrez_up
      entrez_down_hsa <- entrez_down
    }
    }

  # DO enrichment ########################################
  if("DO" %in% GeneSets){
    print("Performing Disease Ontology enrichment")

    if(length(entrez_up)<20){
      print("Too few upregulated genes for DO enrichment (<20)")
      results$DOup <- "Too few upregulated genes for DO enrichment (<20)"
    }else{
      results$DOup <- as.data.frame(enrichDO(gene = entrez_up_hsa,
                                             universe = universe_mouse2human_Entrez,
                                             pAdjustMethod = pCorrection,
                                             pvalueCutoff  = pvalueCutoff,
                                             qvalueCutoff = qvalueCutoff,
                                             minGSSize     = 5,
                                             maxGSSize     = 500,
                                             readable=TRUE))
      if(nrow(results$DOup)>0){results$DOup$Enrichment <- paste("DO enrichment for genes upregulated in ",comparison,sep="")}
    }
    if(length(entrez_down)<20){
      print("Too few downregulated genes for DO enrichment (<20)")
      results$DOdown <- "Too few downregulated genes for DO enrichment (<20)"
    } else{
      results$DOdown <- as.data.frame(enrichDO(gene = entrez_down_hsa,
                                               universe = universe_mouse2human_Entrez,
                                               pAdjustMethod = pCorrection,
                                               pvalueCutoff  = pvalueCutoff,
                                               qvalueCutoff = qvalueCutoff,
                                               minGSSize     = 5,
                                               maxGSSize     = 500,
                                               readable=TRUE))
      if(nrow(results$DOdown)>0){results$DOdown$Enrichment <- paste("DO enrichment for genes downregulated in ",comparison,sep="")}
    }
  }

  # Hallmark enrichment ################################
  if("Hallmark" %in% GeneSets){
    print("Performing Hallmark enrichment")
    if(length(entrez_up_hsa)<20){
      print("Too few upregulated genes for Hallmark enrichment (<20)")
      results$Hallmarkup <- "Too few upregulated genes for Hallmark enrichment (<20)"
    }else{
      results$HALLMARKup <- as.data.frame(enricher(entrez_up_hsa,
                                                   TERM2GENE=hallmark_genes,
                                                   universe = universe_mouse2human_Entrez,
                                                   pAdjustMethod = pCorrection,
                                                   pvalueCutoff  = pvalueCutoff,
                                                   qvalueCutoff = qvalueCutoff))
      if(nrow(results$HALLMARKup)>0){results$HALLMARKup$Enrichment <- paste("HALLMARK enrichment for genes upregulated in ",comparison,sep="")}
    }
    if(length(entrez_down_hsa)<20){
      print("Too few downregulated genes for Hallmark enrichment (<20)")
      results$Hallmarkdown <- "Too few downregulated genes for Hallmark enrichment (<20)"
    }else{
      results$HALLMARKdown <- as.data.frame(enricher(entrez_down_hsa,
                                                     TERM2GENE=hallmark_genes,
                                                     universe = universe_mouse2human_Entrez,
                                                     pAdjustMethod = pCorrection,
                                                     pvalueCutoff  = pvalueCutoff,
                                                     qvalueCutoff = qvalueCutoff))
      if(nrow(results$HALLMARKdown)>0){results$HALLMARKdown$Enrichment <- paste("HALLMARK enrichment for genes downregulated in ",comparison,sep="")}
    }
  }

  # Canonical Pathway enrichment #############################
  if("canonicalPathways" %in% GeneSets){
    print("Performing Canonical Pathway (C2) enrichment")
    if(length(entrez_up_hsa)<20){
      print("Too few upregulated genes for Canonical Pathway enrichment (<20)")
      results$canonicalPathwaysup <- "Too few upregulated genes for Motif enrichment (<20)"
    }else{
      results$canonicalPathwaysup <- as.data.frame(enricher(entrez_up_hsa,
                                                             TERM2GENE=canonicalPathway_genes,
                                                             universe = universe_mouse2human_Entrez,
                                                             pAdjustMethod = pCorrection,
                                                             pvalueCutoff  = pvalueCutoff,
                                                             qvalueCutoff = qvalueCutoff))
      if(nrow(results$canonicalPathwaysup)>0){results$canonicalPathwaysup$Enrichment <- paste("Canonical pathway enrichment for genes upregulated in ",comparison,sep="")}
    }
    if(length(entrez_down_hsa)<20){
      print("Too few downregulated genes for canonical pathway  enrichment (<20)")
      results$canonicalPathwaysdown <- "Too few downregulated genes for canonical pathway enrichment (<20)"
    }else{
      results$canonicalPathwaysdown <- as.data.frame(enricher(entrez_down_hsa,
                                                               TERM2GENE=canonicalPathway_genes,
                                                               universe = universe_mouse2human_Entrez,
                                                               pAdjustMethod = pCorrection,
                                                               pvalueCutoff  = pvalueCutoff,
                                                               qvalueCutoff = qvalueCutoff))
      if(nrow(results$canonicalPathwaysdown)>0){results$canonicalPathwaysdown$Enrichment <- paste("Canonical pathway enrichment for genes downregulated in ",comparison,sep="")}
    }
  }

  # Motif enrichment
  if("Motifs" %in% GeneSets){
    print("Performing Motif enrichment")
    if(length(entrez_up_hsa)<20){
      print("Too few upregulated genes for Motif enrichment (<20)")
      results$Motifup <- "Too few upregulated genes for Motif enrichment (<20)"
    }else{
      results$Motifup <- as.data.frame(enricher(entrez_up_hsa,
                                                TERM2GENE=motifs,
                                                universe = universe_mouse2human_Entrez,
                                                pAdjustMethod = pCorrection,
                                                pvalueCutoff  = pvalueCutoff,
                                                qvalueCutoff = qvalueCutoff))
      if(nrow(results$Motifup)>0){results$Motifup$Enrichment <- paste("TF binding motif enrichment for genes upregulated in ",comparison,sep="")}
    }
    if(length(entrez_down_hsa)<20){
      print("Too few downregulated genes for Motif enrichment (<20)")
      results$Motifdown <- "Too few downregulated genes for Motif enrichment (<20)"
    }else{
      results$Motifdown <- as.data.frame(enricher(entrez_down_hsa,
                                                  TERM2GENE=motifs,
                                                  universe = universe_mouse2human_Entrez,
                                                  pAdjustMethod = pCorrection,
                                                  pvalueCutoff  = pvalueCutoff,
                                                  qvalueCutoff = qvalueCutoff))
      if(nrow(results$Motifdown)>0){results$Motifdown$Enrichment <- paste("TF binding motif enrichment for genes downregulated in ",comparison,sep="")}
    }
  }

  # Immunosignatures enrichment
  if("ImmunoSignatures" %in% GeneSets){
    print("Performing Immunosignature enrichment")
    if(length(entrez_up_hsa)<20){
      print("Too few upregulated genes for Immunosignature enrichment (<20)")
      results$ImmSigup <- "Too few upregulated genes for Immunosignature enrichment (<20)"
    }else{
      results$ImmSigup <- as.data.frame(enricher(entrez_up_hsa,
                                                 TERM2GENE=immuno_genes,
                                                 universe = universe_mouse2human_Entrez,
                                                 pAdjustMethod = pCorrection,
                                                 pvalueCutoff  = pvalueCutoff,
                                                 qvalueCutoff = qvalueCutoff))
      if(nrow(results$ImmSigup)>0){results$ImmSigup$Enrichment <- paste("Immunosignature enrichment for genes upregulated in ",comparison,sep="")}
    }
    if(length(entrez_down_hsa)<20){
      print("Too few downregulated genes for Immunosignature enrichment (<20)")
      results$ImmSigdown <- "Too few downregulated genes for Immunosignature enrichment (<20)"
    }else{
      results$ImmSigdown <- as.data.frame(enricher(entrez_down_hsa,
                                                   TERM2GENE=immuno_genes,
                                                   universe = universe_mouse2human_Entrez,
                                                   pAdjustMethod = pCorrection,
                                                   pvalueCutoff  = pvalueCutoff,
                                                   qvalueCutoff = qvalueCutoff))
      if(nrow(results$ImmSigdown)>0){results$ImmSigdown$Enrichment <- paste("Immunosignature enrichment for genes downregulated in ",comparison,sep="")}
    }
  }
  results
}

#' GO & KEGG enrichment across comparisons
#' @export
compareGSEA <- function(comparisons,
                        DE_results = DEresults,
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
    res <- DE_results[names(DE_results) %in% comparisons]
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

#' GSEA dotplot
#' @export
dotplotGSEA <- function(x,
                        show=25,
                        font.size=10,
                        title.size=10,
                        title.width=100,
                        order="count"){
  if(nrow(x)<1){
    print("No enrichment found.")
  }else{
    x <- if(nrow(x)>show){x[c(1:show),]}else{x}
    if(order=="padj"){
      x <- x[order(x$Count,decreasing=FALSE),]
      x$GeneRatio <- factor(x$GeneRatio, levels = unique(x$GeneRatio))
      x <- x[order(x$p.adjust,decreasing=TRUE),]
      x$Description <- factor(x$Description, levels = unique(x$Description))
    }
    if(order=="count"){
      x <- x[order(x$Count,decreasing=FALSE),]
      x$Description <- factor(x$Description, levels = unique(x$Description))
      x$GeneRatio <- factor(x$GeneRatio, levels = unique(x$GeneRatio))
    }
    ggplot(x, aes(x = GeneRatio, y = Description, color = p.adjust)) +
      geom_point(aes(size = Count)) +
      scale_colour_gradientn(colours=c('red',
                                       'orange',
                                       'darkblue',
                                       'darkblue'),
                             limits=c(0,1),
                             values   = c(0,0.05,0.2,0.5,1),
                             breaks   = c(0.05,0.2,1),
                             labels = format(c(0.05,0.2,1))) +
      ylab(NULL) +
      ggtitle(paste(strwrap(unique(x$Enrichment), width=title.width), collapse = "\n"))+
      theme_bw() +
      theme(text = element_text(size=font.size),
            plot.title = element_text(size=title.size))
  }
}


