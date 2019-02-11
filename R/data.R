#' List of GO terms (human) - Part 1
#'
#' A dataset containing GO term description, the definition and the domain for human gene IDs.
#' Needs to be combined with GO_hs_2 to get the full dataframe.
#'
#' @name GO_hs_1
#' @references original file GO_hg38p12_ensembl181121.txt

#' @keywords gene ontology
#' @format A data frame
#'  \describe{
#'   \item{GENEID}{Ontology term}
#'   \item{SYMBOL}{Entrez ID}
#'   \item{TERM}{GO term description}
#'   \item{DEFINITION}{GO term definition}
#'   \item{DOMAIN}{GO domain. biological_process, cellular_component or molecular_function}
#'   }
"GO_hs_1"

#' List of GO terms (human) - Part 2
#'
#' A dataset containing GO term description, the definition and the domain for human gene IDs.
#' Needs to be combined with GO_hs_1 to get the full dataframe.
#'
#' @name GO_hs_2
#' @references original file GO_hg38p12_ensembl181121.txt

#' @keywords gene ontology
#' @format A data frame
#'  \describe{
#'   \item{GENEID}{Ontology term}
#'   \item{SYMBOL}{Entrez ID}
#'   \item{TERM}{GO term description}
#'   \item{DEFINITION}{GO term definition}
#'   \item{DOMAIN}{GO domain. biological_process, cellular_component or molecular_function}
#'   }
"GO_hs_2"


#' List of GO terms (mouse) - Part 1
#'
#' A dataset containing GO term description, the definition and the domain for mouse gene IDs.
#' Needs to be combined with GO_mm_2 to get the full dataframe.
#'
#' @name GO_mm_1
#' @references original file GO_mm38p12_ensembl181121.txt
#' @keywords gene ontology
#' @format A data frame
#'  \describe{
#'   \item{GENEID}{Ontology term}
#'   \item{SYMBOL}{Entrez ID}
#'   \item{TERM}{GO term description}
#'   \item{DEFINITION}{GO term definition}
#'   \item{DOMAIN}{GO domain. biological_process, cellular_component or molecular_function}
#'   }
#"GO_mm_1"


#' List of GO terms (mouse) - Part 2
#'
#' A dataset containing GO term description, the definition and the domain for mouse gene IDs.
#' Needs to be combined with GO_mm_1 to get the full dataframe.
#'
#' @name GO_mm_2
#' @references original file GO_mm38p12_ensembl181121.txt
#' @keywords gene ontology
#' @format A data frame
#'  \describe{
#'   \item{GENEID}{Ontology term}
#'   \item{SYMBOL}{Entrez ID}
#'   \item{TERM}{GO term description}
#'   \item{DEFINITION}{GO term definition}
#'   \item{DOMAIN}{GO domain. biological_process, cellular_component or molecular_function}
#'   }
#"GO_mm_2"


#' Biomart annotation (human)
#'
#' Gene annotations as downloaded brom biomart database on 2019-02-19, Genome version GRCh38.p12.
#'
#' @name biomart_human
#' @references \url{https://www.ensembl.org/biomart/}
#' @keywords annotation
#' @format A data frame with 64914 obs. of  7 variables. The data frame contains the gene ID, the gene name, gene description, the chromosome scaffold name, strand, gene start and stop basepairs.
"biomart_human"

#' Biomart annotation (mouse)
#'
#' Gene annotations as downloaded brom biomart database on 2018-09-14, Genome version GRCm38.p6.
#'
#' @name biomart_human
#' @references \url{https://www.ensembl.org/biomart/}
#' @keywords annotation
#' @format A data frame with 55029 obs. of  7 variables. The data frame contains the gene ID, the gene name, gene description, the chromosome scaffold name, strand, gene start and stop basepairs.
"biomart_mouse"


#' Canonical pathway genes
#'
#' List with KEGG canonical pathways from the MiSigDB gene set.
#' @name canonicalPathway_genes
#' @references This datasets corresponds to c2.cp.v6.2.entrez.gmt \url{http://software.broadinstitute.org/gsea/msigdb/collections.jsp}
#' @format A data frame
#'   \describe{
#'   \item{ont}{Pathway term}
#'   \item{gene}{Entrez ID}
#'   }
"canonicalPathway_genes"

#' Immuno genes
#'
#' @name immuno_genes
#' @references This datasets corresponds to the GMT file c7.all.v6.2.entrez.gmt \url{http://software.broadinstitute.org/gsea/msigdb/collections.jsp}
#' @format A data frame
#'  \describe{
#'   \item{ont}{Pathway}
#'   \item{gene}{Entrez ID}
#'   }
"immuno_genes"

#' Hallmark genes
#'
#' @name hallmark_genes
#' @references This datasets corresponds to the GMT file h.all.v6.2.entrez.gmt \url{http://software.broadinstitute.org/gsea/msigdb/collections.jsp}
#' @keywords pathway genes
#' @format A data frame
#'  \describe{
#'   \item{ont}{Pathway}
#'   \item{gene}{Entrez ID}
#'   }
"hallmark_genes"

#' Motifs
#'
#' @name Motifs
#' @references This datasets corresponds to the GMT file c3.all.v6.2.entrez.gmt \url{http://software.broadinstitute.org/gsea/msigdb/collections.jsp}
#' @keywords pathway genes
#' @format A data frame
"motifs"


#' KEGG terms human
#'
#' @name KEGG_hs
#' @references original filename KEGG_hg38_clusterProfiler181121.txt
#' @format A data frame
"KEGG_hs"

#' KEGG terms mouse
#'
#' KEGG terms and corresponding genes.
#'
#' @name KEGG_mm
#' @references Original filename KEGG_mm10_clusterProfiler181121.txt
#' @format A data frame
"KEGG_mm"


#' TF list human
#'
#' @name TFlist_hs
#' @references \url{}
#' @format A data frame .
"TFlist_hs"

#' TF list mouse
#'
#' @name TFlist_mm
#' @references \url{}
#' @format A data frame
"TFlist_mm"


#' TX annotation human v27
#'
#' This gene annotation file is used
#' 1) to map the Ensembl transcript IDs to Ensembl gene IDs during the tximport function and
#' 2) annotate the Ensembl IDs with additional information such as gene symbol or type.
#' For consistency, this file should have been produced from the .gtf file used for building the kallisto index. It needs to consist
#' of four columns: Gene ENSEMBL ID, Transcript ID, Gene Symbol, Gene Type.
#'
#' @name tx_annotation_v27_human
#' @references This file corresponds to ID2SYMBOL_gencode_v27_transcript.txt, which has been generated from the GTF file gencode.v27.primary_assembly.annotation.gtf, as downloaded from \url{https://www.gencodegenes.org/human/release_27.html}.
#' @keywords annotation
#' @format A data frame with 200468 obs. of  4 variables:
#'  \describe{
#'   \item{GENEID}{Ensembl Gene ID}
#'   \item{TXNAME}{Transcript ID}
#'   \item{SYMBOL}{Gene Symbol}
#'   \item{GENETYPE}{Genetype description. See \url{https://www.gencodegenes.org/pages/biotypes.html} for a list of all gene types.
#'   }
#'   }
"tx_annotation_v27_human"


#' TX annotation mouse vM13
#'
#' This gene annotation file is used
#' 1) to map the Ensembl transcript IDs to Ensembl gene IDs during the tximport function and
#' 2) annotate the Ensembl IDs with additional information such as gene symbol or type.
#' For consistency, this file should have been produced from the .gtf file used for building the kallisto index. It needs to consist
#' of four columns: Gene ENSEMBL ID, Transcript ID, Gene Symbol, Gene Type.
#'
#' @name tx_annotation_vM13_mouse
#' @references This file corresponds to ID2SYMBOL_gencode_vM13_transcript.txt, which has been generated from the GTF file gencode.vM13.primary_assembly.annotation.gtf, as downloaded from \url{ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M13/gencode.vM13.primary_assembly.annotation.gtf.gz}.
#' @keywords annotation
#' @format A data frame with 125570 obs. of  4 variables:
#'  \describe{
#'   \item{GENEID}{Ensembl Gene ID}
#'   \item{TXNAME}{Transcript ID}
#'   \item{SYMBOL}{Gene Symbol}
#'   \item{GENETYPE}{Genetype description. See \url{https://www.gencodegenes.org/pages/biotypes.html} for a list of all gene types.
#'   }
#'   }
"tx_annotation_vM13_mouse"


#' TX annotation mouse vM16
#'
#' This gene annotation file is used
#' 1) to map the Ensembl transcript IDs to Ensembl gene IDs during the tximport function and
#' 2) annotate the Ensembl IDs with additional information such as gene symbol or type.
#' For consistency, this file should have been produced from the .gtf file used for building the kallisto index. It needs to consist
#' of four columns: Gene ENSEMBL ID, Transcript ID, Gene Symbol, Gene Type.
#' @name tx_annotation_vM16_mouse
#' @references This file corresponds to ID2SYMBOL_gencode_vM16_transcript.txt, which has been generated from the GTF file gencode.vM16.primary_assembly.annotation.gtf, as downloaded from \url{ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M16/gencode.vM16.primary_assembly.annotation.gtf.gz}.
#' @keywords annotation
#' @format A data frame with 133944 obs. of  4 variables:
#'  \describe{
#'   \item{GENEID}{Ensembl Gene ID}
#'   \item{TXNAME}{Transcript ID}
#'   \item{SYMBOL}{Gene Symbol}
#'   \item{GENETYPE}{Genetype description. See \url{https://www.gencodegenes.org/pages/biotypes.html} for a list of all gene types.
#'   }
#'   }
"tx_annotation_vM16_mouse"

#' Example sample table referring to the example dataset
#'
#' @name example_sampletable
#' @references \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92852}
#' @format
"example_sampletable"

