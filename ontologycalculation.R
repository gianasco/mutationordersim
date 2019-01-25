comomineR <- function(gseries.1 = "sample1", gplatform.1 = "plat1", 
                      gsms.1 = "XXXXXXXXXXXXXXX00000000001111111XXXXXX", 
                      affyLib.1 = "hgu133a2.db",
                      gseries.2 = "sample2", gplatform.2 = "plat1", 
                      gsms.2 = "1111111111111110000000000XXXXXXXXXXXXX",
                      affyLib.2 = "hgu133a2.db",
                      top_genes = 500, top_GO_terms = 50, GO.node.size = 5) {
  
  # load appropriate libraries
  library(Biobase)
  library(GEOquery)
  library(limma)
  library(topGO)
  
  #########     Load Series 1     #########
  # load series and platform data
  gset.1 <- getGEO(gseries.1, GSEMatrix =TRUE, AnnotGPL=TRUE)
  if (length(gset.1) > 1) {
    idx <- grep(gplatform.1, attr(gset, "names")) 
  } else { 
    idx <- 1
  }
  gset.1 <- gset.1[[idx]]
  
  # make proper column names
  fvarLabels(gset.1) <- make.names(fvarLabels(gset.1))
  
  
  #########     Load Series 2     #########
  # load series and platform data
  gset.2 <- getGEO(gseries.2, GSEMatrix =TRUE, AnnotGPL=TRUE)
  if (length(gset.2) > 1) {
    idx <- grep(gplatform.2, attr(gset, "names")) 
  } else { 
    idx <- 1
  }
  gset.2 <- gset.2[[idx]]
  
  # make proper column names
  fvarLabels(gset.2) <- make.names(fvarLabels(gset.2))
  
  
  #########     Series Preprocessing     #########
  # find common features/genes between the two GEO Series
  features.1 = featureNames(gset.1)
  features.2 = featureNames(gset.2)
  #common.features = features.1[1:10]
  common.features = intersect(features.1, features.2)
  
  cat("Common features: ", length(common.features))
  
  gset.1 <- gset.1[common.features,]
  gset.2 <- gset.2[common.features,]
  
  # group names for all samples of the first set
  sml.1 <- c()
  for (i in 1:nchar(gsms.1)) { 
    sml.1[i] <- substr(gsms.1,i,i) 
  }
  
  # eliminate samples marked as "X"
  sel.1 <- which(sml.1 != "X")
  sml.1 <- sml.1[sel.1]
  gset.1 <- gset.1[ ,sel.1]
  
  # group names for all samples of the second set
  sml.2 <- c()
  for (i in 1:nchar(gsms.2)) { 
    sml.2[i] <- substr(gsms.2,i,i) 
  }
  
  # eliminate samples marked as "X"
  sel.2 <- which(sml.2 != "X")
  sml.2 <- sml.2[sel.2]
  gset.2 <- gset.2[ ,sel.2]
  
  # log2 transform for first series
  ex.1 <- exprs(gset.1)
  qx.1 <- as.numeric(quantile(ex.1, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC.1 <- (qx.1[5] > 100) ||
    (qx.1[6]-qx.1[1] > 50 && qx.1[2] > 0) ||
    (qx.1[2] > 0 && qx.1[2] < 1 && qx.1[4] > 1 && qx.1[4] < 2)
  if (LogC.1) { ex.1[which(ex.1 <= 0)] <- NaN
  exprs(gset.1) <- log2(ex.1) }
  
  # log2 transform for second series
  ex.2 <- exprs(gset.2)
  qx.2 <- as.numeric(quantile(ex.2, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC.2 <- (qx.2[5] > 100) ||
    (qx.2[6]-qx.2[1] > 50 && qx.2[2] > 0) ||
    (qx.2[2] > 0 && qx.2[2] < 1 && qx.2[4] > 1 && qx.2[4] < 2)
  if (LogC.2) { ex.2[which(ex.2 <= 0)] <- NaN
  exprs(gset.2) <- log2(ex.2) }
  
  #########     Set Up and Analysis     #########
  # set up the data and proceed with analyzing the first Series
  sml.1 <- paste("G", sml.1, sep="")    # set group names
  fl.1 <- as.factor(sml.1)
  gset.1$description <- fl.1
  design.1 <- model.matrix(~ description + 0, gset.1)
  colnames(design.1) <- levels(fl.1)
  fit.1 <- lmFit(gset.1, design.1)
  
  # set up the data and proceed with analyzing the second Series
  sml.2 <- paste("G", sml.2, sep="")    # set group names
  fl.2 <- as.factor(sml.2)
  gset.2$description <- fl.2
  design.2 <- model.matrix(~ description + 0, gset.2)
  colnames(design.2) <- levels(fl.2)
  fit.2 <- lmFit(gset.2, design.2)
  
  # Analysis: fit one model per Series and then get 
  #           the top_genes genes for each Series
  cont.matrix.1 <- makeContrasts(G1-G0, levels=design.1)
  cont.matrix.2 <- makeContrasts(G1-G0, levels=design.2)
  fit2.1 <- contrasts.fit(fit.1, cont.matrix.1)
  fit2.2 <- contrasts.fit(fit.2, cont.matrix.2)
  fit2.1 <- eBayes(fit2.1, 0.01)
  fit2.2 <- eBayes(fit2.2, 0.01)
  tT.1 <- topTable(fit2.1, adjust="fdr", 
                   sort.by="B", number=top_genes)
  tT.2 <- topTable(fit2.2, adjust="fdr", 
                   sort.by="B", number=top_genes)
  
  # Find the number of the common genes
  tT.1 <- subset(tT.1, select=c("ID","adj.P.Val",
                                "P.Value","t","B",
                                "logFC","Gene.symbol",
                                "Gene.title"))
  tT.2 <- subset(tT.2, select=c("ID","adj.P.Val",
                                "P.Value","t","B",
                                "logFC","Gene.symbol",
                                "Gene.title"))
  common.genes <- intersect(tT.1$ID, tT.2$ID)
  
  
  #########     TopGO     #########
  
  source("http://bioconductor.org/biocLite.R")
  #library(package = affyLib.1, character.only = TRUE)
  #library(package = affyLib.2, character.only = TRUE)
  biocLite(affyLib.1, suppressUpdates = TRUE)
  biocLite(affyLib.2, suppressUpdates = TRUE)
  
  allCommonGenes <- rownames(gset.1)
  
  
  #########     Analysis with GO for Biological Processes     #########
  
  # a list with all the genes where 1 indicates an interesting gene
  geneList.1 <- factor(as.integer(allCommonGenes %in% tT.1$ID))
  names(geneList.1) <- allCommonGenes
  BPGOdata.1 <- new("topGOdata", ontology = "BP", 
                    allGenes = geneList.1,
                    nodeSize = GO.node.size,
                    annot = annFUN.db, affyLib = affyLib.1)
  
  BP.res.fisher.1 <- runTest(BPGOdata.1, algorithm = "classic", 
                             statistic = "fisher")
  allRes.1 <- GenTable(BPGOdata.1, classicFisher = BP.res.fisher.1, 
                       orderBy = "classicFisher", 
                       topNodes = top_GO_terms)
  
  geneList.2 <- factor(as.integer(allCommonGenes %in% tT.2$ID))
  names(geneList.2) <- allCommonGenes
  BPGOdata.2 <- new("topGOdata", ontology = "BP",
                    allGenes = geneList.2,
                    nodeSize = GO.node.size,
                    annot = annFUN.db, affyLib = affyLib.2)
  
  BP.res.fisher.2 <- runTest(BPGOdata.2, algorithm = "classic", 
                             statistic = "fisher")
  allRes.2 <- GenTable(BPGOdata.2, classicFisher = BP.res.fisher.2, 
                       orderBy = "classicFisher", 
                       topNodes = top_GO_terms)
  
  # count the number of common significant ontology terms
  common.BP.GO.terms = intersect(allRes.1$GO.ID, allRes.2$GO.ID)
  
  
  #########     Analysis with GO for Molecular Function     #########
  
  MFGOdata.1 <- new("topGOdata", ontology = "MF", 
                    allGenes = geneList.1,
                    nodeSize = GO.node.size,
                    annot = annFUN.db, affyLib = affyLib.1)
  
  MF.res.fisher.1 <- runTest(MFGOdata.1, algorithm = "classic", 
                             statistic = "fisher")
  allRes.1 <- GenTable(MFGOdata.1, classicFisher = MF.res.fisher.1, 
                       orderBy = "classicFisher", 
                       topNodes = top_GO_terms)
  
  MFGOdata.2 <- new("topGOdata", ontology = "MF",
                    allGenes = geneList.2,
                    nodeSize = GO.node.size,
                    annot = annFUN.db, affyLib = affyLib.2)
  
  MF.res.fisher.2 <- runTest(MFGOdata.2, algorithm = "classic", 
                             statistic = "fisher")
  allRes.2 <- GenTable(MFGOdata.2, classicFisher = MF.res.fisher.2, 
                       orderBy = "classicFisher", 
                       topNodes = top_GO_terms)
  
  # count the number of common significant ontology terms
  common.MF.GO.terms = intersect(allRes.1$GO.ID, allRes.2$GO.ID)
  
  
  #########     Analysis with GO for Cellural Component     #########
  
  CCGOdata.1 <- new("topGOdata", ontology = "CC", 
                    allGenes = geneList.1,
                    nodeSize = GO.node.size,
                    annot = annFUN.db, affyLib = affyLib.1)
  
  CC.res.fisher.1 <- runTest(CCGOdata.1, algorithm = "classic", 
                             statistic = "fisher")
  allRes.1 <- GenTable(CCGOdata.1, classicFisher = CC.res.fisher.1, 
                       orderBy = "classicFisher", 
                       topNodes = top_GO_terms)
  
  CCGOdata.2 <- new("topGOdata", ontology = "CC",
                    allGenes = geneList.2,
                    nodeSize = GO.node.size,
                    annot = annFUN.db, affyLib = affyLib.2)
  
  CC.res.fisher.2 <- runTest(CCGOdata.2, algorithm = "classic", 
                             statistic = "fisher")
  allRes.2 <- GenTable(CCGOdata.2, classicFisher = CC.res.fisher.2, 
                       orderBy = "classicFisher", 
                       topNodes = top_GO_terms)
  
  # count the number of common significant ontology terms
  common.CC.GO.terms = intersect(allRes.1$GO.ID, allRes.2$GO.ID)
  
  
  #########     Output     #########
  
  cat(length(common.genes), "of the", top_genes, 
      "genes (", (length(common.genes)/top_genes)*100 ,
      "% ) are the same in the two datasets.\n")
  
  cat(length(common.BP.GO.terms), "of the", top_GO_terms, 
      "GO-BP terms (", 
      (length(common.BP.GO.terms)/top_GO_terms)*100 ,
      "% ) are the same in the two datasets.\n")
  
  
  cat(length(common.MF.GO.terms), "of the", top_GO_terms, 
      "GO-MF terms (", 
      (length(common.MF.GO.terms)/top_GO_terms)*100 ,
      "% ) are the same in the two datasets.\n")
  
  
  cat(length(common.CC.GO.terms), "of the", top_GO_terms, 
      "GO-CC terms (", 
      (length(common.CC.GO.terms)/top_GO_terms)*100 ,
      "% ) are the same in the two datasets.\n")
  
  cat("Common genes: ", length(common.features))
}