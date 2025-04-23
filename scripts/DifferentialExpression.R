# Function to get a list of RangedSummarizedExperiment and SummarizedExperiment objects based on provided paths
# Removes subjects with missing values for inputted list of features
# Get experiment out of list using syntax: experimentList[[index]]
library(edgeR)
library(SummarizedExperiment)

getExperimentByLayerAndCancer <- function(layerName, cancerName, paths=NULL){
  if (is.null(paths)){ paths <- getPathsByLayerAndCancer() }
  revisedLayerName <- layerName
  if (layerName == "RNAseq"){ revisedLayerName <- "RNAseq_filtered" }
  else if (layerName == "RNAseq_Normal") (revisedLayerName <- "RNAseq_NORMAL_filtered")
  else if (layerName == "BinaryMutation") (revisedLayerName <- "binary-mutation")
  
  layerPathList <- list() # Only contains one item but needs to be this format
  if (revisedLayerName %in% names(paths)){ layerPathList[[layerName]] <- grep(cancerName, paths[[revisedLayerName]], value=TRUE) }
  else{ return(NULL) }
  if (length(layerPathList[[1]]) == 0){ return(NULL) }
  
  experiment <- getExperimentsList(layerPathList)
  return(experiment)
}

reduceExperimentBySubjects <- function(exp1, exp2){
  return(exp1[, exp1$submitter_id %in% exp2$submitter_id])
}

reduceExperimentByFeatures <- function(exp1, exp2, assayName1, assayName2=NULL){
  if (is.null(assayName2)) { assayName2 <- assayName1 }
  assay1 <- getTransposedFrame(SummarizedExperiment::assay(exp1, assayName1))
  assay2 <- getTransposedFrame(SummarizedExperiment::assay(exp2, assayName2))
  assay1Filtered <- assay1[, names(assay1) %in% names(assay2)]
  commonFeatures <- names(assay1)[names(assay1) %in% names(assay2)]
  return(exp1[commonFeatures, ])
}

getTransposedFrame <- function(matrix){
  return(as.data.frame(t(matrix)))
}

differentialExpression <- function(experimentList, survDF, anno, label="age", subtype=TRUE){
  
  print("Preprocessing...")
  experiment <- experimentList[[1]]
  experiment <- experiment[, experiment$submitter_id %in% survDF$submitter_id]
  mmFormula <- "~0 + group"
  if (subtype) { mmFormula <- paste(mmFormula, "+ subtype")}
  
  if (label == "age"){ 
    counts <- assay(experiment, "unstranded") 
    group <- survDF$Age > 60
    group[group] <- "Older"
    group[group != "Older"] <- "Younger"
    subtype <- experiment$subtype_m_rna
    voomFormula <- "groupOlder - groupYounger"
  }
  
  else if (label == "normal"){
    # Get normal and tumor experiments to have same subjects and genes
    experimentNormal <- getExperimentByLayerAndCancer("RNAseq_Normal", "BRCA")[[1]]
    reducedRNAseq <- reduceExperimentBySubjects(experiment, experimentNormal)
    reducedRNAseqNormal <- reduceExperimentBySubjects(experimentNormal, reducedRNAseq)
    reducedRNAseq <- reduceExperimentByFeatures(reducedRNAseq, reducedRNAseqNormal, "unstranded")
    reducedRNAseqNormal <- reduceExperimentByFeatures(reducedRNAseqNormal, reducedRNAseq, "unstranded")
    
    # Mark where samples originate from
    reducedRNAseq$group <- "Tumor"
    reducedRNAseqNormal$group <- "Normal"
    tumorCounts <- assay(reducedRNAseq, "unstranded")
    colnames(tumorCounts) <- paste0(colnames(tumorCounts), "_Tumor")
    normalCounts <- assay(reducedRNAseqNormal, "unstranded")
    colnames(normalCounts) <- paste0(colnames(normalCounts), "_Normal")
    
    # Combine assays and metadata
    counts <- cbind(tumorCounts, normalCounts)
    group <- c(reducedRNAseq$group, reducedRNAseqNormal$group)
    subtype <- c(reducedRNAseq$subtype_m_rna, reducedRNAseqNormal$subtype_m_rna)
    age <- c(reducedRNAseq$Age, reducedRNAseqNormal$Age)
    mmFormula <- paste(mmFormula, "+ age")
    voomFormula <- "groupTumor - groupNormal"
  }
  
  # Normalize
  d0 <- DGEList(counts)
  d0 <- calcNormFactors(d0)
  metadata <- data.frame(sample = colnames(counts), group=group, subtype=subtype)
  if (label == "normal") {metadata$age <- age}
  
  print("Voom transforming...")
  mm <- model.matrix(as.formula(mmFormula), data = metadata)
  keep <- filterByExpr(d0, mm)
  d <- d0[keep,]
  logcpm <- cpm(d, log=TRUE)
  y <- voom(d, mm, plot = T)
  tmp <- voom(d0, mm, plot = T)
  fit <- lmFit(y, mm)
  
  print("Contrasting groups...")
  contr <- makeContrasts(voomFormula, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)
  top.table$Ensembl_ID <- rownames(top.table)
  ord <- match(top.table$Ensembl_ID, anno$Ensembl_ID)
  top.table$Gene.name <- anno$Gene_Name[ord]
  return(top.table)
}
