# Function to get a list of RangedSummarizedExperiment and SummarizedExperiment objects based on provided paths
# Removes subjects with missing values for inputted list of features
# Get experiment out of list using syntax: experimentList[[index]]
getExperimentsList <- function(layerPaths, featuresToEnsure=c("Age"), removeHPVPositive=FALSE, removeVitalStatusUnreported=TRUE) {
  
  experimentList <- list()
  if (length(layerPaths) == 0){ return(NULL) }
  
  # Add each experiment
  for (i in 1:length(layerPaths)){
    name <- names(layerPaths)[i]
    if (!file.exists(layerPaths[[name]])) { return(NULL) }
    experiment <- readRDS(layerPaths[[name]])
    
    # Remove subjects with certain missing metadata
    for (feature in featuresToEnsure){
      experiment <- experiment[, !is.na(experiment[[feature]])]
    }
    if (removeHPVPositive){ experiment <- experiment[, experiment$HPV.status=="Negative"] }
    if (removeVitalStatusUnreported){ experiment <- experiment[, experiment$vital_status != "Not Reported"]}
    experiment$Age <- as.numeric(experiment$Age)
    experiment$Survival_Time <- pmax(experiment$Survival_Time, 1e-6)
    experimentList[[name]] <- experiment
  }
  return(experimentList)
}

# Function to get a list of vectors of file paths for each cancer and layer
getPathsByLayerAndCancer <- function() {
  
  layerNames <- c("miRNA", "RNAseq_filtered", "methylation", "RPPA", "RNAseq_NORMAL_filtered", "binary-mutation", "SCNA")
  layerPartialPaths <- c("miR", "RNAseq", "methylation", "RPPA", "normal/RNAseq", "mutation", "scna")
  # filteredOrAll <- c(1,1,1,1,1,2,2)
  
  commonPath <- "/restricted/projectnb/agedisease/CBMrepositoryData/TCGA-GDC/"
  filteredPath <- "/processed_data/filtered/TCGA-*_"
  allPath <- "/processed_data/all/TCGA-*"
  
  pathsByLayerAndCancer <- list()
  for (i in 1:length(layerNames)){
    
    name <- layerNames[i]
    
    # Concatenate main path and specific layer
    filePath <- paste(commonPath, layerPartialPaths[i], sep="") 
    
    # Add filtered data to the path if possible, all if not
    if (i < 6){ filePath <- paste(filePath, filteredPath, sep="") }
    else{ filePath <- paste(filePath, allPath, sep="") }
    
    # Get each cancer type for the layer and add to the list
    filePath <- paste(filePath, name, ".rds", sep="")
    pathsByLayerAndCancer[[name]] <- Sys.glob(filePath)
  }
  
  # Return list of vectors of path names, indexed by layer names
  return(pathsByLayerAndCancer)
}

getLayerPathsForSpecificCancer <- function(paths, cancerName){
  layerPaths <- list()
  for (layer in names(paths)){
    layerPaths[[layer]] <- grep(cancerName, paths[[layer]], value=TRUE)
  }
  return(layerPaths)
}

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

getAllCancers <- function(paths){
  PATH <- file.path( Sys.getenv("PCANAGE","/restricted/projectnb/agedisease/projects/pancancer_aging_pbock")) 
  simplifiedData <- list()
  cancerPaths <- list()
  
  for (layer in c("RNAseq_filtered")){
    stringBefore <- "RNAseq/processed_data/filtered/TCGA-"
    stringAfter <- "_RNAseq_filtered."
    
    for (i in 1:length(paths[[layer]])){
      currentPath <- paths[[layer]][i]
      cancerName <- stringr::str_split_i(currentPath, stringBefore, 2)
      cancerName <- stringr::str_split_i(cancerName, stringAfter, 1)
      cancerPaths[["RNAseq"]][cancerName] <- currentPath
      if (cancerName == "CHOL" || cancerName == "DLBC") { simplifiedData[[cancerName]] <- "" }
      else { simplifiedData[[cancerName]] <- list(data=readRDS(paste(PATH,"/results/RNAseq/DGE_analysis/TCGA-",cancerName, "_filtered_DESeq2_results.rds", sep="")), metric="padj") }
    }
  }
  return(list(cancerPaths=cancerPaths, simplifiedData=simplifiedData))
}

isRelativeVarianceHigh <- function(column, threshold=0.25){
  if (!is.numeric(column)){ return(column) }
  return(abs(var(column) / mean(column)) >= threshold)
}

checkIfNameMatch <- function(name, experiment=NULL){
  if (is.null(experiment)){ return(FALSE) } 
  return(any(experiment$submitter_id==name))
}

graphAgeModel <- function(tests, title, caption, biasCorrection=FALSE) {
  if (biasCorrection) {
    predictions <- tests$predicted_corrected_age
  } else {
    predictions <- tests$predicted_age
  }
  
  DF <- data.frame(
    sampleID = rownames(tests),
    true_age = tests$Age,
    pred_age = predictions
  )
  
  g1 <- ggplot2::ggplot(DF, aes(x = true_age, y = pred_age)) +
    geom_point(color = "blue") +
    geom_smooth(method = "lm", se = FALSE, aes(color = "Regression Line")) +  # Adding linear model fit line with legend
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", aes(color = "Identity Line")) +  # Identity line with legend
    scale_color_manual(values = c("Regression Line" = "red")) +  # Define custom colors
    labs(
      title = paste("True Age vs Predicted Age (", title, ")", sep=""),
      x = "True Age",
      y = "Predicted Age",
      color = "Legend"
    ) +
    labs(caption = caption) +
    ylim(10, 100) +  # Fix y-axis range
    theme(
      plot.caption.position = "plot",
      plot.caption = element_text(hjust = 0)
    )
  
  print(g1)
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

getFold <- function(data, fold){
  meta_trn = data[-fold, ]
  meta_tst = data[fold, ]
  return(list(train = meta_trn, test = meta_tst))
}

getNormalSize <- function(cancerName){
  normExp <- getExperimentByLayerAndCancer("RNAseq_Normal", cancerName)
  tumorExp <- getExperimentByLayerAndCancer("RNAseq", cancerName)
  experimentsList <- list(RNAseq = reduceExperimentBySubjects(tumorExp[[1]], normExp[[1]]))
  needsImputation <- list(RNAseq = TRUE, RNAseq_Normal = FALSE)
  assayNames <- list(RNAseq = "DESeq2_log", RNAseq_Normal = "DESeq2_log")
  imputedExperiments <- imputeMissingValues(experimentsList, assayNames, needsImputation, cancerName)
  return(imputedExperiments)
}

getTransposedFrame <- function(matrix){
  return(as.data.frame(t(matrix)))
}

getCancerLayersPrefix <- function(cancerType, layers){
  return(paste0(cancerType, "_", paste(layers, collapse = "_")))
}

nameMapping <- function(layerName){
  if (layerName == "RNAseq"){ return("RNAseq_filtered") }
  if (layerName == "RNAseq_Normal"){ return("RNAseq_NORMAL_filtered") }
  return(layerName)
}

makeSampleCountsTable <- function(layers){
  
  paths <- getPathsByLayerAndCancer()
  cancers <- getAllCancers(paths)
  cancerNames <- names(cancers$cancerPaths$RNAseq)
  nFrame <- data.frame(cancer=cancerNames)
  
  for (layer in layers){
    n_vec <- c()
    for (name in cancerNames){
      experiment <- getExperimentByLayerAndCancer(layer, name, paths=paths)[[1]]
      if (is.null(experiment)){ n_vec <- c(n_vec, 0) }
      else{ n_vec <- c(n_vec, length(experiment@colData@rownames)) }
    }
    nFrame[[layer]] <- n_vec
  }
  return(nFrame)
}

getSampleCount <- function(cancerName, layerName, sampleCountsTable=NULL){
  if (is.null(sampleCountsTable)){
    sampleCountsTable <- read.csv("/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/results/allLayersSampleCounts.csv")
    rownames(sampleCountsTable) <- sampleCountsTable$cancer
  }
  return(sampleCountsTable[cancerName, layerName])
}

getEligibleCancers <- function(layerName, cutoff=300){
  eligibleCancers <- c()
  sampleCountsTable <- read.csv("/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/results/allLayersSampleCounts.csv")
  rownames(sampleCountsTable) <- sampleCountsTable$cancer
  for (cancer in sampleCountsTable$cancer){
    if (getSampleCount(cancer, layerName, sampleCountsTable=sampleCountsTable) > cutoff){
      eligibleCancers <- c(eligibleCancers, cancer)
    }
  }
  return(eligibleCancers)
}

getZScores <- function(data){
  return((data - mean(data)) / sd(data))
}

differentialExpression <- function(experimentList, survDF, anno, label="age", subtype=TRUE){
  library(edgeR)
  
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
