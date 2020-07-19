set.seed(1)
library(CellID)
library(Seurat)
library(SingleCellExperiment)
library(magrittr)
library(scater)  
library(xgboost) 
library(igraph)  
library(SingleR)
library(scmap)
library(batchelor)
library(scPred)
library(CHETAH)
library(scID)
library(singleCellNet)
library(data.table)
library(r2excel)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(plotly)
library(irr)
library(pbmcapply)
library(biomaRt)
library(BiocParallel)
library(patchwork)
library(foreach)
library(fastmatch)



RunSeuratBasics <-function(x) {
  y <- FindVariableFeatures(x)
  y <- NormalizeData(y)
  y <- ScaleData(y)
  y <- RunPCA(y)
  y <- RunTSNE(y, dims = 1:30, num_threads = 16)
  y <- RunUMAP(y, dims = 1:30)
  y <- SetIdent(y, value = "cell_type1")
  return(y)
}

pred_cellidc <- function(x, ref_seurat){
  set.seed(1)
  ref_seurat <- RunMCA(ref_seurat, nmcs = 50)
  cell_gs <- GetCellGeneSet(ref_seurat,dims = 1:50, n.features = 200)
  x <- RunMCA(x,  nmcs = 50)
  HGT_cells <- RunCellHGT(x, pathways = cell_gs, n.features = 200, reduction = "mca", dims = 1:50, log.trans = T)
  cell_match <- rownames(HGT_cells)[apply(HGT_cells, 2, which.max)]
  prediction <- as.vector(ref_seurat$cell_type1[cell_match])
  prediction_significant <- ifelse(apply(HGT_cells, 2, max) >= 2,prediction, "unassigned")
  return(prediction_significant)
}

pred_cellidg <- function(x, ref_seurat){
  set.seed(1)
  ref_seurat <- RunMCA(ref_seurat, nmcs = 50)
  x <- RunMCA(x,  nmcs = 50)
  group_gs <- GetGroupGeneSet(ref_seurat,dims = 1:50, group.by = "cell_type1", n.features = 50)
  HGT_groups <- RunCellHGT(x, pathways = group_gs, n.features = 200, reduction = "mca", dims = 1:50, log.trans = T)
  group_match <- rownames(HGT_groups)[apply(HGT_groups, 2, which.max)]
  prediction_significant <- ifelse(apply(HGT_groups, 2, max) >= 2, group_match, "unassigned")
  return(prediction_significant)
}

pred_SingleR <- function(x, ref_seurat) {
  set.seed(1)
  sce <- as.SingleCellExperiment(ref_seurat)
  sce2 <- as.SingleCellExperiment(x)
  sce <- sce[intersect(rownames(sce), rownames(sce2)), ]
  sce2 <- sce2[intersect(rownames(sce), rownames(sce2)), ]
  SingleR_train <-
    trainSingleR(ref = sce, labels = sce$cell_type1)
  pred_SingleR <-
    SingleR::classifySingleR(
      sce2,
      trained = SingleR_train,
      fine.tune = TRUE,
      BPPARAM = MulticoreParam(workers = 16)
    )
  pred_SingleR$pruned.labels[is.na(pred_SingleR$pruned.labels)] <-
    "unassigned"
  pred_SingleR <- pred_SingleR$pruned.labels
  return(pred_SingleR)
}

pred_CaSTLe <- function(x, ref_seurat) {
  set.seed(1)
  BREAKS = c(-1, 0, 1, 6, Inf)
  nFeatures = 100
  source = as.SingleCellExperiment(ref_seurat)
  target = as.SingleCellExperiment(x)
  ds1 = t(exprs(source))
  ds2 = t(exprs(target))
  sourceCellTypes = as.factor(colData(source)[, "cell_type1"])
  
  # 2. Unify sets, excluding low expressed genes
  source_n_cells_counts = apply(exprs(source), 1, function(x) {
    sum(x > 0)
  })
  target_n_cells_counts = apply(exprs(target), 1, function(x) {
    sum(x > 0)
  })
  common_genes = intersect(rownames(source)[source_n_cells_counts > 5],
                           rownames(target)[target_n_cells_counts > 5])
  remove(source_n_cells_counts, target_n_cells_counts)
  ds1 = ds1[, colnames(ds1) %in% common_genes]
  ds2 = ds2[, colnames(ds2) %in% common_genes]
  ds = rbind(ds1[, common_genes], ds2[, common_genes])
  isSource = c(rep(TRUE, nrow(ds1)), rep(FALSE, nrow(ds2)))
  remove(ds1, ds2)
  
  # 3. Highest mean in both source and target
  topFeaturesAvg = colnames(ds)[order(apply(ds, 2, mean), decreasing = T)]
  
  # 4. Highest mutual information in source
  topFeaturesMi = names(sort(apply(ds[isSource, ], 2, function(x) {
    compare(cut(x, breaks = BREAKS), sourceCellTypes, method = "nmi")
  }), decreasing = T))
  
  # 5. Top n genes that appear in both mi and avg
  selectedFeatures = union(head(topFeaturesAvg, nFeatures) ,
                           head(topFeaturesMi, nFeatures))
  
  # 6. remove correlated features
  tmp = cor(as.matrix(ds[, selectedFeatures]), method = "pearson")
  tmp[!lower.tri(tmp)] = 0
  selectedFeatures = selectedFeatures[apply(tmp, 2, function(x)
    any(x < 0.9))]
  remove(tmp)
  
  # 7,8. Convert data from continous to binned dummy vars
  # break datasets to bins
  dsBins = apply(ds[, selectedFeatures], 2, cut, breaks = BREAKS)
  # use only bins with more than one value
  nUniq = apply(dsBins, 2, function(x) {
    length(unique(x))
  })
  # convert to dummy vars
  ds = model.matrix( ~ . , as.data.frame(dsBins[, nUniq > 1]))
  remove(dsBins, nUniq)
  
  # 9. Classify
  train = runif(nrow(ds[isSource, ])) < 0.8
  # slightly different setup for multiclass and binary classification
  if (length(unique(sourceCellTypes)) > 2) {
    xg = xgboost(
      data = ds[isSource, ][train,] ,
      label = as.numeric(sourceCellTypes[train]) - 1,
      objective = "multi:softmax",
      num_class = length(unique(sourceCellTypes)),
      eta = 0.7 ,
      nthread = 5,
      nround = 20,
      verbose = 0,
      gamma = 0.001,
      max_depth = 5,
      min_child_weight = 10
    )
  } else {
    xg = xgboost(
      data = ds[isSource, ][train,] ,
      label = as.numeric(sourceCellTypes[train]) - 1,
      eta = 0.7 ,
      nthread = 5,
      nround = 20,
      verbose = 0,
      gamma = 0.001,
      max_depth = 5,
      min_child_weight = 10
    )
  }
  
  # 10. Predict
  predictedClasses = predict(xg, ds[!isSource,])
  pred_CaSTLe <- levels(sourceCellTypes)[predictedClasses + 1]
  return(pred_CaSTLe)
}

pred_SCN <- function(x, ref_seurat) {
  set.seed(1)
  seuratref <-
    extractSeurat(ref_seurat, exp_slot_name = "counts")
  stTM <- rownames_to_column(seuratref$sampTab, var = "cell")
  expTMraw <- seuratref$expDat
  seuratquery <- extractSeurat(x, exp_slot_name = "counts")
  stQuery <- seuratquery$sampTab
  expQuery <-  seuratquery$expDat
  inters <- intersect(rownames(expQuery), rownames(expTMraw))
  stList <- splitCommon(sampTab = stTM,  dLevel = "cell_type1")
  stTrain <- stList[[1]]
  expTrain <- expTMraw[inters, stTrain$cell]
  expQuery <- expQuery[inters,]
  class_info <- scn_train(
    stTrain = stTrain,
    expTrain = expTrain,
    dLevel = "cell_type1",
    colName_samp = "cell"
  )
  crPBMC <-
    scn_predict(class_info[['cnProc']], expQuery, nrand = 0)
  predictions <-
    assign_cate(classRes = crPBMC,
                sampTab = stQuery,
                cThresh = 0.99)$category
  pred_SCN <- as.vector(predictions)
  pred_SCN[pred_SCN == "rand"] <- "unassigned"
  return(pred_SCN)
}

pred_MNN <- function(x, ref_seurat) {
  set.seed(1)
  matrix <- x@assays$RNA@data
  integration <-
    intersect(VariableFeatures(ref_seurat) , VariableFeatures(x))
  matrix <- matrix[integration, ]
  MNN <-
    findMutualNN(
      data1 = t(as.matrix(matrix)),
      data2 = t(as.matrix(ref_seurat@assays$RNA@data[integration, ])),
      k1 = 50,
      k2 = 50
    )
  MNN_pair <- dplyr::arrange(as.data.frame(MNN), second)
  MNN_pair$cell_type1 <-
    ref_seurat$cell_type1[MNN_pair$second]
  MNN_prediction <-
    summarise(group_by(MNN_pair, first), prediction = names(which.max(table(cell_type1))))
  MNN_prediction$cell <-
    names(x$cell_type1[MNN_prediction$first])
  MNN_prediction <-
    left_join(data.frame(cell = colnames(x)), MNN_prediction)
  MNN_prediction$prediction[is.na(MNN_prediction$prediction)] <-
    "unassigned"
  pred_MNN <- MNN_prediction$prediction
  return(pred_MNN)
}

pred_CHETAH <- function(x, ref_seurat) {
  set.seed(1)
  ref <- as.SingleCellExperiment(ref_seurat)
  sce <- as.SingleCellExperiment(x)
  CHETAHpred <-
    CHETAH::CHETAHclassifier(input = sce,
                             ref_cells = ref,
                             ref_ct = "cell_type1")
  pred_CHETAH <-
    str_replace(CHETAHpred$celltype_CHETAH, "Node.*", "unassigned")
  pred_CHETAH <-
    str_replace(pred_CHETAH, "Unassigned", 'unassigned')
  return(pred_CHETAH)
}

pred_scPred <- function(x, ref_seurat) {
  set.seed(1)
  scp <- eigenDecompose(as.matrix(ref_seurat@assays$RNA@data))
  scPred::metadata(scp) <-
    set_colnames(as.data.frame(ref_seurat$cell_type1), "cell_type1")
  scp <- getFeatureSpace(scp, pVar = "cell_type1")
  scp <- trainModel(scp, seed = 1)
  scp <-
    scPredict(scp,
              newData = as.matrix(x@assays$RNA@data),
              threshold = 0.9)
  pred_scPred <- getPredictions(scp)$predClass
  return(pred_scPred)
}

pred_scID <- function(x, ref_seurat) {
  set.seed(1)
  SCID <-
    scID::scid_multiclass(
      target_gem = as.data.frame(as.matrix(x@assays$RNA@data)),
      reference_gem = as.data.frame(as.matrix(ref_seurat@assays$RNA@data)),
      reference_clusters = ref_seurat$cell_type1,
      normalize_reference = F
    )
  pred_scID <- SCID$labels
  return(pred_scID)
}

pred_Seurat <- function(x, ref_seurat){
  set.seed(1)
  ref_seurat <- FindVariableFeatures(ref_seurat)
  x <- FindVariableFeatures(x)
  transfer_anchor <- FindTransferAnchors(reference = ref_seurat, query = x, dims = 1:30)
  prediction <- TransferData(anchorset = transfer_anchor, refdata = ref_seurat@meta.data$cell_type1,dims = 1:30)
  prediction_vec <- ifelse(prediction$prediction.score.max > 0.5, prediction$predicted.id, "unassigned")
  pred_Seurat <- prediction_vec
  return(pred_Seurat)
}


pred_scmap_cluster <- function(x, ref_seurat) {
  set.seed(1)
  ref                                   <- as.SingleCellExperiment(ref_seurat)
  rowData(ref)$feature_symbol           <- rownames(ref)
  counts(ref)                           <- as.matrix(counts(ref))
  logcounts(ref)                        <- as.matrix(logcounts(ref))
  ref                                   <- selectFeatures(ref)
  ref                                   <- indexCluster(ref)
  index_cluster <- list((ref@metadata)$scmap_cluster_index)
  sce <- as.SingleCellExperiment(x)
  rowData(sce)$feature_symbol           <- rownames(sce)
  counts(sce)                           <- as.matrix(counts(sce))
  logcounts(sce)                        <- as.matrix(logcounts(sce))
  scmap_cluster <- scmapCluster(sce, index_cluster)$combined_labs
  return(scmap_cluster)
}

# Get Index
pred_scmap_cell <- function(x, ref_seurat) {
  set.seed(1)
  ref                                   <- as.SingleCellExperiment(ref_seurat)
  rowData(ref)$feature_symbol           <- rownames(ref)
  counts(ref)                           <- as.matrix(counts(ref))
  logcounts(ref)                        <- as.matrix(logcounts(ref))
  ref                                   <- selectFeatures(ref)
  ref                                   <- indexCell(ref)
  index_cell <- list(ref@metadata$scmap_cell_index)    
  sce <- as.SingleCellExperiment(x)
  rowData(sce)$feature_symbol           <- rownames(sce)
  counts(sce)                           <- as.matrix(counts(sce))
  logcounts(sce)                        <- as.matrix(logcounts(sce))
  res <- scmapCell(sce, index_cell)
  pred_scmap_cell <- scmapCell2Cluster(scmapCell_results = res,
                                       cluster_list = list(as.vector(colData(ref)$cell_type1)))$combined_labs
  return(pred_scmap_cell)
}


pred_func <- c(pred_cellidc,pred_cellidg,pred_CaSTLe,pred_CHETAH,pred_scmap_cell,pred_scmap_cluster,pred_Seurat,pred_MNN,pred_SCN,pred_scPred,pred_scID,pred_SingleR)
pred_func_name <- c("pred_CellID_C","pred_CellID_G","pred_CaSTLe","pred_CHETAH","pred_scmap_cell","pred_scmap_cluster","pred_Seurat","pred_MNN","pred_SCN","pred_scPred","pred_scID","pred_SingleR")

convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}
