library(Matrix)
library(data.table)
library(fastmatch)
library(AUCell)


#   ____________________________________________________________________________
#   Question 2 reviewer 2                                                   ####

# functions ---------------------------------------------------------------

AUCellID <- function (X, pathways, reduction = "mca", n.features = 200, features = NULL, dims = seq(50), minSize = 10, log.trans = TRUE, p.adjust = TRUE) 
{
  DT <- getRanking(AUCell_buildRankings(X@assays$RNA@counts, nCores = 16))
  message("ranking genes")
  features <- rownames(DT)
  cells <- colnames(DT)
  i <- pbapply::pbapply(DT, 2, order)[seq(n.features), ]
  j <- rep(seq(ncol(DT)), each = n.features)
  TargetMatrix <- sparseMatrix(i, j, x = 1, dims = c(length(features), 
                                                     length(cells)), dimnames = list(features, cells))
  pathways <- lapply(pathways, function(x) x[x %fin% features])
  pathways <- pathways[sapply(pathways, function(x) length(x) >= 
                                minSize)]
  message("calculating number of success\n")
  PathwayMat <- pbapply::pbsapply(pathways, function(x) which(features %fin% 
                                                                x), simplify = F)
  PathwayLen <- unlist(lapply(PathwayMat, length))
  j <- rep(seq(length(PathwayMat)), times = PathwayLen)
  PathwayMatrix <- sparseMatrix(unlist(PathwayMat), j, x = 1, 
                                dims = c(length(features), length(PathwayMat)), dimnames = list(features, 
                                                                                                names(PathwayMat)))
  q <- as.data.frame((t(TargetMatrix) %*% PathwayMatrix) - 
                       1)
  m <- sapply(pathways, function(x) sum(x %fin% features))
  n <- sapply(m, function(x) length(features) - x)
  k <- n.features
  message("performing hypergeometric test\n")
  A <- pbapply::pbmapply(FUN = function(q, m, n, k) {
    listhyper <- phyper(seq(-1, max(q)), m, n, k, lower.tail = F)[q + 
                                                                    2]
    return(listhyper)
  }, q = q, m = m, n = n, k = k)
  rownames(A) <- rownames(q)
  A <- t(A)
  if (p.adjust) {
    A <- apply(A, 2, function(x) p.adjust(x, "BH"))
  }
  if (log.trans) {
    A <- as.sparse(-log10(A))
  }
  return(A)
}

NaiveCellID <- function (X, pathways, reduction = "mca", n.features = 200, features = NULL, dims = seq(50), minSize = 10, log.trans = TRUE, p.adjust = TRUE) 
{
  RM <- (rowMeans(X@assays$RNA@counts)+1)
  DT <- apply(X@assays$RNA@counts, 2, function(x) frank(log((RM+1)/(x+1)), ties.method = "random"))
  message("ranking genes")
  features <- rownames(X@assays$RNA@counts)
  cells <- colnames(X@assays$RNA@counts)
  i <- pbapply::pbapply(DT, 2, order)[seq(n.features), ]
  j <- rep(seq(ncol(DT)), each = n.features)
  TargetMatrix <- sparseMatrix(i, j, x = 1, dims = c(length(features), 
                                                     length(cells)), dimnames = list(features, cells))
  pathways <- lapply(pathways, function(x) x[x %fin% features])
  pathways <- pathways[sapply(pathways, function(x) length(x) >= 
                                minSize)]
  message("calculating number of success\n")
  PathwayMat <- pbapply::pbsapply(pathways, function(x) which(features %fin% 
                                                                x), simplify = F)
  PathwayLen <- unlist(lapply(PathwayMat, length))
  j <- rep(seq(length(PathwayMat)), times = PathwayLen)
  PathwayMatrix <- sparseMatrix(unlist(PathwayMat), j, x = 1, 
                                dims = c(length(features), length(PathwayMat)), dimnames = list(features, 
                                                                                                names(PathwayMat)))
  q <- as.data.frame((t(TargetMatrix) %*% PathwayMatrix) - 
                       1)
  m <- sapply(pathways, function(x) sum(x %fin% features))
  n <- sapply(m, function(x) length(features) - x)
  k <- n.features
  message("performing hypergeometric test\n")
  A <- pbapply::pbmapply(FUN = function(q, m, n, k) {
    listhyper <- phyper(seq(-1, max(q)), m, n, k, lower.tail = F)[q + 
                                                                    2]
    return(listhyper)
  }, q = q, m = m, n = n, k = k)
  rownames(A) <- rownames(q)
  A <- t(A)
  if (p.adjust) {
    A <- apply(A, 2, function(x) p.adjust(x, "BH"))
  }
  if (log.trans) {
    A <- as.sparse(-log10(A))
  }
  return(A)
}


# Cite Seq ----------------------------------------------------------------

AUCellID_Cite <- AUCellID(X = cite_seurat, reduction = "mca",n.features = 200, pathways = immune_signature, log.trans = T, minSize = 5)
AUCellID_Cite_pred <- rownames(AUCellID_Cite)[apply(AUCellID_Cite, 2, which.max)]
cite_seurat@meta.data$pred_CellID_AUC <- ifelse(2 < apply(AUCellID_Cite, 2, max), AUCellID_Cite_pred, "unassigned")

NaiveCite <- NaiveCellID(X = cite_seurat, reduction = "mca", n.features = 200, pathways = immune_signature, log.trans = T, minSize = 5)
AUcite_HGT_naive <- rownames(NaiveCite)[apply(NaiveCite, 2, which.max)]
cite_seurat@meta.data$pred_Naive <-ifelse(2 < apply(NaiveCite, 2, max), AUcite_HGT_naive, "unassigned")


# Reap Seq ----------------------------------------------------------------

AUCellID_Reap <- AUCellID(X = reap_seurat, reduction = "mca", n.features = 200, pathways = immune_signature, log.trans = T, minSize = 5)
AUCellID_Reap_pred <- rownames(AUCellID_Reap)[apply(AUCellID_Reap, 2, which.max)]
reap_seurat@meta.data$pred_CellID_AUC <- ifelse(2 < apply(AUCellID_Reap, 2, max), AUCellID_Reap_pred, "unassigned")

NaiveReap <- NaiveCellID(X = reap_seurat, reduction = "mca",n.features = 200, pathways = immune_signature, log.trans = T, minSize = 5)
AUreap_HGT_naive <- rownames(NaiveReap)[apply(NaiveReap, 2, which.max)]
reap_seurat@meta.data$pred_Naive <- ifelse(2 < apply(NaiveReap, 2, max), AUreap_HGT_naive, "unassigned")


cite_seurat@assays$ImSigAUC <- CreateAssayObject(AUCellID_Cite)
cite_seurat@assays$ImSigNaive <- CreateAssayObject(NaiveCite)

# SupFig3C AUcell ranking version -----------------------------------------
citeCol <- palCol[unique(sort(ImmuneBenchCITE$cell_type))]
names(citeCol) <- unique(sort(ImmuneBenchCITE$cell_type))

SupFig3C_AUC <-
  Seurat::DoHeatmap(
    cite_seurat,
    features = rownames(cite_HGT),
    slot = "data",
    assay = "ImSigAUC",
    raster = T,
    group.by = "cell_type1" ,
    group.colors = citeCol,
    disp.max = 10,
    size = 3.5,lines.width = 110
  ) + scale_fill_gradientn(values=c(1, .2, .1, 0), colours=c("#D73027",  "#FEDF8F", "#4E7DB8", "#4575B4"), breaks = seq(0,10,2), na.value = 'lightgrey') + theme(plot.margin = margin(0.8,0.5,0,0.5, "cm")) + theme(
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    legend.position = "bottom"
  ) 

ggsave(plot = SupFig3C_AUC, filename = "../FinalFigure/SupFig3C_review_AUC.pdf", width = 14, height = 6, dpi = 600)
ggsave(plot = SupFig3C_AUC, filename = "../FinalFigure/SupFig3C_review_AUC.png", width = 14, height = 6, dpi = 600)

# SupFig3C Naive version --------------------------------------------------
citeCol <- palCol[unique(sort(ImmuneBenchCITE$cell_type))]
names(citeCol) <- unique(sort(ImmuneBenchCITE$cell_type))


SupFig3C_Naive <-
  Seurat::DoHeatmap(
    cite_seurat,
    features = rownames(cite_HGT),
    slot = "data",
    assay = "ImSigDUM",
    raster = T,
    group.by = "cell_type1" ,
    group.colors = citeCol,
    disp.max = 10,
    size = 3.5,lines.width = 110
  ) + scale_fill_gradientn(values=c(1, .2, .1, 0), colours=c("#D73027",  "#FEDF8F", "#4E7DB8", "#4575B4"), breaks = seq(0,10,2), na.value = 'lightgrey') + theme(plot.margin = margin(0.8,0.5,0,0.5, "cm")) + theme(
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    legend.position = "bottom"
  )

ggsave(plot = SupFig3C_Naive, filename = "../FinalFigure/SupFig3C_review_Naive.pdf", width = 14, height = 6, dpi = 600)
ggsave(plot = SupFig3C_Naive, filename = "../FinalFigure/SupFig3C_review_Naive.png", width = 14, height = 6, dpi = 600)

SupFig3C_review <- ggarrange(plotlist = list(SupFig3C_AUC, SupFig3C_Naive), labels = c("AUCell Ranking", "Naive Ranking"), nrow = 2, ncol = 1, common.legend = T, legend = "right")
ggsave(plot = SupFig3C_review, filename = "../FinalFigure/SupFig3C_review.pdf", width = 14, height = 12, dpi = 600)
ggsave(plot = SupFig3C_review, filename = "../FinalFigure/SupFig3C_review.png", width = 14, height = 12, dpi = 600)

# Calculate Metrics -------------------------------------------------------

map_cite <- list(c("B-cells"), c("CD14 Monocytes"), c("CD16 Monocytes"), c("HSC","MPP","GMP","MEP","CMP"),c("CD4 T-cells"), c("CD8 T-cells"),  c("cDC","DC"), c("Erythrocytes","MEP"),c("Megakaryocytes", "Platelets"),c("NK cells"), c("pDC")) %>% set_names(unique(sort(as.vector(cite_seurat@active.ident))))
map_reap <- list(c("B-cells"), c("CD14 Monocytes"), c("CD16 Monocytes"), c("CD4 T-cells"), c("CD8 T-cells"),  c("cDC","DC"), c("Megakaryocytes", "Platelets"),c("NK cells"), c("pDC")) %>% set_names(head(unique(sort(as.vector(reap_seurat@active.ident))),-1))


ImmuneBenchCell <-
  foreach(
    Seurat = list(cite_seurat, subset(
      reap_seurat, idents = "unknown", invert = T
    )),
    mapping = list(map_cite, map_reap),
    data = c("CITE", "REAP"),
    .combine = rbind
  ) %:% foreach(
    predictions = c("pred_CellID", "pred_SCINA", "pred_AUCell", "pred_Naive", "pred_CellID_AUC"),
    .combine = rbind
  ) %do% {
    TruthDF <-
      tibble(
        Author = as.vector(Seurat$cell_type1),
        Predictions = as.vector(unlist(Seurat[[predictions]])),
        map_label = as.vector(sapply(mapping[as.vector(Seurat$cell_type1)], function(x)
          paste0(x, collapse = ", ")))
      )
    TruthDF <- TruthDF %>% filter(map_label != "removed")
    TruthDF <-
      TruthDF %>% mutate(positive = mapply(
        x = .$map_label,
        y = .$Predictions,
        FUN = function(x, y) {
          x %like% y
        }
      ))
    TruthDF <-
      TruthDF %>% mutate(final_map = ifelse(positive, Author, Predictions))
    PredDF <- sapply(unique(TruthDF$Author), function(x, truth) {
      TruthDFpos <- TruthDF %>% filter(Author == x)
      TruthDFneg <- TruthDF %>% filter(Author != x)
      TP <- sum(TruthDFpos$Author == TruthDFpos$final_map)
      FN <- sum(TruthDFpos$Author != TruthDFpos$final_map)
      FP <- sum(TruthDFneg$final_map %in% mapping[[x]])
      Recall <- TP / (TP + FN)
      Precision <- TP / (FP + TP)
      F1 <- 2 * (Recall * Precision) / (Recall + Precision)
      return(c(
        Precision = Precision,
        Recall = Recall,
        F1 = F1
      ))
    },
    truth = TruthDF) %>% as.data.frame() %>%  rownames_to_column(var = "metrics") %>%  gather("cell_type", "value",-1) %>%  mutate(methods = predictions) %>%  mutate(data = data) %>%  mutate(value = ifelse(is.na(value), 0, value))
  } %>% dplyr::select(data, methods, cell_type, everything()) %>%  mutate(methods = str_remove(methods, "pred_")) %>% mutate(methods = str_replace(methods, "CellID_G", "CellID(G)")) %>% mutate(methods = str_replace(methods, "CellID_C", "CellID(C)")) %>%  mutate(cell_type = factor(cell_type, immune_levels))

# Overall assessment
ImmuneBenchOverall <- (ImmuneBenchCell  %>%  group_by(methods, data, metrics) %>%  summarise(`standard deviation` = sd(value), value = mean(value)))[,c(1,2,3,5,4)] %>%  mutate(value = round(value, digits = 2), `standard deviation` = round(`standard deviation`, digits = 2))


# Split Table -------------------------------------------------------------

ImmuneBenchCITE <-ImmuneBenchCell %>% filter(data == "CITE") %>%  dplyr::select(data, everything()) %>%  mutate(value = round(value, digits = 2))
ImmuneBenchREAP <-ImmuneBenchCell %>% filter(data == "REAP") %>%  dplyr::select(data, everything()) %>%  mutate(value = round(value, digits = 2))
ImmuneBenchOverallCITE <- ImmuneBenchOverall %>% filter(data == "CITE") %>%  dplyr::select(data, everything())
ImmuneBenchOverallREAP <- ImmuneBenchOverall %>% filter(data == "REAP") %>%  dplyr::select(data, everything())

# write supplementary table 7 ---------------------------------------------

xlsx::write.xlsx(x = as.data.frame(ImmuneBenchCITE), sheetName = "AssessmentPerCellType_CITE", file = "../FinalTable/SupTable7_review.xlsx", append = F, row.names = F)
xlsx::write.xlsx(x = as.data.frame(ImmuneBenchOverallCITE), sheetName = "AssessmentPerCellType_CITE", file = "../FinalTable/SupTable7_review.xlsx", append = T, row.names = F)
xlsx::write.xlsx(x = as.data.frame(ImmuneBenchREAP), sheetName = "AssessmentPerCellType_CITE", file = "../FinalTable/SupTable7_review.xlsx", append = T, row.names = F)
xlsx::write.xlsx(x = as.data.frame(ImmuneBenchOverallREAP), sheetName = "OverallAssessment_REAP", file = "../FinalTable/SupTable7_review.xlsx", append = T, row.names = F)


#   ____________________________________________________________________________
#   Question 3 reviewer 2                                                   ####

predictions_ery <- lapply(simulation_data_ery, function(x) ifelse(apply(x, 2, max) >2, rownames(x)[apply(x, 2, which.max)], "unassigned"))

EryTP <- sapply(predictions_ery, function(x) {x[1] == "Erythrocytes"})
EryFP <- sapply(predictions_ery, function(x) {x[-1]== "Erythrocytes"}) %>%  colSums()
EryFN <- 1-EryTP
EryPrecision <- mean(na.omit(EryTP/(EryTP+EryFP)))
EryRecall <- mean(EryTP/(EryFN+EryTP))
EryF1 <- (2 * EryPrecision * EryRecall)/(EryPrecision + EryRecall)

SupTable8 <- do.call(rbind,list(c(pDCPrecision,pDCRecall,pDCF1),
                                c(CD34Precision,CD34Recall,CD34F1),
                                c(EryPrecision,EryRecall,EryF1))) %>%  
  set_rownames(c("pDC", "CD34+", "Erythrocytes")) %>%  
  set_colnames(c("Precision", "Recall", "F1")) %>%  
  as.data.frame() %>%  
  rownames_to_column("cell_type") %>%  
  gather("metric", 'values', 2:4) %>%  
  arrange(cell_type)

simulation_data_cd34 <- foreach(y= 1:100, .final = function(x) setNames(x,paste0("data",1:100))) %dopar% {
  set.seed(y)
  pDCs <- names(cbmc_seurat@active.ident[cbmc_seurat@active.ident == "CD34"])
  nonpDCs <- names(cbmc_seurat@active.ident[cbmc_seurat@active.ident != "CD34"])
  A <- sample(pDCs,1)
  B <- sample(nonpDCs, 999)
  cells <- c(A,B)
  sub <- SubsetData(cbmc_seurat,cells = cells)
  HGT <- RunCellHGT(sub, pathways = immune_signature, log.trans = T, p.adjust = T, minSize =5)
}

predictions_CD34 <- lapply(simulation_data_cd34, function(x) ifelse(apply(x, 2, max) >2, rownames(x)[apply(x, 2, which.max)], "unassigned"))
CD34cells <- names(cbmc_seurat@active.ident[cbmc_seurat@active.ident %in% "CD34"])
Erycells <- names(cbmc_seurat@active.ident[cbmc_seurat@active.ident %in% "Eryth"])
CD34TP <- sapply(predictions_CD34, function(x) {x[names(x) %in% CD34cells] %in% c("HSC", "CMP", "GMP", "MPP", "MEP")}) 
CD34FP <- sapply(predictions_CD34, function(x) {(x[!(names(x) %in% CD34cells)] %in% c("HSC", "CMP", "GMP", "MPP", "MEP")) & (!x[(names(x) %in% Erycells)] %in% c("MEP"))}) %>% colSums()
CD34FN <- 1 - CD34TP
CD34Precision <- mean(na.omit(CD34TP/(CD34TP+CD34FP)))
CD34Recall <- mean(CD34TP/(CD34FN+CD34TP))
CD34F1 <- (2 * CD34Precision * CD34Recall)/(CD34Precision + CD34Recall)

predictions_pdc <- lapply(simulation_data_pdc, function(x) ifelse(apply(x, 2, max) >2, rownames(x)[apply(x, 2, which.max)], "unassigned"))
pDCcells <- names(cbmc_seurat@active.ident[cbmc_seurat@active.ident %in% "pDC"])

pDCTP <- sapply(predictions_pdc, function(x) {x[names(x) %in% pDCcells] %in% c("pDC")})
pDCFP <- sapply(predictions_pdc, function(x){x[!names(x) %in% pDCcells] %in% c("pDC")}) %>%  colSums()
pDCFN <- 1-pDCTP
pDCPrecision <- mean(na.omit(pDCTP/(pDCTP+pDCFP)))
pDCRecall <- mean(pDCTP/(pDCFN+pDCTP))
pDCF1 <- (2 * pDCPrecision * pDCRecall)/(pDCPrecision + pDCRecall)


SupTable8 <- do.call(rbind,list(c(pDCPrecision,pDCRecall,pDCF1),
                                c(CD34Precision,CD34Recall,CD34F1),
                                c(EryPrecision,EryRecall,EryF1))) %>%  
  set_rownames(c("pDC", "CD34+", "Erythrocytes")) %>%  
  set_colnames(c("Precision", "Recall", "F1")) %>%  
  as.data.frame() %>%  
  rownames_to_column("cell_type") %>%  
  gather("metric", 'values', 2:4) %>%  
  arrange(cell_type) %>% mutate(values = round(values, digits = 3)) %>%  tibble()


xlsx::write.xlsx(SupTable8, file = "../../../Definitive/CellIDPaperScript/FinalTable/SupTable8_review.xlsx")



# R1Q4 --------------------------------------------------------------------
library(CellID)
library(reticulate)
library(phateR)
library(cowplot)
library(ggpubr)
library(tidyverse)
library(patchwork)

reticulate::use_virtualenv("r-reticulate")
reticulate::use_python("/home/acortal/.conda/envs/r-reticulate/bin/python")

write_rds(cite_seurat_dr, "data/cite_seurat_dr.rds")
cite_seurat_dr <- readRDS("data/cite_seurat_dr.rds")
cite_seurat_dr@active.ident
cite_seurat_dr <- RunMCA(cite_seurat_dr, features = VariableFeatures(cite_seurat_dr))
cite_seurat_dr <- RunUMAP(cite_seurat_dr, reduction = "mca", dims = 1:10, reduction.name = "UMAP_mca", key = "UMAPmca_", reduction.key = "UMAP_")
cite_seurat_dr <- RunUMAP(cite_seurat_dr, reduction = "pca", dims = 1:50, reduction.name = "UMAP_pca", key = "UMAPpca_", reduction.key = "UMAP_")
cite_seurat_dr <- RunTSNE(cite_seurat_dr, reduction = "mca", dims = 1:10, reduction.name = "TSNE_mca", num_threads = 16, key = "TSNEmca_", reduction.key = "tSNE_")
cite_seurat_dr <- RunTSNE(cite_seurat_dr, reduction = "pca", dims = 1:10, reduction.name = "TSNE_pca", num_threads = 16, key = "TSNEpca_", reduction.key = "tSNE_")
phate_pca <- phate(cite_seurat_dr@reductions$pca@cell.embeddings[,1:10])
phate_mca <- phate(cite_seurat_dr@reductions$mca@cell.embeddings[,1:10])
cite_seurat_dr@reductions$PHATE_pca <- CreateDimReducObject(phate_pca$embedding, key = "PHATE_")
cite_seurat_dr@reductions$PHATE_mca <- CreateDimReducObject(phate_mca$embedding, key = "PHATE_")

pca <- DimPlot(cite_seurat_dr, reduction = "pca") + ggtitle("pca") + theme_R1Q4
mca <- DimPlot(cite_seurat_dr, reduction = "mca") + ggtitle("mca") + theme_R1Q4
UMAP_pca <- DimPlot(cite_seurat_dr, reduction = "UMAP_pca") + ggtitle("UMAP_pca") + theme_R1Q4
UMAP_mca <- DimPlot(cite_seurat_dr, reduction = "UMAP_mca") + ggtitle("UMAP_mca") + theme_R1Q4
TSNE_pca <- DimPlot(cite_seurat_dr, reduction = "TSNE_pca") + ggtitle("TSNE_pca") + theme_R1Q4
TSNE_mca <- DimPlot(cite_seurat_dr, reduction = "TSNE_mca") + ggtitle("TSNE_mca") + theme_R1Q4
PHATE_pca <- DimPlot(cite_seurat_dr, reduction = "PHATE_pca") + ggtitle("PHATE_pca")  + theme_R1Q4
PHATE_mca <- DimPlot(cite_seurat_dr, reduction = "PHATE_mca") + ggtitle("PHATE_mca")  + theme_R1Q4
dev.off()
ggImmune <- (pca/mca)|(UMAP_pca/UMAP_mca)|(TSNE_pca/TSNE_mca)|(PHATE_pca/PHATE_mca)|ggdraw(get_legend(DimPlot(cite_seurat_dr, reduction = "pca")))

Pancreas <- readRDS("../SupNote3/data/SeuratPancreas.rds")
Baron <- Pancreas[[1]]
Baron <- RunMCA(Baron, features = VariableFeatures(Baron))
Baron <- RunPCA(Baron, features = VariableFeatures(Baron), npcs = 50)

Baron <- RunUMAP(Baron, reduction = "mca", dims = 1:10, reduction.name = "UMAP_mca", key = "UMAPmca_", reduction.key = "UMAP_")
Baron <- RunUMAP(Baron, reduction = "pca", dims = 1:50, reduction.name = "UMAP_pca", key = "UMAPpca_", reduction.key = "UMAP_")
Baron <- RunTSNE(Baron, reduction = "mca", dims = 1:10, reduction.name = "TSNE_mca", num_threads = 16, key = "TSNEmca_", reduction.key = "tSNE_")
Baron <- RunTSNE(Baron, reduction = "pca", dims = 1:10, reduction.name = "TSNE_pca", num_threads = 16, key = "TSNEpca_", reduction.key = "tSNE_")
phate_pca <- phate(Baron@reductions$pca@cell.embeddings[,1:10])
phate_mca <- phate(Baron@reductions$mca@cell.embeddings[,1:10])
Baron@reductions$PHATE_pca <- CreateDimReducObject(phate_pca$embedding, key = "PHATE_")
Baron@reductions$PHATE_mca <- CreateDimReducObject(phate_mca$embedding, key = "PHATE_")

Baron_pca <- DimPlot(Baron, reduction = "pca") + ggtitle("pca") + theme_R1Q4  
Baron_mca <- DimPlot(Baron, reduction = "mca") + ggtitle("mca") + theme_R1Q4  
Baron_UMAP_pca <- DimPlot(Baron, reduction = "UMAP_pca") + ggtitle("UMAP_pca") + theme_R1Q4 
Baron_UMAP_mca <- DimPlot(Baron, reduction = "UMAP_mca") + ggtitle("UMAP_mca") + theme_R1Q4 
Baron_TSNE_pca <- DimPlot(Baron, reduction = "TSNE_pca") + ggtitle("TSNE_pca") + theme_R1Q4 
Baron_TSNE_mca <- DimPlot(Baron, reduction = "TSNE_mca") + ggtitle("TSNE_mca") + theme_R1Q4 
Baron_PHATE_pca <- DimPlot(Baron, reduction = "PHATE_pca") + ggtitle("PHATE_pca")  + theme_R1Q4 
Baron_PHATE_mca <- DimPlot(Baron, reduction = "PHATE_mca") + ggtitle("PHATE_mca")  + theme_R1Q4 
ggPancreas <- (Baron_pca/Baron_mca)|(Baron_UMAP_pca/Baron_UMAP_mca)|(Baron_TSNE_pca/Baron_TSNE_mca)|(Baron_PHATE_pca/Baron_PHATE_mca)|ggdraw(get_legend(DimPlot(Baron, reduction = "pca")))

Epithelium <- readRDS("../SupNote3/data/SeuratEpithelial.rds")
Plass <- Epithelium[[1]]
Plass <- RunMCA(Plass, features = VariableFeatures(Plass))
Plass <- RunPCA(Plass, features = VariableFeatures(Plass), npcs = 50)

Plass <- RunUMAP(Plass, reduction = "mca", dims = 1:10, reduction.name = "UMAP_mca", key = "UMAPmca_", reduction.key = "UMAP_")
Plass <- RunUMAP(Plass, reduction = "pca", dims = 1:50, reduction.name = "UMAP_pca", key = "UMAPpca_", reduction.key = "UMAP_")
Plass <- RunTSNE(Plass, reduction = "mca", dims = 1:10, reduction.name = "TSNE_mca", num_threads = 16, key = "TSNEmca_", reduction.key = "tSNE_")
Plass <- RunTSNE(Plass, reduction = "pca", dims = 1:10, reduction.name = "TSNE_pca", num_threads = 16, key = "TSNEpca_", reduction.key = "tSNE_")
phate_pca <- phate(Plass@reductions$pca@cell.embeddings[,1:10])
phate_mca <- phate(Plass@reductions$mca@cell.embeddings[,1:10])
Plass@reductions$PHATE_pca <- CreateDimReducObject(phate_pca$embedding, key = "PHATE_")
Plass@reductions$PHATE_mca <- CreateDimReducObject(phate_mca$embedding, key = "PHATE_")

theme_R1Q4 <- theme(aspect.ratio = 1, legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
Plass_pca <- DimPlot(Plass, reduction = "pca") + ggtitle("pca") + theme_R1Q4
Plass_mca <- DimPlot(Plass, reduction = "mca") + ggtitle("mca") + theme_R1Q4
Plass_UMAP_pca <- DimPlot(Plass, reduction = "UMAP_pca") + ggtitle("UMAP_pca") + theme_R1Q4
Plass_UMAP_mca <- DimPlot(Plass, reduction = "UMAP_mca") + ggtitle("UMAP_mca") + theme_R1Q4
Plass_TSNE_pca <- DimPlot(Plass, reduction = "TSNE_pca") + ggtitle("TSNE_pca") + theme_R1Q4
Plass_TSNE_mca <- DimPlot(Plass, reduction = "TSNE_mca") + ggtitle("TSNE_mca") + theme_R1Q4
Plass_PHATE_pca <- DimPlot(Plass, reduction = "PHATE_pca") + ggtitle("PHATE_pca")  +theme_R1Q4
Plass_PHATE_mca <- DimPlot(Plass, reduction = "PHATE_mca") + ggtitle("PHATE_mca")  +theme_R1Q4

ggEpithelium <- (Plass_pca/Plass_mca)|(Plass_UMAP_pca/Plass_UMAP_mca)|(Plass_TSNE_pca/Plass_TSNE_mca)|(Plass_PHATE_pca/Plass_PHATE_mca)|ggdraw(get_legend(DimPlot(Plass, reduction = "pca")))


ggsave(plot = ggImmune, filename = "../FinalFigure/R1Q4_Immune.pdf")
ggsave(plot = ggEpithelium, filename = "../FinalFigure/R1Q4_Epithelium.pdf")
ggsave(plot = ggPancreas, filename = "../FinalFigure/R1Q4_Pancreas.pdf")

ggsave(plot = ggImmune, filename = "../FinalFigure/R1Q4_Immune.png")
ggsave(plot = ggEpithelium, filename = "../FinalFigure/R1Q4_Epithelium.png")
ggsave(plot = ggPancreas, filename = "../FinalFigure/R1Q4_Pancreas.png")
